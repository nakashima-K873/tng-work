"""Mock observing experiments and descendant-mass densities for notebook 11.

No simulation is read and no output is written on import. The observing model
is deliberately phenomenological: noisy SED masses, Gaussian redshift errors,
mass-dependent completeness, and spectroscopic sampling, not synthetic JWST
fluxes, radiative transfer, or a calibrated instrument selection function.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import asdict, dataclass
from itertools import product
from pathlib import Path
import hashlib
import json
import os
import sqlite3
import tempfile

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.spatial import cKDTree
from scipy.special import expit, logsumexp, ndtr
from sklearn.ensemble import ExtraTreesRegressor

VERSION = 'observable-environment-v1'
C_KMS = 299792.458
POS = ['x_cMpc', 'y_cMpc', 'z_cMpc']
VEL = ['vx_kms', 'vy_kms', 'vz_kms']


@dataclass(frozen=True)
class Config:
    snapshots: tuple = (8,)
    seed: int = 42
    split_seed: int = 317
    archive_logmstar_min: float = 8.0
    anchor_logmstar_min: float = 9.0
    shallow_logmstar_min: float = 9.0
    deep_logmstar_min: float = 8.5
    mass_error_dex: float = 0.2
    sfr_error_dex: float = 0.3
    completeness_mid: float = 8.6
    deep_completeness_mid: float = 8.1
    completeness_width: float = 0.25
    photoz_sigma_over_1pz: float = 0.02
    spec_sigma_kms: float = 50.0
    spec_fraction: float = 0.5
    deep_spec_fraction: float = 0.8
    velocity_halfwidth_kms: float = 1200.0
    radii_cMpc: tuple = (1., 2., 5., 10., 15.)
    los_axis: int = 2
    test_fraction: float = 0.2
    cv_folds: int = 3
    max_anchors_per_snapshot: int = 2000
    cv_cap_per_snapshot: int = 150
    test_cap_per_snapshot: int = 300
    n_trees: int = 64
    max_depth: int = 10
    leaf_grid: tuple = (10, 20, 50)
    ywidth_grid: tuple = (0.15, 0.30)
    mwidth_grid: tuple = (0.10, 0.25)
    bridge_grid_size: int = 101
    min_train: int = 40
    min_validation: int = 8
    min_hosts: int = 10
    bootstrap_draws: int = 200
    region_audit_n: int = 30
    region_radius_cMpc: float = 10.0

    def validate(self, boxsize=None):
        if not self.snapshots or len(set(self.snapshots)) != len(self.snapshots):
            raise ValueError('Use distinct snapshots')
        if self.los_axis not in (0, 1, 2):
            raise ValueError('los_axis must be 0, 1 or 2')
        radii = np.asarray(self.radii_cMpc, float)
        if not len(radii) or not np.isfinite(radii).all() or np.any(radii <= 0) or np.any(np.diff(radii) <= 0):
            raise ValueError('Radii must be finite, positive, and strictly increasing')
        if self.deep_logmstar_min < self.archive_logmstar_min:
            raise ValueError('Deep selection must not extend below the stored truth floor')
        if self.deep_logmstar_min > self.shallow_logmstar_min:
            raise ValueError('Deep mass threshold must be at least as faint as shallow')
        if self.deep_completeness_mid > self.completeness_mid:
            raise ValueError('Deep imaging must have at least shallow completeness')
        if not 0 < self.spec_fraction <= self.deep_spec_fraction <= 1:
            raise ValueError('Invalid spectroscopic fractions')
        nonnegative = [self.mass_error_dex, self.sfr_error_dex,
                       self.photoz_sigma_over_1pz, self.spec_sigma_kms]
        if not np.isfinite(nonnegative).all() or min(nonnegative) < 0:
            raise ValueError('Errors must be finite and nonnegative')
        if self.completeness_width <= 0 or self.velocity_halfwidth_kms <= 0:
            raise ValueError('Positive completeness width and velocity window required')
        if not 0 < self.test_fraction < .5 or self.cv_folds < 2:
            raise ValueError('Invalid spatial split')
        if self.n_trees < 1 or min(self.leaf_grid) < 2 or min(self.ywidth_grid) <= 0 or min(self.mwidth_grid) <= 0:
            raise ValueError('Invalid density model parameters')
        if self.bridge_grid_size < 21 or self.bootstrap_draws < 1:
            raise ValueError('Increase integration grid / bootstrap draws')
        if min(self.min_train, self.min_validation, self.min_hosts,
               self.max_anchors_per_snapshot, self.cv_cap_per_snapshot,
               self.test_cap_per_snapshot) < 1:
            raise ValueError('Sample sizes must be positive')
        if boxsize is not None and 4*radii[-1] >= boxsize*(1-self.test_fraction):
            raise ValueError('Box too small for the transverse spatial embargo')


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def signature(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def atomic_csv(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name+'.', suffix='.tmp', dir=str(path.parent))
    try:
        with os.fdopen(fd, 'w', encoding='utf-8', newline='') as stream:
            frame.to_csv(stream, index=False)
        Path(name).replace(path)
    finally:
        if Path(name).exists():
            Path(name).unlink()


def discover_base(explicit=None):
    value = explicit or os.environ.get('TNG_OUTPUT_DIR')
    if value:
        base = Path(value).expanduser().resolve()
        if not base.is_dir():
            raise FileNotFoundError(base)
        return base
    for root in [Path.cwd(), *Path.cwd().parents, Path('/home/tnguser')]:
        for name in ('TNG300-1', 'L205n2500TNG'):
            base = root/'sims.TNG'/name/'output'
            if base.is_dir():
                return base.resolve()
    raise FileNotFoundError('Set BASE_PATH or TNG_OUTPUT_DIR to TNG300-1/output')


def periodic_delta(a, b, box):
    return (np.asarray(a)-np.asarray(b)+box/2) % box-box/2


def rtag(radius):
    return ('%g' % radius).replace('.', 'p')


def validate_raw(frame, meta):
    required = set(POS+VEL+['SubhaloID', 'snapshot', 'redshift', 'logMstar',
                           'SFR', 'logMz', 'GroupID', 'Nstar'])
    if not required.issubset(frame) or frame.empty:
        raise ValueError('Incomplete raw galaxy table')
    if frame.duplicated(['snapshot', 'SubhaloID']).any():
        raise ValueError('Duplicate source galaxies')
    if frame.snapshot.nunique() != 1 or frame.redshift.nunique() != 1:
        raise ValueError('A raw table must contain exactly one snapshot')
    if not np.isfinite(frame[POS+VEL+['logMstar', 'SFR']].to_numpy(float)).all() or (frame.SFR < 0).any():
        raise ValueError('Nonfinite positions, velocities, stellar masses, or SFR')
    if (frame[POS].to_numpy() < 0).any() or (frame[POS].to_numpy() >= meta['boxsize_cMpc']).any():
        raise ValueError('Coordinates must lie inside the periodic box')


def load_tng_snapshot(base, snapshot, cfg, cache_root):
    """Read full physical tracers, with no descendant or train/test restriction."""
    import glob
    import h5py
    import illustris_python as il
    from illustris_python.snapshot import snapPath
    head = il.groupcat.loadHeader(str(base), int(snapshot))
    h = float(head['HubbleParam'])
    if not np.isclose(h, .6774, atol=1e-4) or not np.isclose(head['BoxSize'], 205000):
        raise ValueError('Expected TNG300 cosmology and geometry')
    with h5py.File(snapPath(str(base), int(snapshot)), 'r') as stream:
        attrs = stream['Header'].attrs
        ndm = int(attrs['NumPart_Total'][1]) + (int(attrs.get('NumPart_Total_HighWord', [0]*6)[1]) << 32)
    if ndm != 2500**3:
        raise ValueError('Expected TNG300-1 resolution')
    z = float(head['Redshift'])
    if not 5 <= z <= 9.01:
        raise ValueError('This notebook is configured for 5 <= z <= 9 (snapshot 6 brackets z=9)')
    hz = 100*h*np.sqrt(float(head['Omega0'])*(1+z)**3+float(head['OmegaLambda']))
    paths = sorted(glob.glob(str(Path(il.groupcat.gcPath(str(base), int(snapshot))).parent/'*.hdf5')))
    if not paths:
        raise FileNotFoundError('No group catalog chunks')
    meta = dict(version=VERSION, simulation='TNG300-1', base=str(base), snapshot=int(snapshot),
                redshift=z, h=h, H_kms_Mpc=float(hz), boxsize_cMpc=float(head['BoxSize'])/(1000*h),
                archive_logmstar_min=cfg.archive_logmstar_min,
                source_files=[[p, Path(p).stat().st_size, Path(p).stat().st_mtime_ns] for p in paths])
    directory = Path(cache_root)/('raw_'+signature(meta)[:16])
    path = directory/'galaxies.csv'
    if path.exists():
        frame = pd.read_csv(path)
        validate_raw(frame, meta)
        return frame, meta
    fields = ['SubhaloFlag', 'SubhaloMassType', 'SubhaloPos', 'SubhaloVel',
              'SubhaloSFR', 'SubhaloLenType', 'SubhaloGrNr']
    sub = il.groupcat.loadSubhalos(str(base), int(snapshot), fields=fields)
    group = il.groupcat.loadHalos(str(base), int(snapshot), fields=['Group_M_Crit200', 'GroupFirstSub'])
    mass = np.asarray(sub['SubhaloMassType'][:, 4], float)*1e10/h
    keep = (np.asarray(sub['SubhaloFlag']) == 1) & np.isfinite(mass) & (mass >= 10**cfg.archive_logmstar_min)
    ids = np.flatnonzero(keep)
    gids = np.asarray(sub['SubhaloGrNr'][ids], np.int64)
    halo_mass = np.asarray(group['Group_M_Crit200'], float)*1e10/h
    if ((gids < 0) | (gids >= len(halo_mass))).any():
        raise ValueError('Invalid high-z group ID')
    mz = halo_mass[gids]
    frame = pd.DataFrame(dict(SubhaloID=ids, snapshot=int(snapshot), redshift=z,
                              GroupID=gids, logMstar=np.log10(mass[ids]),
                              SFR=np.asarray(sub['SubhaloSFR'][ids], float),
                              Nstar=np.asarray(sub['SubhaloLenType'][ids, 4], np.int64),
                              logMz=np.log10(np.where(mz > 0, mz, np.nan)),
                              is_central=np.asarray(group['GroupFirstSub'])[gids] == ids))
    for j, name in enumerate(POS):
        frame[name] = np.asarray(sub['SubhaloPos'][ids, j], float)/(1000*h) % meta['boxsize_cMpc']
    # SubhaloVel is the physical peculiar velocity in km/s (no sqrt(a) factor).
    for j, name in enumerate(VEL):
        frame[name] = np.asarray(sub['SubhaloVel'][ids, j], float)
    validate_raw(frame, meta)
    directory.mkdir(parents=True, exist_ok=True)
    atomic_csv(frame, path)
    (directory/'manifest.json').write_text(json.dumps(meta, indent=2))
    return frame, meta


class DescendantStore:
    """Resumable endpoint tracking; no I/O exception is called 'untracked'."""
    def __init__(self, base, cache_root, legacy_csv=None):
        import glob
        import illustris_python as il
        self.il, self.base = il, str(base)
        paths = glob.glob(il.sublink.treePath(self.base, 'SubLink', '*'))
        for sn in (0, 99):
            paths.extend(glob.glob(str(Path(il.groupcat.gcPath(self.base, sn)).parent/'*.hdf5')))
        # Per-snapshot offsets are catalog dependencies of treeOffsets.
        for sn in range(100):
            p = Path(il.groupcat.offsetPath(self.base, sn))
            if p.exists():
                paths.append(str(p))
        if not paths:
            raise FileNotFoundError('Missing SubLink / group catalog files')
        self.fingerprint = dict(version=VERSION, base=self.base,
            target='z0 FoF Group_M_Crit200 of SubLink RootDescendantID',
            files=[[p, Path(p).stat().st_size, Path(p).stat().st_mtime_ns] for p in sorted(set(paths))])
        directory = Path(cache_root)/('traces_'+signature(self.fingerprint)[:16])
        directory.mkdir(parents=True, exist_ok=True)
        (directory/'manifest.json').write_text(json.dumps(self.fingerprint, indent=2))
        self.conn = sqlite3.connect(str(directory/'traces.sqlite'))
        self.conn.execute('CREATE TABLE IF NOT EXISTS traces (snap INTEGER, sid INTEGER, payload TEXT, PRIMARY KEY(snap,sid))')
        self.handles, self.memory, self.pending = OrderedDict(), {}, 0
        self.offsets = np.asarray(il.sublink.subLinkOffsets(self.base, 'SubLink'), np.int64)
        self.g0 = np.asarray(il.groupcat.loadSubhalos(self.base, 99, fields='SubhaloGrNr'), np.int64)
        h = float(il.groupcat.loadHeader(self.base, 99)['HubbleParam'])
        self.m0 = np.asarray(il.groupcat.loadHalos(self.base, 99, fields='Group_M_Crit200'), float)*1e10/h
        if legacy_csv is not None:
            # Existing 09 data are only reused after provenance and source-stat validation.
            legacy_csv = Path(legacy_csv)
            meta = json.loads(legacy_csv.with_name('manifest.json').read_text())
            if meta.get('version') != 'latent-host-v1' or Path(meta['base']).resolve() != Path(base).resolve():
                raise ValueError('Legacy catalog provenance differs from notebook 09')
            for p, size, mtime in meta['files']:
                stat = Path(p).stat()
                if (stat.st_size, stat.st_mtime_ns) != (size, mtime):
                    raise ValueError('Legacy catalog source files changed')
            old = pd.read_csv(legacy_csv)
            if old.duplicated(['snapshot', 'SubhaloID']).any():
                raise ValueError('Duplicate legacy trace keys')
            cols = ['TrackingStatus', 'SubhaloID_z0', 'TerminalSnap', 'GroupID_z0', 'logM0']
            for row in old.itertuples(index=False):
                payload = {k: getattr(row, k) for k in cols}
                self.conn.execute('INSERT OR IGNORE INTO traces VALUES (?,?,?)',
                                  (int(row.snapshot), int(row.SubhaloID), json.dumps(payload)))
            self.conn.commit()

    def trace(self, snapshot, sid):
        key = (int(snapshot), int(sid))
        if key in self.memory:
            return dict(self.memory[key])
        saved = self.conn.execute('SELECT payload FROM traces WHERE snap=? AND sid=?', key).fetchone()
        if saved:
            payload = json.loads(saved[0])
        else:
            import h5py
            row, _, tree_id = map(int, self.il.sublink.treeOffsets(self.base, *key, 'SubLink'))
            payload = dict(TrackingStatus='not_in_tree', SubhaloID_z0=-1, TerminalSnap=-1,
                           GroupID_z0=-1, logM0=np.nan)
            if row >= 0:
                chunk = int(np.searchsorted(self.offsets, row, side='right')-1)
                if chunk < 0:
                    raise ValueError('Invalid SubLink file offset')
                if chunk not in self.handles:
                    self.handles[chunk] = h5py.File(self.il.sublink.treePath(self.base, 'SubLink', chunk), 'r')
                    if len(self.handles) > 4:
                        self.handles.popitem(last=False)[1].close()
                stream = self.handles[chunk]
                local = row-int(self.offsets[chunk])
                if not 0 <= local < len(stream['SubhaloID']):
                    raise ValueError('Seed outside tree file')
                if (int(stream['SubhaloID'][local]), int(stream['SnapNum'][local]), int(stream['SubfindID'][local])) != (tree_id, *key):
                    raise ValueError('SubLink seed mismatch')
                root = int(stream['RootDescendantID'][local])
                endpoint = local+root-tree_id
                if not 0 <= endpoint < len(stream['SubhaloID']) or int(stream['SubhaloID'][endpoint]) != root:
                    raise ValueError('SubLink endpoint mismatch')
                end_snap = int(stream['SnapNum'][endpoint])
                if not key[0] <= end_snap <= 99:
                    raise ValueError('Invalid terminal snapshot')
                payload['TerminalSnap'] = end_snap
                payload['TrackingStatus'] = 'terminated_before_z0'
                if end_snap == 99:
                    s0 = int(stream['SubfindID'][endpoint])
                    if not 0 <= s0 < len(self.g0):
                        raise ValueError('Invalid z0 subhalo ID')
                    g0 = int(self.g0[s0])
                    if not 0 <= g0 < len(self.m0):
                        raise ValueError('Invalid z0 FoF ID')
                    mass = self.m0[g0]
                    payload.update(SubhaloID_z0=s0, GroupID_z0=g0,
                                   logM0=float(np.log10(mass)) if mass > 0 else np.nan,
                                   TrackingStatus='success' if np.isfinite(mass) and mass > 0 else 'invalid_host_mass')
            self.conn.execute('INSERT INTO traces VALUES (?,?,?)', (*key, json.dumps(payload)))
            self.pending += 1
            if self.pending >= 100:
                self.conn.commit()
                self.pending = 0
        self.memory[key] = payload
        return dict(payload)

    def close(self):
        self.conn.commit()
        self.conn.close()
        for stream in self.handles.values():
            stream.close()


def mock_observe(raw, meta, cfg):
    """Common random numbers across shallow/deep and photo/spec experiments."""
    validate_raw(raw, meta)
    d = raw.sort_values('SubhaloID').reset_index(drop=True).copy()
    rng = np.random.default_rng(np.random.SeedSequence([cfg.seed, int(d.snapshot.iloc[0])]))
    n = len(d)
    d['logMstar_obs'] = d.logMstar+cfg.mass_error_dex*rng.normal(size=n)
    d['logSFR_proxy'] = np.log10(d.SFR+.01)+cfg.sfr_error_dex*rng.normal(size=n)
    udet, uspec, ez = rng.random(n), rng.random(n), rng.normal(size=n)
    d['p_detect_shallow'] = expit((d.logMstar_obs-cfg.completeness_mid)/cfg.completeness_width)
    d['p_detect_deep'] = expit((d.logMstar_obs-cfg.deep_completeness_mid)/cfg.completeness_width)
    d['detected_shallow'] = udet < d.p_detect_shallow
    d['detected_deep'] = udet < d.p_detect_deep
    d['spec_shallow'] = uspec < cfg.spec_fraction
    d['spec_deep'] = uspec < cfg.deep_spec_fraction
    scale = meta['H_kms_Mpc']/(1+meta['redshift'])
    d['v_systemic_true'] = scale*d[POS[cfg.los_axis]]+d[VEL[cfg.los_axis]]
    d['v_spec'] = d.v_systemic_true+cfg.spec_sigma_kms*ez
    d['v_photo'] = d.v_systemic_true+C_KMS*cfg.photoz_sigma_over_1pz*ez
    # Every selected anchor is assumed to have separately targeted spectroscopy.
    # spec_shallow/deep refer to the environmental tracers only.
    return d


def select_anchors(observed, cfg):
    selected = observed.loc[observed.detected_shallow & observed.logMstar_obs.ge(cfg.anchor_logmstar_min)].copy()
    if len(selected) > cfg.max_anchors_per_snapshot:
        selected = selected.sample(cfg.max_anchors_per_snapshot,
                                   random_state=cfg.seed+int(observed.snapshot.iloc[0]))
    return selected.sort_values('SubhaloID').reset_index(drop=True)


def attach_labels(anchors, store):
    rows = [dict(snapshot=int(r.snapshot), SubhaloID=int(r.SubhaloID),
                 **store.trace(int(r.snapshot), int(r.SubhaloID))) for r in anchors.itertuples(index=False)]
    if not rows:
        raise ValueError('No observed anchors; revise selection without using descendant mass')
    return anchors.merge(pd.DataFrame(rows), on=['snapshot', 'SubhaloID'], validate='one_to_one')


def wrapped_window_probability(delta, sigma, halfwidth, period):
    """Gaussian-error membership in a periodic velocity interval, flat prior.

    This is a phenomenological redshift-space membership weight. It is NOT
    probability of a shared z0 descendant, or a density-prior-calibrated p(z).
    """
    delta = periodic_delta(delta, 0., period)
    if not 0 < halfwidth < period/2 or sigma < 0 or not np.isfinite(sigma):
        raise ValueError('Invalid velocity window, uncertainty, or periodic extent')
    if sigma == 0:
        return (np.abs(delta) < halfwidth).astype(float)
    kmax = int(np.ceil(8*sigma/period))+1
    value = np.zeros_like(delta, float)
    for k in range(-kmax, kmax+1):
        value += ndtr((halfwidth-delta+k*period)/sigma)-ndtr((-halfwidth-delta+k*period)/sigma)
    return np.clip(value, 0., 1.)


def environment_features(observed, anchors, meta, cfg, scenario):
    """Self excluded by ID; full periodic LOS queried for photo-z weights."""
    box, z = meta['boxsize_cMpc'], meta['redshift']
    rmax = max(cfg.radii_cMpc)
    truth = scenario == 'truth3d'
    if truth:
        use = observed.logMstar.gt(cfg.archive_logmstar_min)
    elif scenario in ('photo', 'spec'):
        use = observed.detected_shallow & observed.logMstar_obs.ge(cfg.shallow_logmstar_min)
        if scenario == 'spec':
            use &= observed.spec_shallow
    elif scenario == 'deep_spec':
        use = observed.detected_deep & observed.logMstar_obs.ge(cfg.deep_logmstar_min) & observed.spec_deep
    else:
        raise ValueError('Unknown observing scenario')
    tracers = observed.loc[use].reset_index(drop=True)
    perp = [a for a in range(3) if a != cfg.los_axis]
    xycols = [POS[a] for a in perp]
    positions = tracers[POS if truth else xycols].to_numpy(float)
    query = anchors[POS if truth else xycols].to_numpy(float)
    tree = cKDTree(positions, boxsize=box)
    nearby = tree.query_ball_point(query, rmax)
    lookup = observed.set_index('SubhaloID', verify_integrity=True)
    period = meta['H_kms_Mpc']*box/(1+z)
    if cfg.velocity_halfwidth_kms >= period/2:
        raise ValueError('Velocity window exceeds half the periodic LOS')
    rows = []
    for anchor, candidates in zip(anchors.itertuples(index=False), nearby):
        t = tracers.iloc[candidates]
        t = t.loc[t.SubhaloID.ne(anchor.SubhaloID)]
        delta = periodic_delta(t[POS if truth else xycols].to_numpy(),
                               [getattr(anchor, c) for c in POS if truth] if truth else
                               [getattr(anchor, c) for c in xycols], box)
        radius = np.linalg.norm(delta, axis=1)
        if truth:
            weight = np.ones(len(t))
        elif scenario == 'photo':
            sigma = np.hypot(C_KMS*cfg.photoz_sigma_over_1pz, cfg.spec_sigma_kms)
            weight = wrapped_window_probability(t.v_photo.to_numpy()-anchor.v_spec,
                                                sigma, cfg.velocity_halfwidth_kms, period)
        else:
            dv = periodic_delta(t.v_spec.to_numpy(), anchor.v_spec, period)
            weight = (np.abs(dv) < cfg.velocity_halfwidth_kms).astype(float)
        row = dict(snapshot=int(anchor.snapshot), SubhaloID=int(anchor.SubhaloID))
        mstar = t.logMstar.to_numpy() if truth else t.logMstar_obs.to_numpy()
        for r in cfg.radii_cMpc:
            inside = radius < r
            w = weight*inside
            tag = rtag(r)
            row['lnN_R'+tag] = float(np.log1p(w.sum()))
            row['log1pMstarSum_R'+tag] = float(np.log10(1+np.dot(w, 10**mstar)))
        within = (radius < rmax) & (weight > .5)
        row['has_neighbor'] = int(within.any())
        row['nearest_cMpc'] = float(radius[within].min()) if within.any() else float(rmax)
        anchor_mass = lookup.loc[anchor.SubhaloID, 'logMstar' if truth else 'logMstar_obs']
        row['mass_gap'] = float(anchor_mass-mstar[within].max()) if within.any() else 0.
        w = weight*(radius < rmax)
        row['mean_logSFR_proxy'] = float(np.dot(w, t.logSFR_proxy)/w.sum()) if w.sum() else 0.
        row['has_sfr_neighbors'] = int(w.sum() > 0)
        xy = periodic_delta(t[xycols].to_numpy(), [getattr(anchor, c) for c in xycols], box)
        neff = w.sum()**2/np.dot(w, w) if np.dot(w, w) else 0.
        row['has_shape'] = int(neff >= 3)
        row['ellipticity'] = 0.
        if neff >= 3:
            center = np.average(xy, weights=w, axis=0)
            cov = ((xy-center)*w[:, None]).T@(xy-center)/w.sum()
            eig = np.linalg.eigvalsh(cov)
            row['ellipticity'] = float((eig[-1]-eig[0])/eig.sum()) if eig.sum() else 0.
        rows.append(row)
    return pd.DataFrame(rows).add_prefix(scenario+'__').rename(columns={
        scenario+'__snapshot': 'snapshot', scenario+'__SubhaloID': 'SubhaloID'})


def feature_sets(cfg):
    counts = lambda s, radii: [s+'__lnN_R'+rtag(r) for r in radii]
    core = ['logMstar_obs']
    spec = core+counts('spec', cfg.radii_cMpc)
    geometry = spec+['spec__log1pMstarSum_R'+rtag(max(cfg.radii_cMpc)),
                     'spec__nearest_cMpc', 'spec__mass_gap', 'spec__has_neighbor',
                     'spec__ellipticity', 'spec__has_shape']
    population = geometry+['logSFR_proxy', 'spec__mean_logSFR_proxy', 'spec__has_sfr_neighbors']
    narrow = [r for r in cfg.radii_cMpc if r <= 5]
    if not narrow:
        narrow = [min(cfg.radii_cMpc)]
    return OrderedDict([
        ('halo_oracle', ['logMz']),
        ('halo_plus_environment', ['logMz']+counts('spec', cfg.radii_cMpc)),
        ('core', core), ('photo_environment', core+counts('photo', cfg.radii_cMpc)),
        ('spec_narrow', core+counts('spec', narrow)), ('spec_wide', spec),
        ('deep_spec', core+counts('deep_spec', cfg.radii_cMpc)),
        ('geometry', geometry), ('population_proxy', population),
        ('truth3d_reference', ['logMstar']+counts('truth3d', cfg.radii_cMpc)),
        ('two_stage_population', population),
    ])


def spatial_split(frame, meta, cfg):
    """One global test slab; transverse 2*Rmax embargo and z0-host purging.

    Since the split axis is perpendicular to the LOS, even the full-depth
    photo-z cylinders of training and test anchors cannot overlap. Tracers
    themselves are observed data, and are never filtered by future labels.
    """
    box = meta['boxsize_cMpc']
    cfg.validate(box)
    axis = next(a for a in range(3) if a != cfg.los_axis)
    origin = np.random.default_rng(cfg.split_seed).uniform(0, box)
    d = frame.copy()
    d['split_coordinate'] = (d[POS[axis]]-origin) % box
    edge, buffer = cfg.test_fraction*box, 2*max(cfg.radii_cMpc)
    test = d.split_coordinate.lt(edge)
    train = d.split_coordinate.ge(edge+buffer) & d.split_coordinate.lt(box-buffer)
    test_hosts = set(d.loc[test, 'GroupID_z0'])
    host_purge = train & d.GroupID_z0.isin(test_hosts)
    d['partition'] = np.where(test, 'test', np.where(train & ~host_purge, 'train', 'embargo'))
    d.loc[host_purge, 'partition'] = 'host_purged'
    d['cv_fold'] = -1
    span = box-buffer-(edge+buffer)
    fold = np.floor((d.split_coordinate-(edge+buffer))/span*cfg.cv_folds).astype(int)
    d.loc[d.partition.eq('train'), 'cv_fold'] = fold.clip(0, cfg.cv_folds-1)
    if not set(d.loc[d.partition.eq('train'), 'GroupID_z0']).isdisjoint(test_hosts):
        raise AssertionError('Final host leakage')
    spec = dict(axis=axis, origin_cMpc=float(origin), boxsize_cMpc=float(box),
                test_edge_cMpc=float(edge), embargo_cMpc=float(buffer),
                cv_lower_cMpc=float(edge+buffer), cv_upper_cMpc=float(box-buffer))
    return d, spec


def cv_splits(train, split_spec, cfg):
    """Validation slabs are shared by all models, with spatial/host purging."""
    width = (split_spec['cv_upper_cMpc']-split_spec['cv_lower_cMpc'])/cfg.cv_folds
    splits, audit = [], []
    for fold in range(cfg.cv_folds):
        lo = split_spec['cv_lower_cMpc']+fold*width
        hi = lo+width
        val = train.cv_fold.eq(fold)
        # Linear distance is sufficient inside the already-embargoed train slab.
        x = train.split_coordinate
        far = (x < lo-split_spec['embargo_cMpc']) | (x >= hi+split_spec['embargo_cMpc'])
        held_hosts = set(train.loc[val, 'GroupID_z0'])
        fit = ~val & far & ~train.GroupID_z0.isin(held_hosts)
        a, b = train.loc[fit].copy(), train.loc[val].copy()
        if not set(a.GroupID_z0).isdisjoint(b.GroupID_z0):
            raise AssertionError('CV final-host leakage')
        splits.append((a, b))
        for sn in sorted(train.snapshot.unique()):
            aa, bb = a.loc[a.snapshot.eq(sn)], b.loc[b.snapshot.eq(sn)]
            audit.append(dict(fold=fold, snapshot=int(sn), N_train=len(aa), N_validation=len(bb),
                              hosts_train=aa.GroupID_z0.nunique(), hosts_validation=bb.GroupID_z0.nunique(),
                              **cv_requirements(len(aa), len(bb), aa.GroupID_z0.nunique(), cfg)))
    audit = pd.DataFrame(audit)
    usable = audit.groupby('snapshot').eligible.all()
    epochs = usable.index[usable].tolist()
    audit['used_by_all_models'] = audit.snapshot.isin(epochs)
    return splits, audit, epochs


def cv_requirements(n_train, n_validation, hosts_train, cfg):
    """Sample-size gate, independent of targets, scores, and held-out test.

    In addition to the explicit floor, at least the smallest candidate must
    have enough rows for two leaves. Larger leaves are checked separately
    against the smallest training fold shared by all selected epochs.
    """
    required = max(cfg.min_train, 2*min(cfg.leaf_grid))
    failures = []
    if n_train < required:
        failures.append('N_train=%d < required=%d' % (n_train, required))
    if n_validation < cfg.min_validation:
        failures.append('N_validation=%d < required=%d' % (n_validation, cfg.min_validation))
    if hosts_train < cfg.min_hosts:
        failures.append('hosts_train=%d < required=%d' % (hosts_train, cfg.min_hosts))
    return dict(eligible=not failures, required_train=int(required),
                max_splittable_leaf=int(n_train//2), reason='; '.join(failures) if failures else 'included')


class GaussianMixture:
    """Normalized conditional density per dex in log10(M/Msun)."""
    def __init__(self, means, weights, width, hosts=None):
        self.means, self.weights = np.asarray(means, float), np.asarray(weights, float)
        self.width = float(width)
        if (self.means.ndim != 1 or self.means.shape != self.weights.shape or not len(self.means)
                or not np.isfinite(self.means).all() or not np.isfinite(self.weights).all()
                or (self.weights < 0).any() or self.weights.sum() <= 0
                or not np.isfinite(width) or width <= 0):
            raise ValueError('Invalid Gaussian mixture')
        self.weights = self.weights/self.weights.sum()
        self.hosts = np.arange(len(self.means)) if hosts is None else np.asarray(hosts)
        if len(self.hosts) != len(self.means):
            raise ValueError('Host length mismatch')

    def cdf(self, value):
        x = np.asarray(value, float)
        result = ndtr((x[..., None]-self.means)/self.width)@self.weights
        return float(result) if x.ndim == 0 else result

    def logpdf(self, value):
        with np.errstate(divide='ignore'):
            lw = np.log(self.weights)
        x = np.asarray(value, float)
        result = logsumexp(lw-.5*((x[..., None]-self.means)/self.width)**2, axis=-1)
        result -= np.log(self.width)+.5*np.log(2*np.pi)
        return float(result) if x.ndim == 0 else result

    def quantile(self, q):
        if not 0 < q < 1:
            raise ValueError('Use 0 < q < 1')
        return brentq(lambda y: self.cdf(y)-q, self.means.min()-12*self.width,
                      self.means.max()+12*self.width)

    def diagnostics(self):
        w = pd.Series(self.weights).groupby(self.hosts).sum().to_numpy()
        return dict(ESS_galaxy=float(1/np.dot(self.weights, self.weights)),
                    ESS_host=float(1/np.dot(w, w)))


class ForestDensity:
    """Leaf empirical response distributions, smoothed by Gaussian kernels.

    Every tree is fitted without bootstrap. A response gets 1/n_leaf weight
    in each query leaf, averaged across trees. This is an empirical density
    estimator, not a Bayesian posterior over tree parameters.
    """
    def __init__(self, cols, leaf, ywidth, cfg, target='logM0'):
        self.cols, self.width, self.target = list(cols), float(ywidth), target
        self.leaf = int(leaf)
        self.forest = ExtraTreesRegressor(n_estimators=cfg.n_trees, min_samples_leaf=int(leaf),
            max_depth=cfg.max_depth, max_features=1., bootstrap=False, n_jobs=1, random_state=cfg.seed)

    def fit(self, frame):
        x = frame[self.cols].to_numpy(float)
        y = frame[self.target].to_numpy(float)
        if not len(frame) or not np.isfinite(x).all() or not np.isfinite(y).all():
            raise ValueError('Training features and response must be finite')
        if len(frame) < 2*self.leaf:
            raise ValueError('N_train=%d cannot support two leaves of size %d; inspect cv_candidate_audit.csv'
                             % (len(frame), self.leaf))
        self.means, self.hosts = y, frame.GroupID_z0.to_numpy()
        self.minimum, self.maximum = x.min(axis=0), x.max(axis=0)
        self.forest.fit(x, y)
        leaves = self.forest.apply(x)
        self.lookup = []
        for j in range(leaves.shape[1]):
            column = leaves[:, j]
            order = np.argsort(column, kind='stable')
            unique, starts, sizes = np.unique(column[order], return_index=True, return_counts=True)
            self.lookup.append({int(k): order[start:start+size] for k, start, size in zip(unique, starts, sizes)})
        return self

    def weights(self, query):
        x = query[self.cols].to_numpy(float)
        if not np.isfinite(x).all():
            raise ValueError('Query has missing / infinite features')
        leaves = self.forest.apply(x)
        result = np.zeros((len(x), len(self.means)))
        for i in range(len(x)):
            for j, leaf in enumerate(leaves[i]):
                ids = self.lookup[j][int(leaf)]
                result[i, ids] += 1./(len(self.lookup)*len(ids))
        return result

    def predict(self, query):
        weights = self.weights(query)
        x = query[self.cols].to_numpy(float)
        outside = ((x < self.minimum) | (x > self.maximum)).any(axis=1)
        return [GaussianMixture(self.means, w, self.width, self.hosts) for w in weights], outside


class TwoStageDensity:
    """p(logM0|D)= integral p(logM0|logMz) p(logMz|D) dlogMz.

    Both factors and their smoothing are tuned using end-to-end validation
    NLL of logM0, on the same training folds as direct prediction.
    """
    def __init__(self, cols, leaf, ywidth, mwidth, cfg):
        self.first = ForestDensity(cols, leaf, mwidth, cfg, target='logMz')
        self.second = ForestDensity(['logMz'], leaf, ywidth, cfg)
        self.grid_size, self.mwidth = cfg.bridge_grid_size, float(mwidth)

    def fit(self, frame):
        self.first.fit(frame)
        self.second.fit(frame)
        edges = np.linspace(frame.logMz.min()-8*self.mwidth,
                            frame.logMz.max()+8*self.mwidth, self.grid_size+1)
        self.centers = (edges[:-1]+edges[1:])/2
        # Analytically integrate each first-stage Gaussian within each mass bin.
        cdf = ndtr((edges[:, None]-self.first.means[None, :])/self.mwidth)
        self.mass_bin_probability = np.diff(cdf, axis=0).T
        if not np.allclose(self.mass_bin_probability.sum(axis=1), 1., atol=1e-10):
            raise ValueError('Integration grid loses first-stage probability mass')
        self.transition = self.second.weights(pd.DataFrame({'logMz': self.centers}))
        self.bridge_outside = ((self.centers < frame.logMz.min()) | (self.centers > frame.logMz.max())).astype(float)
        return self

    def predict(self, query):
        weights = self.first.weights(query)
        mass_prob = weights@self.mass_bin_probability
        final = mass_prob@self.transition
        x = query[self.first.cols].to_numpy(float)
        outside = ((x < self.first.minimum) | (x > self.first.maximum)).any(axis=1)
        self.last_bridge_extrapolated_mass = mass_prob@self.bridge_outside
        return [GaussianMixture(self.second.means, w, self.second.width, self.second.hosts) for w in final], outside


def fit_density(name, cols, params, cfg, frame):
    if name.startswith('two_stage'):
        return TwoStageDensity(cols, **params, cfg=cfg).fit(frame)
    return ForestDensity(cols, **params, cfg=cfg).fit(frame)


def density_grid(name, cfg):
    values = product(cfg.leaf_grid, cfg.ywidth_grid, cfg.mwidth_grid if name.startswith('two_stage') else (None,))
    return [dict(leaf=int(leaf), ywidth=float(yw), **({'mwidth': float(mw)} if mw is not None else {}))
            for leaf, yw, mw in values]


def density_candidate_audit(splits, epochs, cfg):
    """One common candidate set; never skip an individual difficult fold."""
    if not epochs or not splits:
        raise ValueError('At least one eligible epoch and CV split is required')
    counts = [len(fit.loc[fit.snapshot.eq(sn)]) for fit, _ in splits for sn in epochs]
    minimum = min(counts)
    rows = []
    for name in feature_sets(cfg):
        for candidate, params in enumerate(density_grid(name, cfg)):
            eligible = minimum >= 2*params['leaf']
            rows.append(dict(model=name, candidate=candidate, minimum_cv_train=minimum,
                required_train=2*params['leaf'], eligible=eligible,
                reason='included' if eligible else 'too few rows for two leaves in the smallest CV training fold',
                **params))
    return pd.DataFrame(rows)


def tune_models(train, splits, epochs, cfg, progress=print, checkpoint=None):
    """No test access; equal fold and epoch weights, common validation caps."""
    candidate_audit = density_candidate_audit(splits, epochs, cfg)
    active = candidate_audit.loc[candidate_audit.eligible]
    for name in feature_sets(cfg):
        if not active.model.eq(name).any():
            raise ValueError('No candidate supports every CV training fold for '+name)
    allowed_keys = {(r.model, int(r.candidate), fold, int(sn))
                    for r in active.itertuples(index=False)
                    for fold in range(len(splits)) for sn in epochs}
    records, selected = [], {}
    if checkpoint is not None and Path(checkpoint).exists():
        saved = pd.read_csv(checkpoint, float_precision='round_trip')
        if saved.duplicated(['model', 'candidate', 'fold', 'snapshot']).any():
            raise ValueError('Duplicate CV checkpoint keys')
        records = saved.to_dict(orient='records')
    completed = {(r['model'], int(r['candidate']), int(r['fold']), int(r['snapshot'])) for r in records}
    if not completed.issubset(allowed_keys):
        raise ValueError('CV checkpoint contains ineligible candidates or different folds; use a new run directory')
    for row in records:
        expected = density_grid(row['model'], cfg)[int(row['candidate'])]
        if any(row.get(k) != v for k, v in expected.items()):
            raise ValueError('CV checkpoint parameters differ from the current candidate grid')
    for name, cols in feature_sets(cfg).items():
        for candidate, params in enumerate(density_grid(name, cfg)):
            if not ((active.model == name) & (active.candidate == candidate)).any():
                continue
            for fold, (all_fit, all_val) in enumerate(splits):
                for sn in epochs:
                    if (name, candidate, fold, int(sn)) in completed:
                        continue
                    fit = all_fit.loc[all_fit.snapshot.eq(sn)]
                    val = all_val.loc[all_val.snapshot.eq(sn)]
                    val = val.sample(min(len(val), cfg.cv_cap_per_snapshot), random_state=cfg.seed+fold+int(sn))
                    model = fit_density(name, cols, params, cfg, fit)
                    post, _ = model.predict(val)
                    records.append(dict(model=name, candidate=candidate, fold=fold, snapshot=int(sn),
                        N=len(val), nll=float(np.mean([-p.logpdf(y) for p, y in zip(post, val.logM0)])), **params))
                    if checkpoint is not None:
                        atomic_csv(pd.DataFrame(records), checkpoint)
        scores = pd.DataFrame([r for r in records if r['model'] == name]).groupby('candidate').nll.mean()
        best = int(scores.idxmin())
        selected[name] = dict(params=density_grid(name, cfg)[best], cv_nll=float(scores.loc[best]), features=cols)
        progress('CV complete: '+name+' '+str(selected[name]['params']))
    return selected, pd.DataFrame(records)


def score_predictions(posteriors, query, model_name, outside=None, bridge_mass=None):
    rows = []
    outside = np.zeros(len(query), bool) if outside is None else outside
    bridge_mass = np.full(len(query), np.nan) if bridge_mass is None else bridge_mass
    for p, row, out, bridge in zip(posteriors, query.itertuples(index=False), outside, bridge_mass):
        y = row.logM0
        pit = p.cdf(y)
        q16, med, q84 = [p.quantile(v) for v in (.16, .5, .84)]
        prob = 1-p.cdf(14.)
        diagnostics = p.diagnostics()
        rows.append(dict(model=model_name, snapshot=int(row.snapshot), redshift=float(row.redshift),
            SubhaloID=int(row.SubhaloID), GroupID_z0=int(row.GroupID_z0), logM0=float(y),
            median=med, q16=q16, q84=q84, W68_dex=q84-q16, pit=pit, nll=-p.logpdf(y),
            coverage68=.16 <= pit <= .84, coverage95=.025 <= pit <= .975,
            p_cluster=prob, brier=(prob-float(y >= 14))**2, outside_range=bool(out),
            support_ok=bool(not out and diagnostics['ESS_host'] >= 5),
            bridge_extrapolated_mass=float(bridge), **diagnostics))
    return pd.DataFrame(rows)


def evaluate_models(train, test, selected, epochs, cfg):
    """Fit frozen models; unsupported queries remain in the headline scores."""
    output, models = [], {}
    for sn in epochs:
        fit = train.loc[train.snapshot.eq(sn)]
        query = test.loc[test.snapshot.eq(sn)]
        if query.empty:
            raise ValueError('No test anchors at snapshot '+str(sn))
        query = query.sample(min(len(query), cfg.test_cap_per_snapshot), random_state=10000+cfg.seed+int(sn))
        for name, spec in selected.items():
            model = fit_density(name, spec['features'], spec['params'], cfg, fit)
            models[(int(sn), name)] = model
            post, outside = model.predict(query)
            bridge = model.last_bridge_extrapolated_mass if name.startswith('two_stage') else None
            output.append(score_predictions(post, query, name, outside, bridge))
        prior_width = selected['core']['params']['ywidth']
        prior = GaussianMixture(fit.logM0, np.ones(len(fit)), prior_width, fit.GroupID_z0)
        output.append(score_predictions([prior]*len(query), query, 'training_prior'))
    raw = pd.concat(output, ignore_index=True)
    summary = raw.groupby(['model', 'snapshot', 'redshift']).agg(
        N=('nll', 'size'), N_hosts=('GroupID_z0', 'nunique'), NLL=('nll', 'mean'),
        Brier=('brier', 'mean'), coverage68=('coverage68', 'mean'), coverage95=('coverage95', 'mean'),
        W68_dex=('W68_dex', 'mean'), supported_fraction=('support_ok', 'mean'),
        outside_fraction=('outside_range', 'mean'), ESS_host_median=('ESS_host', 'median'),
        bridge_extrapolated_mass=('bridge_extrapolated_mass', 'mean')).reset_index()
    return raw, summary, models


COMPARISONS = (
    ('halo_plus_environment', 'halo_oracle', 'environment_at_fixed_true_Mz'),
    ('photo_environment', 'core', 'photometric_environment'),
    ('spec_wide', 'photo_environment', 'spectroscopic_environment'),
    ('spec_wide', 'spec_narrow', 'wider_environment'),
    ('deep_spec', 'spec_wide', 'deeper_and_more_complete_spectroscopy'),
    ('geometry', 'spec_wide', 'mass_profile_and_geometry'),
    ('population_proxy', 'geometry', 'SFR_proxy_information'),
    ('population_proxy', 'two_stage_population', 'direct_vs_two_stage_same_data'),
)


def paired_differences(predictions, cfg):
    """Paired resampling of test z0 hosts; epochs receive equal weight.

    Inter-host correlations can remain inside the held-out slab, so intervals
    are descriptive host-bootstrap intervals, not cosmic-variance errors.
    """
    key = ['snapshot', 'SubhaloID', 'GroupID_z0']
    rng = np.random.default_rng(cfg.seed+908)
    rows = []
    for new, old, label in COMPARISONS:
        a = predictions.loc[predictions.model.eq(new)]
        b = predictions.loc[predictions.model.eq(old)]
        joined = a.merge(b, on=key, suffixes=('_new', '_old'), validate='one_to_one')
        if len(joined) != len(a) or len(joined) != len(b):
            raise ValueError('Models were not scored on identical test anchors')
        hosts = joined.GroupID_z0.unique()
        for metric in ('nll', 'brier', 'W68_dex'):
            tab = joined[key].copy()
            tab['delta'] = joined[metric+'_new']-joined[metric+'_old']
            estimate = tab.groupby('snapshot').delta.mean().mean()
            boot = []
            for _ in range(cfg.bootstrap_draws):
                multiplicity = pd.Series(rng.choice(hosts, len(hosts), replace=True)).value_counts()
                w = tab.GroupID_z0.map(multiplicity).fillna(0).to_numpy()
                vals = []
                for _, ids in tab.groupby('snapshot').groups.items():
                    ix = np.asarray(list(ids))
                    if w[ix].sum():
                        vals.append(np.average(tab.delta.to_numpy()[ix], weights=w[ix]))
                if len(vals) == tab.snapshot.nunique():
                    boot.append(np.mean(vals))
            if not boot:
                raise ValueError('No valid paired bootstrap replicates')
            lo, hi = np.quantile(boot, [.025, .975])
            rows.append(dict(comparison=label, new=new, reference=old, metric=metric,
                             delta=float(estimate), lo95=float(lo), hi95=float(hi),
                             N=len(tab), N_hosts=len(hosts), valid_draws=len(boot)))
    return pd.DataFrame(rows)


def calibration_table(predictions, bins=5):
    tab = predictions.copy()
    tab['probability_bin'] = np.minimum((tab.p_cluster*bins).astype(int), bins-1)
    tab['is_cluster'] = tab.logM0.ge(14).astype(float)
    return tab.groupby(['model', 'snapshot', 'probability_bin']).agg(
        N=('is_cluster', 'size'), mean_probability=('p_cluster', 'mean'),
        cluster_fraction=('is_cluster', 'mean')).reset_index()


def region_audit(observed, anchors, meta, cfg, store):
    """Fixed truth 3D reference apertures; labels never used as predictors.

    Missing endpoints stay in denominators. Upper/lower bounds describe the
    unresolved allocation. No maximum-descendant-mass label is constructed.
    """
    if cfg.region_audit_n <= 0:
        return pd.DataFrame()
    chosen = anchors.sample(min(cfg.region_audit_n, len(anchors)), random_state=cfg.seed+201)
    tree = cKDTree(observed[POS].to_numpy(), boxsize=meta['boxsize_cMpc'])
    rows = []
    for anchor in chosen.itertuples(index=False):
        ids = tree.query_ball_point([getattr(anchor, c) for c in POS], cfg.region_radius_cMpc)
        members = observed.iloc[ids].copy()
        delta = periodic_delta(members[POS].to_numpy(), [getattr(anchor, c) for c in POS], meta['boxsize_cMpc'])
        members = members.loc[np.linalg.norm(delta, axis=1) < cfg.region_radius_cMpc]
        # The reference parent is all stored physical galaxies, including anchor.
        labels = [store.trace(int(anchor.snapshot), int(s)) for s in members.SubhaloID]
        label = pd.DataFrame(labels, index=members.index)
        known = label.TrackingStatus.eq('success') & np.isfinite(label.logM0) & label.GroupID_z0.ge(0)
        anchor_label = store.trace(int(anchor.snapshot), int(anchor.SubhaloID))
        for weighting in ('number', 'stellar_mass'):
            w = np.ones(len(members)) if weighting == 'number' else 10**members.logMstar.to_numpy()
            total = w.sum()
            grouped = pd.Series(w[known.to_numpy()], index=label.loc[known, 'GroupID_z0']).groupby(level=0).sum().sort_values(ascending=False)
            missing = w[~known.to_numpy()].sum()
            main = int(grouped.index[0]) if len(grouped) else -1
            largest = float(grouped.iloc[0]) if len(grouped) else 0.
            second = float(grouped.iloc[1]) if len(grouped) > 1 else 0.
            stable = bool(len(grouped) and largest > second+missing)
            anchor_known = anchor_label['TrackingStatus'] == 'success'
            anchor_weight = float(grouped.get(anchor_label['GroupID_z0'], 0)) if anchor_known else 0.
            # If the anchor endpoint is unknown, its host may be any of the
            # known hosts, not only an entirely untracked host.
            anchor_upper = (anchor_weight+missing)/total if anchor_known else (largest+missing)/total
            rows.append(dict(snapshot=int(anchor.snapshot), SubhaloID=int(anchor.SubhaloID), weighting=weighting,
                N_members=len(members), tracking_weight_fraction=float(1-missing/total),
                known_descendant_hosts=len(grouped), dominant_known_host=main, dominant_identity_stable=stable,
                f_main_lower=largest/total, f_main_upper=(largest+missing)/total,
                f_anchor_lower=anchor_weight/total, f_anchor_upper=anchor_upper,
                anchor_is_dominant_known=bool(len(grouped) and anchor_label['TrackingStatus'] == 'success'
                                               and main == anchor_label['GroupID_z0'])))
    return pd.DataFrame(rows)


class MemoryStore:
    """Endpoint lookup for deliberately synthetic software checks only."""
    def __init__(self, labels):
        self.labels = labels
        self.fingerprint = {'mode': 'synthetic_demo', 'version': VERSION}

    def trace(self, snapshot, sid):
        return dict(self.labels[(int(snapshot), int(sid))])

    def close(self):
        pass


def synthetic_inputs(cfg, n_groups=150, per_group=25, boxsize=300.):
    """Artificial clustering and labels; never a substitute for TNG results."""
    rng = np.random.default_rng(cfg.seed+710)
    centers = rng.uniform(0, boxsize, size=(n_groups, 3))
    group_mass = rng.uniform(12., 15., size=n_groups)
    raw, metadata, labels = {}, {}, {}
    for snapshot in cfg.snapshots:
        n = n_groups*per_group
        gids = np.repeat(np.arange(n_groups), per_group)
        pos = (centers[gids]+rng.normal(0, 3, (n, 3))) % boxsize
        logmass = np.maximum(cfg.archive_logmstar_min, rng.normal(8.9+.15*(group_mass[gids]-13.5), .45))
        sfr = 10**(logmass-8.2+rng.normal(0, .3, n))
        z, hz = 8.01217294886593, 1018.
        frame = pd.DataFrame(dict(SubhaloID=np.arange(n), snapshot=int(snapshot), redshift=z,
            GroupID=np.arange(n), logMstar=logmass, SFR=sfr, Nstar=np.maximum(1, (10**(logmass-7.1)).astype(int)),
            logMz=11+.3*(logmass-9)+rng.normal(0, .2, n), is_central=np.ones(n, bool)))
        for j, c in enumerate(POS):
            frame[c] = pos[:, j]
        for c in VEL:
            frame[c] = rng.normal(0, 180, n)
        for sid, host in enumerate(gids):
            labels[(int(snapshot), sid)] = dict(TrackingStatus='success', SubhaloID_z0=int(host),
                TerminalSnap=99, GroupID_z0=int(host), logM0=float(group_mass[host]))
        raw[int(snapshot)] = frame
        metadata[int(snapshot)] = dict(simulation='SYNTHETIC SOFTWARE DEMO', snapshot=int(snapshot),
            redshift=z, H_kms_Mpc=hz, boxsize_cMpc=boxsize, h=.6774,
            archive_logmstar_min=cfg.archive_logmstar_min,
            synthetic_note='Artificial coordinates, galaxies, and descendant labels; no scientific result')
    return raw, metadata, MemoryStore(labels)


def build_experiment(raw_catalogs, metadata, cfg, store, progress=print):
    """Freeze observed anchors before looking up endpoints; retain audits."""
    frames, quality, observed_by_snap, all_anchors = [], [], {}, {}
    for sn in cfg.snapshots:
        meta = metadata[sn]
        cfg.validate(meta['boxsize_cMpc'])
        observed = mock_observe(raw_catalogs[sn], meta, cfg)
        anchors = select_anchors(observed, cfg)
        labeled = attach_labels(anchors, store)
        tracked = labeled.TrackingStatus.eq('success') & np.isfinite(labeled.logM0) & labeled.GroupID_z0.ge(0)
        valid = tracked & np.isfinite(labeled.logMz)
        # Diagnostic bins include successes AND failures; no completeness correction.
        labeled['mass_bin'] = np.floor(labeled.logMstar_obs*2)/2
        labeled['tracked'] = tracked
        labeled['valid_label'] = valid
        report = labeled.groupby('mass_bin').agg(N_selected=('SubhaloID', 'size'),
            N_tracked=('tracked', 'sum'), tracking_fraction=('tracked', 'mean'),
            N_valid=('valid_label', 'sum'), model_fraction=('valid_label', 'mean'),
            median_Nstar=('Nstar', 'median')).reset_index()
        report['snapshot'], report['redshift'] = sn, meta['redshift']
        quality.append(report)
        observed_by_snap[sn], all_anchors[sn] = observed, labeled
        for scenario in ('truth3d', 'photo', 'spec', 'deep_spec'):
            features = environment_features(observed, anchors, meta, cfg, scenario)
            labeled = labeled.merge(features, on=['snapshot', 'SubhaloID'], validate='one_to_one')
        labeled = labeled.loc[valid.to_numpy()].copy()
        if labeled.empty:
            raise ValueError('No successfully tracked common anchors')
        frames.append(labeled)
        progress('Features complete: snapshot '+str(sn)+'; common tracked anchors='+str(len(labeled)))
    data = pd.concat(frames, ignore_index=True)
    if data.duplicated(['snapshot', 'SubhaloID']).any():
        raise ValueError('Duplicate anchors')
    # Each final host must have exactly one final M200c at every source epoch.
    if (data.groupby('GroupID_z0').logM0.agg(lambda x: x.max()-x.min()) > 1e-8).any():
        raise ValueError('Inconsistent descendant mass for the same z0 FoF')
    for cols in feature_sets(cfg).values():
        if not np.isfinite(data[cols].to_numpy(float)).all():
            raise ValueError('A model input is missing or infinite')
    return data, pd.concat(quality, ignore_index=True), observed_by_snap, all_anchors


def experiment_manifest(cfg, metadata, store, mode, implementation_sha256):
    return dict(version=VERSION, implementation_sha256=implementation_sha256,
                mode=mode, config=asdict(cfg), inputs=metadata, traces=store.fingerprint,
                target='anchor descendant z0 FoF M200c; galaxy-weighted per epoch',
                direct_inputs='mock noisy stellar mass and mock environmental/population summaries',
                selection='observed shallow detection and logMstar_obs cut, before descendant lookup',
                limitations=['single-snapshot periodic mock, not a lightcone',
                    'Gaussian errors and mass-dependent selection are uncalibrated design assumptions',
                    'SFR feature is a noisy instantaneous-SFR proxy, not a measured Balmer/UV SFR',
                    'no UV, line-flux, dust, metallicity, ALMA, Lyalpha or radiative-transfer likelihood',
                    'posterior is conditional on successful tracing and finite halo mass',
                    'two-stage control is learned on the same sample; not a numerical reproduction of Witten SMHM',
                    'host-bootstrap intervals omit inter-host spatial covariance and cosmic variance'])


def analytic_self_checks():
    """Cheap deterministic checks suitable for every TNG kernel startup."""
    assert np.isclose(periodic_delta(.2, 9.8, 10.), .4)
    assert wrapped_window_probability(np.array([0., 2.]), 0., 1., 10.).tolist() == [1., 0.]
    grid = np.linspace(-5, 5, 2001)
    probability = wrapped_window_probability(grid, 2., 1., 10.)
    assert np.all((probability >= 0) & (probability <= 1))
    assert np.isclose(np.trapz(probability, grid), 2., atol=1e-6)
    mixture = GaussianMixture([13., 14.], [1., 2.], .2)
    assert np.isclose(mixture.cdf(mixture.quantile(.5)), .5, atol=1e-8)
    assert np.isclose(mixture.cdf(20)-mixture.cdf(0), 1.)
    return 'Passed: periodic boundary, velocity weights, and normalized probability / quantile checks'
