"""Scientific invariants and executable smoke checks for notebook 11."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from protohalo import observational_environment as oe


def small_config(**kwargs):
    cfg = oe.Config(n_trees=6, leaf_grid=(4,), ywidth_grid=(.2,), mwidth_grid=(.15,),
                    min_train=15, min_validation=3, min_hosts=3, max_anchors_per_snapshot=400,
                    cv_cap_per_snapshot=10, test_cap_per_snapshot=12, bootstrap_draws=12,
                    region_audit_n=3, bridge_grid_size=61)
    return replace(cfg, **kwargs)


@pytest.fixture(scope='module')
def experiment():
    cfg = small_config()
    raw, meta, store = oe.synthetic_inputs(cfg, n_groups=180, per_group=20)
    data, quality, observed, anchors = oe.build_experiment(raw, meta, cfg, store, progress=lambda _: None)
    split, spec = oe.spatial_split(data, meta[8], cfg)
    return cfg, raw, meta, store, data, quality, observed, anchors, split, spec


def test_analytic_and_wrapped_photoz():
    assert oe.analytic_self_checks().startswith('Passed')
    p = oe.wrapped_window_probability(np.array([0., 2., 4.9]), 1., .8, 10.)
    assert p[0] > p[1] > p[2]
    assert np.allclose(p, oe.wrapped_window_probability(np.array([10., 12., 14.9]), 1., .8, 10.))
    # Very broad photo-z errors tend to the fractional periodic LOS window.
    assert np.allclose(oe.wrapped_window_probability(np.array([0., 4.]), 50., 1., 10.), .2)


def test_noise_is_row_order_invariant_and_depth_nested(experiment):
    cfg, raw, meta = experiment[:3]
    a = oe.mock_observe(raw[8], meta[8], cfg)
    b = oe.mock_observe(raw[8].sample(frac=1., random_state=4), meta[8], cfg)
    pd.testing.assert_frame_equal(a, b)
    assert not (a.detected_shallow & ~a.detected_deep).any()
    assert not (a.spec_shallow & ~a.spec_deep).any()


def test_environment_self_exclusion_and_periodic_geometry():
    cfg = small_config(radii_cMpc=(1., 2.), velocity_halfwidth_kms=30.,
                       mass_error_dex=0., photoz_sigma_over_1pz=0., spec_sigma_kms=0.)
    meta = dict(boxsize_cMpc=100., redshift=8., H_kms_Mpc=1000.)
    d = pd.DataFrame(dict(SubhaloID=[0, 1, 2], snapshot=[8]*3, redshift=[8.]*3,
        logMstar=[10.]*3, SFR=[2.]*3, logMz=[11.]*3, GroupID=[0, 1, 2], Nstar=[50]*3,
        x_cMpc=[.2, 99.8, .3], y_cMpc=[20.]*3, z_cMpc=[30., 30., 45.],
        vx_kms=[0.]*3, vy_kms=[0.]*3, vz_kms=[0.]*3))
    obs = oe.mock_observe(d, meta, cfg)
    obs['spec_shallow'] = True
    obs['detected_shallow'] = True
    anchor = obs.iloc[[0]]
    truth = oe.environment_features(obs, anchor, meta, cfg, 'truth3d')
    spec = oe.environment_features(obs, anchor, meta, cfg, 'spec')
    assert np.isclose(np.expm1(truth.truth3d__lnN_R1.iloc[0]), 1.)
    assert np.isclose(np.expm1(spec.spec__lnN_R1.iloc[0]), 1.)
    assert np.isclose(spec.spec__nearest_cMpc.iloc[0], .4)
    # Moving the LOS neighbor into the velocity window adds it only to cylinder.
    obs.loc[2, 'v_spec'] = obs.loc[0, 'v_spec']
    spec2 = oe.environment_features(obs, anchor, meta, cfg, 'spec')
    assert np.isclose(np.expm1(spec2.spec__lnN_R1.iloc[0]), 2.)
    # Future labels, even if present, must have no effect on any feature.
    obs['logM0'], obs['GroupID_z0'] = [10., 20., 30.], [1, 1, 1]
    pd.testing.assert_frame_equal(spec2, oe.environment_features(obs, anchor, meta, cfg, 'spec'))


def test_full_los_is_not_truncated_in_photoz():
    cfg = small_config(radii_cMpc=(1., 2.), velocity_halfwidth_kms=100.,
                       photoz_sigma_over_1pz=.02)
    meta = dict(boxsize_cMpc=300., redshift=8., H_kms_Mpc=1000.)
    d = pd.DataFrame(dict(SubhaloID=[0, 1], snapshot=[8]*2, redshift=[8.]*2,
        logMstar=[10.]*2, SFR=[2.]*2, logMz=[11.]*2, GroupID=[0, 1], Nstar=[50]*2,
        x_cMpc=[10., 10.2], y_cMpc=[20., 20.], z_cMpc=[30., 100.],
        vx_kms=[0.]*2, vy_kms=[0.]*2, vz_kms=[0.]*2))
    obs = oe.mock_observe(d, meta, cfg)
    obs['detected_shallow'] = True
    result = oe.environment_features(obs, obs.iloc[[0]], meta, cfg, 'photo')
    assert result.photo__lnN_R1.iloc[0] > 0


def test_spatial_and_host_purging(experiment):
    cfg, _, meta, _, _, _, _, _, split, spec = experiment
    tr, te = split.loc[split.partition.eq('train')], split.loc[split.partition.eq('test')]
    assert set(tr.GroupID_z0).isdisjoint(te.GroupID_z0)
    dist = abs(oe.periodic_delta(tr.split_coordinate.to_numpy()[:, None],
                                 te.split_coordinate.to_numpy()[None, :], meta[8]['boxsize_cMpc']))
    assert dist.min() >= 2*max(cfg.radii_cMpc)-1e-8
    folds, audit, epochs = oe.cv_splits(tr, spec, cfg)
    assert epochs == [8]
    for a, b in folds:
        assert set(a.GroupID_z0).isdisjoint(b.GroupID_z0)
        assert np.min(abs(a.split_coordinate.to_numpy()[:, None]-b.split_coordinate.to_numpy()[None, :])) >= spec['embargo_cMpc']-1e-8


def test_saved_tng_fold_counts_pass_without_relaxing_spatial_split():
    # Counts from the user's TNG snapshot-8 execution, not artificial 400-anchor defaults.
    cfg = oe.Config()
    for n_train, n_validation, n_hosts in [(75, 65, 57), (53, 44, 41), (88, 54, 72)]:
        result = oe.cv_requirements(n_train, n_validation, n_hosts, cfg)
        assert result['eligible']
        assert result['required_train'] == 40
    old = oe.cv_requirements(53, 44, 41, replace(cfg, min_train=60))
    assert not old['eligible']
    assert 'N_train=53 < required=60' in old['reason']
    too_small = oe.cv_requirements(39, 44, 41, cfg)
    assert not too_small['eligible']
    # A custom large-leaf grid must not silently pass the new 40-row floor.
    assert oe.cv_requirements(53, 44, 41, replace(cfg, leaf_grid=(50,)))['required_train'] == 100


def test_sparse_tng_sized_folds_exclude_unsplittable_candidates(experiment, tmp_path):
    data = experiment[4]
    cfg = replace(oe.Config(), n_trees=4, cv_cap_per_snapshot=8,
                  ywidth_grid=(.2,), mwidth_grid=(.15,), bridge_grid_size=41)
    folds = []
    for n_train, n_val in [(75, 65), (53, 44), (88, 54)]:
        fit, val = data.iloc[:n_train].copy(), data.iloc[200:200+n_val].copy()
        fit['GroupID_z0'] = np.arange(len(fit))
        val['GroupID_z0'] = np.arange(len(val))+1000
        folds.append((fit, val))
    audit = oe.density_candidate_audit(folds, [8], cfg)
    assert set(audit.loc[audit.eligible, 'leaf']) == {10, 20}
    assert not audit.loc[audit.leaf.eq(50), 'eligible'].any()
    assert audit.minimum_cv_train.eq(53).all()
    selected, scores = oe.tune_models(data, folds, [8], cfg, progress=lambda _: None,
                                     checkpoint=tmp_path/'cv.csv')
    assert set(selected) == set(oe.feature_sets(cfg))
    assert set(scores.leaf) == {10, 20}
    assert set(scores.fold) == {0, 1, 2}
    assert np.isfinite(scores.nll).all()
    # Guard against running a stale checkpoint with a different candidate order.
    with pytest.raises(ValueError, match='checkpoint'):
        oe.tune_models(data, folds, [8], replace(cfg, leaf_grid=(20, 10, 50)), progress=lambda _: None,
                       checkpoint=tmp_path/'cv.csv')


def test_density_rejects_a_leaf_larger_than_the_available_split(experiment):
    data = experiment[4]
    cfg = small_config()
    with pytest.raises(ValueError, match='cannot support two leaves'):
        oe.ForestDensity(['logMstar_obs'], 50, .2, cfg).fit(data.iloc[:53])


def test_density_and_two_stage_preserve_distributions(experiment):
    cfg, _, _, _, data = experiment[:5]
    frame = data.iloc[:180].copy()
    query = data.iloc[180:185]
    direct = oe.ForestDensity(['logMstar_obs', 'spec__lnN_R10'], 4, .2, cfg).fit(frame)
    weights = direct.weights(query)
    assert np.allclose(weights.sum(axis=1), 1)
    post, _ = direct.predict(query)
    for p in post:
        assert np.isclose(p.cdf(30)-p.cdf(0), 1)
        assert p.quantile(.16) < p.quantile(.5) < p.quantile(.84)
    two = oe.TwoStageDensity(['logMstar_obs'], 4, .2, .15, cfg).fit(frame)
    assert np.allclose(two.mass_bin_probability.sum(axis=1), 1)
    second, _ = two.predict(query)
    for p in second:
        assert np.isclose(p.weights.sum(), 1)
        assert np.count_nonzero(p.weights > 0) > 1  # not a plug-in halo-mass point
    assert np.all((two.last_bridge_extrapolated_mass >= 0) & (two.last_bridge_extrapolated_mass <= 1))


def test_region_missing_endpoints_stay_in_denominator(experiment):
    cfg, _, meta, store, _, _, observed, anchors = experiment[:8]
    class IncompleteStore:
        def trace(self, snap, sid):
            if sid % 2:
                return dict(TrackingStatus='not_in_tree', GroupID_z0=-1, logM0=np.nan)
            return store.trace(snap, sid)
    result = oe.region_audit(observed[8], anchors[8], meta[8], cfg, IncompleteStore())
    assert (result.tracking_weight_fraction < 1).any()
    assert (result.f_main_lower <= result.f_main_upper).all()
    assert np.allclose(result.f_main_upper-result.f_main_lower, 1-result.tracking_weight_fraction)
    unknown = result.SubhaloID.mod(2).eq(1)
    assert not result.loc[unknown, 'anchor_is_dominant_known'].any()
    assert np.allclose(result.loc[unknown, 'f_anchor_upper'], result.loc[unknown, 'f_main_upper'])


def test_end_to_end_common_queries_and_cv_resume(experiment, tmp_path):
    cfg, _, _, _, _, _, _, _, split, spec = experiment
    train, test = split.loc[split.partition.eq('train')], split.loc[split.partition.eq('test')]
    folds, audit, epochs = oe.cv_splits(train, spec, cfg)
    checkpoint = tmp_path/'cv.csv'
    selected, cv = oe.tune_models(train, folds, epochs, cfg, progress=lambda _: None, checkpoint=checkpoint)
    again, resumed = oe.tune_models(train, folds, epochs, cfg, progress=lambda _: None, checkpoint=checkpoint)
    assert selected == again
    assert len(cv) == len(resumed)
    pred, summary, models = oe.evaluate_models(train, test, selected, epochs, cfg)
    assert np.isfinite(pred[['nll', 'brier', 'W68_dex']]).all().all()
    keys = [set(zip(g.snapshot, g.SubhaloID)) for _, g in pred.groupby('model')]
    assert all(k == keys[0] for k in keys)
    assert set(summary.model) == set(selected) | {'training_prior'}
    differences = oe.paired_differences(pred, cfg)
    assert differences.comparison.nunique() == len(oe.COMPARISONS)
    assert (differences.lo95 <= differences.hi95).all()
    assert len(oe.calibration_table(pred)) > 0


def test_sublink_root_and_checkpoint(tmp_path, monkeypatch):
    import h5py
    tree = tmp_path/'tree.0.hdf5'
    with h5py.File(tree, 'w') as stream:
        for name, values in dict(SubhaloID=[100, 101], SnapNum=[8, 99], SubfindID=[0, 0],
                                 RootDescendantID=[101, 101]).items():
            stream[name] = values
    for sn in (0, 99):
        directory = tmp_path/str(sn)
        directory.mkdir()
        (directory/'group.0.hdf5').write_bytes(b'fake fingerprint only')
    fake = SimpleNamespace(
        sublink=SimpleNamespace(treePath=lambda base, name, chunk: str(tmp_path/('tree.'+str(chunk)+'.hdf5')),
                               subLinkOffsets=lambda base, name: [0],
                               treeOffsets=lambda base, snap, sid, name: (0, 2, 100)),
        groupcat=SimpleNamespace(gcPath=lambda base, sn: str(tmp_path/str(sn)/'group.0.hdf5'),
            offsetPath=lambda base, sn: str(tmp_path/('offset'+str(sn))),
            loadSubhalos=lambda base, sn, fields: np.array([0]),
            loadHeader=lambda base, sn: {'HubbleParam': .6774},
            loadHalos=lambda base, sn, fields: np.array([1e4])))
    monkeypatch.setitem(sys.modules, 'illustris_python', fake)
    store = oe.DescendantStore(tmp_path, tmp_path/'cache')
    result = store.trace(8, 0)
    assert result['TrackingStatus'] == 'success'
    assert np.isclose(result['logM0'], np.log10(1e14/.6774))
    store.close()
    # Prove the checkpoint is reusable without invoking a tree read.
    fake.sublink.treeOffsets = lambda *args: (_ for _ in ()).throw(AssertionError('Unexpected tree read'))
    store = oe.DescendantStore(tmp_path, tmp_path/'cache')
    assert store.trace(8, 0) == result
    store.close()


def test_standalone_notebook_matches_tested_implementation():
    import ast
    import hashlib
    import json
    root = Path(__file__).resolve().parents[1]
    source = (root/'observational_environment.py').read_text()
    nb = json.loads((root/'notebooks'/'11_tng300_observational_environment_descendant_prediction.ipynb').read_text())
    assert nb['metadata']['implementation_sha256'] == hashlib.sha256(source.encode()).hexdigest()
    code = '\n\n'.join(''.join(c['source']) for c in nb['cells'] if c['cell_type'] == 'code')
    assert "MODE = 'tng'" in code
    assert 'from protohalo' not in code
    original = {n.name: ast.dump(n) for n in ast.parse(source).body if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
    embedded = {n.name: ast.dump(n) for n in ast.parse(code).body if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
    for name, tree in original.items():
        assert tree == embedded[name], name
    assert all(c.get('execution_count') is None and not c.get('outputs') for c in nb['cells'] if c['cell_type'] == 'code')
