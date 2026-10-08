# CMB FR × TNG galaxy tomography

Upload **`tng_cmb_fr_photoz_specz.ipynb`** to IllustrisTNG remote JupyterLab. It is
self-contained; no companion Python file or API key is needed. The existing
`correlation.ipynb` is retained as the original toy experiment.

## Run

1. Review the settings in section 1. The default uses **TNG300-3**, target bins
   `1.5 < z < 2.0`, `2.0 < z < 2.5`, `2.5 < z < 3.0`, and a `1.5° × 1.5°` patch.
2. Set `BASE_PATH` if automatic discovery of `sims.TNG/<run>/output` fails.
3. Run all cells. The first run streams gas fields from five full snapshots and
   can involve tens of GB of I/O. It does not load an entire snapshot into RAM.
4. Inspect the generated tables and figures in `fr_photoz_outputs/`. Gas grids
   are cached separately in `cache/` beside the notebook (`detect-mag/cache/`
   locally); changing galaxy errors/density does not require rereading gas.
   Existing matching caches in `fr_photoz_outputs/cache/` are copied and reused
   automatically. Cache keys include the simulation, data path, grid resolution,
   gas selection and implementation version. The cache directory is git-ignored.

Required packages: numpy, scipy, pandas, matplotlib, h5py, astropy. CAMB and tqdm
are optional. `CMB_E_MODEL='toy'` works without CAMB. For physical primary CMB EE,
choose an installed CAMB or supply a file with `ell, cl_ee` (raw Cl in μK²).

To check the pipeline without TNG data, use `MODE='synthetic'` and a **different
`WORK_DIR`**, so synthetic outputs do not overwrite a TNG run.

## Main outputs

| File | Content |
| --- | --- |
| `tradeoff_summary.csv` | Detection S/N and jointly fitted true-z amplitude errors |
| `tradeoff_all_realizations.csv` | Individual photo-z/thinning realizations, number density, closure tests |
| `tradeoff_noise_*.png` | Redshift-error versus retained-density comparison |
| `photo_spec_comparison.png` | Direct photo-z/spec-z comparison at different galaxy densities |
| `redshift_mixing_and_spectra.png` | Bin contamination, cross spectra, amplitude covariance |
| `mixing_bias_demo.csv` | Injected amplitudes fitted with/without a mixing model |
| `representative_templates.npz` | True-z templates, guards, Gaussian covariance, transfer matrix |
| `cone_rm_layers.npz`, `cone_galaxies.csv` | Projected RM layers and the selected galaxy catalog |
| `cmb_fr_demo.npz`, `cmb_fr_bmode.png` | Initial Q/U and frequency-dependent FR B-mode |
| `cmb_rm_reconstruction_diagnostic.png` | Ideal phase recovery and noisy pixel-fit diagnostic |
| `run_configuration.json` | Parameters, units, snapshots, package versions |

## What the numbers mean

The primary statistic is a **filtered RM²–galaxy cross spectrum**, not signed
B-mode–galaxy correlation. Templates are calibrated with expected photo-z PDFs;
chance correlations from thinning are not treated as a forecast signal. Extra
galaxy assignment/thinning noise, cross-bin covariance, and the `2 RM × noise`
term in squared-RM noise are included.

This is a **conditional Gaussian sensitivity study**, not an end-to-end CMB
experiment forecast. The RM noise grid is in rad/m² per unsmoothed map pixel,
not μK-arcmin, and is not inferred from the CMB demonstration. The latter uses
pixel phase/slope recovery, not an EB quadratic estimator. A white Milky Way
residual proxy is included; a spatial Milky Way model and lens B-mode are not.

The cone is an approximate periodic-box construction over the configured z
range, using nearest full snapshots. It is not a complete lightcone to the CMB.
The nominal 1000 deg² is a covariance area extrapolation, not independent
simulated sky. FFTs use a periodic-patch approximation. Connected non-Gaussian
covariance and real-survey masks remain to be added.

The first comparison controls photo-z precision and number density using the
same stellar-mass-selected parent population. It does not reproduce a named
survey's magnitude/color selection. TNG300-3 is an exploration default;
resolution, mass completeness, gas deposition, smoothing, guard widths, and
snapshot cadence need convergence checks before interpreting a physical S/N.

## Validation performed locally

- Notebook schema and Python syntax checks.
- All-cell execution in synthetic mode and with a small TNG-format HDF5 fixture.
- Electron-number deposition and RM units against an analytic uniform slab;
  proper rotation preserves the LOS sign.
- Squared-field noise power against 1800 independent Monte Carlo draws.
- True-z template closure, injected-amplitude recovery, covariance positivity,
  noiseless multi-frequency RM recovery, and decreasing S/N with added RM noise.

The actual IllustrisTNG data run must be performed in remote JupyterLab; it has
not been executed locally.

## Theoretical fiducial figures for the paper

Upload `tng_fr_basic_figures.ipynb` beside the main notebook. It reuses only the
main `fr_photoz_outputs/run_configuration.json` for cosmology and cone geometry;
it recomputes RM and galaxies directly from full gas snapshots and group catalogs.
Set `MAIN_OUTPUT_DIR` and `BASE_PATH` if the outputs/mount are elsewhere.

This notebook is dedicated to theory: no observing beam, Gaussian smoothing,
RM noise, photo-z scatter, random thinning, or survey-area extrapolation. Both
map displays and spectra use unfiltered fields. By default, all valid subhalos
with at least one star particle are used, without a stellar-mass threshold.
This maximizes the catalog sample but does not guarantee mass completeness.

Outputs in `fr_basic_figures/` (git-ignored), as PDF/300 dpi PNG and CSV/NPZ:

- `01_gas_magnetic_projection`: full-box gas surface density, mass-weighted
  magnetic rms, and signed mass-weighted LOS magnetic field (default snap 33).
- `02_density_magnetic_distribution`: full-snapshot mass-weighted density–B distribution.
- `03_lightcone_and_nz`: cone geometry, snapshot assignment, and true galaxy n(z).
- `04_fiducial_galaxy_maps`: unsmoothed true-z maps of the selected subhalos.
- `05_fiducial_rm_maps`: unfiltered signed RM and centered RM squared.
- `06_fiducial_angular_spectra`: gg, RM²–RM², and RM²–g spectra.
- `theory_cone_galaxies.csv`, `theory_cone_rm_layers.npz`: independently rebuilt cone.

All gas cells are streamed without subsampling; by default star-forming gas is
included. `GAS_SLAB_DEPTH_CMPC=None` uses the full box for the snapshot projection.
Optional finite slab depth and galaxy thresholds are explicit physical/sample
choices, not observing-resolution settings. Gas products, per-snapshot cone projections, and completed cone data are cached
in `cache/` with parameter and input-file metadata keys. Large 3D grids are
temporary working files, not newly retained caches.

`RM_GRID_N=512` and `THEORY_MAP_N=512` specify numerical grids. A three-component
float64 RM grid occupies about 3 GiB on disk per snapshot, and is processed one
snapshot at a time and removed after projection, including on ordinary exceptions.
Each snapshot has a small completed projection checkpoint, so an interrupted run
restarts at the first incomplete snapshot. Old completed 3D grids can still be
read in place, but can be removed to reclaim storage. The initial run reads every
required full gas snapshot; completed projection caches avoid recomputation. `MAKE_GAS_FIGURES=False` skips only the
snapshot gas figures, not construction of the theoretical RM cone.

Full-cell use still employs cell-center NGP deposition, linear interpolation and
nearest-snapshot cone geometry, rather than exact Voronoi ray tracing. Grid size,
angular pixels and LOS integration steps require convergence checks. The spectra
extend to numerical grid/pixel limits; high-L modes can be affected by deposition
and aliasing, especially after squaring RM. No observational L cutoff is imposed.
Raw gg includes shot noise; a separately labeled subtraction is also shown.
Cross spectra retain their signs. Magnetic projections are mass-weighted and
are distinct from electron-weighted RM. Cone positions use original continuous
subhalo positions rather than the main notebook's saved coarse pixels.

## Coeval benchmark and the bridge to lightcone correlations

`tng_fr_coeval_to_lightcone.ipynb` starts with the conditions of Zhang & Lidz
Fig. 1 and Fig. 3: TNG300-3, full coeval boxes, 500² maps and a 0.41 h⁻¹ cMpc
mesh. It produces halo/electron overdensity, projected |Bz| and RM² maps at z=0,
then periodic real-space RM²–halo and |RM|–halo correlations at twelve epochs.
Set `BASE_PATH` to the TNG output mount if discovery fails. The coeval part
does not depend on prior notebooks.

The paper does not specify every assignment-kernel/catalog-selection detail.
The notebook records its volume-weighted NGP definitions and compares separate
ne/B gridding with direct electron-number-weighted B deposition. `HALO_CATALOG`
can be `fof` or `subhalo`; the default `auto` checks both z=0 counts against the
paper's 391,144 halos and warns if that does not uniquely identify the catalog.
This is a controlled benchmark, not a claim of exact numerical reproduction.

It then compares halo and star-bearing galaxy tracers, crops the coeval map
without changing its epoch, and compares true finite-window and periodic
correlations. The lightcone part reads the unfiltered theoretical outputs of
`tng_fr_basic_figures.ipynb` from `fr_basic_figures/` (`LIGHTCONE_DIR` is editable).
It measures real-space angular correlations by zero-padding and subtracting
the uniform random-center expectation, normalized by available pixel pairs.
It also decomposes RM² into the target layer, outside layers and their cross
term, with an exact estimator-linearity closure check. Pair counts are not
independent sample counts or error estimates.

Outputs in `fr_coeval_to_lightcone/` (git-ignored):

- `01_coeval_fig1_fields`: four full-box fields at z=0 (zero pixels use the darkest color).
- `01b_fig1_definition_controls` and `fig1_definition_diagnostics.csv`: same-catalog NGP/CIC and signed LOS mean/sum controls, using cached projections and the group catalogue without rereading gas. The sum hypothesis does not establish the paper normalization or alter RM.
- `02_coeval_fig3_correlations`: real-space coeval curves at twelve redshifts.
- `03_tracer_and_assignment_controls`: fixed-epoch tracer/deposition comparisons.
- `04_finite_window_control`: full box versus finite crop.
- `05_lightcone_real_space_correlations`: edge-corrected w(theta) in each z bin.
- `06_lightcone_layer_decomposition`: within/outside/cross-term contributions.

PDF/PNG figures, CSV correlation values and JSON provenance are saved. Complete
2D projections are cached in `cache/`; temporary 3D memory-mapped work grids
(~3.73 GiB for the default 500³ mesh) are removed after projection or ordinary exceptions. A hard kernel/server crash
can leave work files behind; stop kernels before cleaning them.
All gas chunks/cells are read on the initial run. `RUN_LIGHTCONE=False` runs the
coeval and crop stages alone. No observing beam/noise or angular smoothing is used.

## Cache budget and remote cleanup

All three notebooks have `CACHE_MAX_GIB=8`, `CACHE_FREE_RESERVE_GIB=1`, and
`CACHE_MAX_PERSISTENT_GIB=1` in their configuration cells. These are configurable
local limits, not a statement of the TNG account quota. The total budget includes
temporary grids. Writes are checked conservatively using uncompressed sizes;
over-budget writes stop without evicting existing files. Output directories and
other account files are outside this budget. Avoid simultaneous cache-producing
notebook runs: capacity checks are not a cross-process reservation system.
Full gas usage, numerical grids and float64 accumulation are unchanged. Completed
NPZ caches are written atomically. The main forecast reads legacy caches in place
instead of copying them. A missing 3D grid is rebuilt only if its corresponding
projection/cone checkpoint is also missing or incompatible.

Upload `manage_cache.py` along with the updated notebooks. Stop kernels using the
cache before deletion (an active `.tmp.npy` may be the current computation, not
an abandoned file). On the TNG server:

```bash
python ~/tng-work/detect-mag/manage_cache.py
python ~/tng-work/detect-mag/manage_cache.py --delete
```

The first command only lists files. The second deletes only known 3D intermediates
larger than 1 GiB: `theory_rm_grid_*.npy`, `coeval_work_*.tmp.npy`, and
`theory_work_*.tmp.npy`. It preserves projection/cone NPZs, catalogues and science
outputs, and does not follow symlinks. Unknown large files are listed but retained.
Deleting a grid sacrifices reuse for a new cone geometry, but not the existing
completed cone or projection results. Ordinary exceptions clean up newly created
work files; hard crashes still require manual cleanup. One default coeval work
grid needs about 3.73 GiB; one default theory work grid needs 3 GiB.
