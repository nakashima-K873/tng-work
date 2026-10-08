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

## Fiducial figures for the paper

Upload `tng_fr_basic_figures.ipynb` separately, beside the main notebook. Run the
main notebook first and keep `fr_photoz_outputs/`; the figure notebook loads its
RM layers, galaxy catalog and run configuration without rerunning the forecast.
Set `MAIN_OUTPUT_DIR` if those outputs are elsewhere.

The new notebook writes PDF and 300 dpi PNG figures plus CSV/NPZ data to
`fr_basic_figures/` (git-ignored):

- `01_gas_magnetic_slab`: gas surface density, mass-weighted magnetic rms, and
  signed mass-weighted LOS magnetic field in a periodic slab (default snap 33).
- `02_density_magnetic_distribution`: full-snapshot, mass-weighted density–B distribution.
- `03_lightcone_and_nz`: cone geometry, snapshot assignment, and true galaxy n(z).
- `04_fiducial_galaxy_maps`: full-density, true-z galaxy maps.
- `05_fiducial_rm_maps`: signed filtered RM and centered RM squared.
- `06_fiducial_angular_spectra`: gg, RM²–RM², and RM²–g spectra.

Gas products require one full gas snapshot, streamed in bounded batches, and
are cached in `cache/` using settings and input-file metadata. Existing RM grids
cannot reconstruct gas density or the magnitude of B, so a separate initial
gas read is necessary. Change `BASE_PATH` if the simulation mount moved.
`MAKE_GAS_FIGURES=False` produces the remaining figures from main outputs alone.

Galaxy spectra use unsmoothed count maps; only the displayed galaxy maps are
smoothed. The raw gg spectra include shot noise; a separately labeled subtraction
is also shown. All spectra are finite-patch measurements without observational
noise or survey-area extrapolation. Cross spectra retain their signs. Gas magnetic
projections are mass-weighted and are distinct from electron-weighted RM.
The cone plot uses saved angular pixel centers, not exact continuous galaxy positions.
