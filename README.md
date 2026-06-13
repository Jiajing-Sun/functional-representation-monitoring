# R Code for Online Monitoring of Functional Representations

This repository contains the public R code for the associated manuscript
`When to Refresh a Functional Representation: Online Monitoring of Structural Change in Functional Time Series`.

The repository contains code only. It excludes raw TickData files, generated CSV/RDS outputs, figures, LaTeX tables,
cached Monte Carlo results, and old source snapshots.

## Repository Layout

- `critical_values/`: finite-horizon critical-value simulation code for the main KS/weighted-CvM monitors and the supplementary detector families.
- `simulations/`: Monte Carlo designs, contaminated-training designs, HAC kernel-and-bandwidth sensitivity checks, and diagnostics.
- `empirical_analysis/`: S&P 500 intraday empirical scripts. Licensed TickData input is not included.
- `tests/`: parse checks, private-artifact checks, and a small critical-value smoke test.

## Quick Start

1. Check R package dependencies:

```sh
Rscript requirements.R --check
```

2. Run lightweight repository checks:

```sh
./run_smoke_tests.sh
```

3. Generate finite-horizon main critical values:

```sh
Rscript critical_values/main/run_generate_streamingcurve_cv.R --ncores=4
```

4. Generate finite-horizon supplementary detector critical values when needed:

```sh
Rscript critical_values/alternative_detectors/run_generate_alt_detector_cv.R --ncores=4
```

5. Run simulations after the relevant critical-value tables exist. For example:

```sh
Rscript simulations/main/run_null_streamingcurve_simulation.R --ncores=4
Rscript simulations/main/run_power_streamingcurve_simulation.R --scenario=level_shift --ncores=4
```

6. Run empirical scripts only after placing the licensed TickData file at `empirical_analysis/data/SPX.csv`, or pass an explicit file:

```sh
Rscript empirical_analysis/run_empirical_bundle_refresh.R --data_file=/path/to/SPX.csv
```

## Monitoring Horizons

The main simulations and empirical application use finite monitoring horizons. The public critical-value generators
therefore default to `T = 1, 2, 5, 10`. The appendix also reports open-end KS and weighted-CvM critical values; those
simulations can be run explicitly with `--include_openend=true`.

## Outputs

Generated outputs are ignored by `.gitignore`. The canonical public critical-value locations are:

- `critical_values/main/outputs/critical_values_all.csv`
- `critical_values/main/outputs/critical_values_all_weights.csv`
- `critical_values/alternative_detectors/outputs/critical_values_alt_detectors.csv`

Simulation and empirical scripts look for those paths before falling back to local output folders.

## Data

The empirical application uses licensed 1-minute S&P 500 data from TickData. The data are not redistributed here.
See `DATA_ACCESS.md` and `empirical_analysis/data/README.md`.

## Reproducibility

See `REPRODUCIBILITY.md` for the intended reproduction sequence and for the distinction between lightweight smoke
tests and full Monte Carlo runs.

## License

No public reuse license has been selected for this repository. Unless a license is added later, default copyright
restrictions apply.
