# Adjusted-Range Self-Normalization for Functional Time Series

R code for **Adjusted-Range Self-Normalization for Sequential Monitoring of Functional Principal Component Score Means**.

This repository contains R source files, this README, and Movie S1. Manuscripts, LaTeX sources, numerical outputs, static figures, and replication archives are distributed separately with the article. The R scripts previously inside the public replication archive are now available directly under `reproducibility/code/`; their contents are unchanged.

## Main implementation

- `reproducibility/code/R/dual_channel_monitor.R`: LDL transformation, adjusted-range self-normalization, and full quadratic self-normalization.
- `reproducibility/code/R/intraday_interface.R`: conversion of minute OHLC observations into daily intraday volatility curves.
- `reproducibility/code/simulations/componentwise_range/`: critical-value, empirical-size, power, dimension-selection, and dependence experiments.
- `reproducibility/code/empirical/`: empirical analysis and scripts for producing the corresponding tables and figures.
- `reproducibility/code/tests/test_dual_channel_monitor.R`: implementation checks that do not require proprietary data.

The paper uses the boundary `(1 + s)^2`, corresponding to `gamma = 0`. Older filenames containing `dual_channel` or `shao` are retained for compatibility. The main retained-score implementation is the one under `reproducibility/code/`.

## Quick check

Run from the repository root:

```sh
Rscript reproducibility/code/tests/test_dual_channel_monitor.R
Rscript tests/parse_all_R_files.R
```

The implementation checks use synthetic inputs. They do not rerun the complete Monte Carlo study or verify archived empirical results.

## Simulations and empirical analysis

The scripts retain their original relative paths. Generate their required critical-value inputs before running simulations. For the paired range–quadratic experiment:

```sh
cd reproducibility/code
Rscript simulations/componentwise_range/run_critical_value_audit.R
Rscript simulations/componentwise_range/run_shao_comparison.R
```

Run empirical scripts from the repository root, following the input and command-line options defined in each script. Some empirical and plotting scripts require saved outputs from earlier analysis stages; supply the separate replication materials before running those scripts. The raw minute-level S&P 500 observations are commercially licensed from TickData and are not distributed here. Obtain the data from the provider and set `SPX_DATA_FILE` to the licensed input file when required.

The core implementation uses base R. The intraday interface additionally requires `data.table`; individual simulation and plotting scripts declare their further dependencies. Long simulation runs can be computationally intensive.

## Additional R code

The top-level `critical_values/`, `simulations/`, and `empirical_analysis/` directories retain previously published R code for functional-representation monitoring, including KS/weighted-CvM and alternative detectors. These retain their own method definitions and are distinct from the retained-score adjusted-range implementation above. `requirements.R` checks dependencies for those scripts.

## Supplementary video

[Movie S1: S&P 500 intraday volatility surfaces](Movie_S1_SPX_intraday_surface.mp4) is described in Supplement Section D.2 and the caption to Figure D.2. It shows all 500 monitoring observations in 2019–2020 and a closer view of January–June 2020. The adjusted-range and full quadratic panels use the same curves, viewing angles and scales. The displayed first alarms are 25 March and 4 May 2020, respectively, using asymptotic Brownian critical values; the video is distinct from the bootstrap comparison in the current main text.

The generation script is `reproducibility/code/empirical/build_spx_3d_multimedia.R`. Regeneration requires the separate saved empirical results, licensed minute observations and ffmpeg. The published video itself can be viewed without those inputs.

## Reproduction materials

Saved numerical results, empirical caches, and revision-specific materials are supplied separately with the manuscript. They are not included in this repository. Scripts may create local data files, tables, and figures when executed; those generated outputs are not included in the published source tree.

No public reuse license has been selected; existing copyright restrictions apply.
