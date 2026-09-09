# Adjusted-Range Self-Normalization for Functional Time Series

Replication materials for **Adjusted-Range Self-Normalization for Sequential Monitoring of Functional Principal Component Score Means**, updated 9 September 2026.

## Current manuscript: download and run

The implementation, saved numerical results and Movie S1 are provided below. Download and extract the archives into the **same directory**, preserving their internal paths.

| Archive | Contents |
| --- | --- |
| [Replication_Code_and_Results.zip](Replication_Code_and_Results.zip) | R implementation, simulation and empirical scripts, CSV results, documentation and static figures |
| [Replication_Empirical_Continuous.zip](Replication_Empirical_Continuous.zip) | Saved objects and Brownian maxima for continuous 2019–2020 monitoring |
| [Replication_Empirical_Sensitivity.zip](Replication_Empirical_Sensitivity.zip) | Brownian caches for supplementary reference-window comparisons |
| [Movie_S1_SPX_intraday_surface.mp4](Movie_S1_SPX_intraday_surface.mp4) | Two-panel animation of intraday profiles and first monitoring alarms |

If downloading the movie separately, put it in the extracted package's `multimedia/` directory. The code-and-results archive includes its documentation and generation script.

From the extracted package root, run:

```sh
Rscript reproducibility/code/tests/test_dual_channel_monitor.R
Rscript reproducibility/code/empirical/check_spx_continuous_outputs.R
```

The first command checks the implementation; the second checks saved continuous-monitoring results. Neither requires licensed minute observations. These checks are not a rerun of the full Monte Carlo study. See the extracted `README.md`, `reproducibility/README.md` and `reproducibility/code/README.md` for the file-to-result map, dependencies and commands for full replication.

## Current empirical specification

FPCA is fitted on 2014–2015, reference quantities are estimated from 2016–2018, and monitoring continues through all 500 retained trading days in 2019–2020 without a January reset. Method-specific critical values use 250,000 matched-grid Brownian simulations. First alarms occur on 25 March 2020 for adjusted-range SN and 4 May 2020 for quadratic SN. These are first threshold crossings, not estimates of the date a change began. Earlier 2020-only analyses are retained for supplementary Appendix E.

## Data access

Raw licensed TickData observations are not redistributed. Researchers with the appropriate licence can obtain the source observations directly from TickData. Set `SPX_DATA_FILE` to the licensed input when running scripts that require it. Derived outputs are supplied for checking reported results; regenerating the curves and movie from minute observations requires licensed data access.

## Licence

No public reuse licence has been selected. Existing copyright restrictions remain in place.
