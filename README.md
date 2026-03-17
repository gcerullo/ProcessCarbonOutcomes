# ProcessCarbonOutcomes (Nature Revision 2 pipeline)

This repository contains the scripts and data needed to reproduce the carbon modelling and scenario propagation analyses.

The **publication-ready, reproducible workflow** is in:

- `Scripts/Nature_Revision_2/`

All required input data are included in:

- `RawData/` (raw plot datasets)
- `Inputs/` (scenario tables + parameters + SCC temperature trajectory input)

## Quick start (run end-to-end)

From the repository root, run:

```bash
Rscript Scripts/Nature_Revision_2/run_all.R
```

Outputs are written deterministically to:

- `Outputs/Nature_Revision_Outputs/NR2/current/`

## Software requirements

- R (recent version recommended)
- Packages used across the pipeline include `tidyverse`, `brms`, `cmdstanr`, `tidybayes`, `data.table`, and several plotting/helper packages.
- **CmdStan** is required for the Bayesian model fits (via `cmdstanr`).

### Reproducible package environment (recommended)

This repository is configured with [`renv`](https://rstudio.github.io/renv/) for reproducible package versions.

In R, from the repository root:

```r
install.packages("renv")
renv::restore()
```

### CmdStan (required for model fitting)

If you do not already have CmdStan installed, in R run:

```r
cmdstanr::install_cmdstan()
```

## Pipeline steps (Nature Revision 2)

- **01** `Scripts/Nature_Revision_2/01_one_model_to_rule_them_all.R`
  - Fits the unified Bayesian recovery model and exports fitted model objects and prediction draws.
  - Reads: `RawData/SSB_plot_data.csv`, `RawData/Philipson20_PlotData2.csv`

- **02** `Scripts/Nature_Revision_2/02_prepare_draws.R`
  - Converts the fitted model into habitat-by-age draw tables (including twice-logged inference rules).
  - Reads: outputs from step 01

- **03** `Scripts/Nature_Revision_2/03_calculate_social_discount_rates.R`
  - Computes SCC discount-rate curves (2%, 4%, 6%).
  - Reads: `Inputs/Venmans_Value of an offset 25_1evening.csv`

- **04** `Scripts/Nature_Revision_2/04_propage_carbon_thru_scenarios.R`
  - Propagates carbon uncertainty through scenario schedules and monetises fluxes using SCC discount rates.
  - Reads: `Inputs/FixedScenarioParmams.R`, `Inputs/MasterAllScenarios.rds`, `Inputs/MasterAllScenarios_withHarvestDelays.rds`, plus outputs from steps 02 and 03.

## Input data notes

See:
- `RawData/README.md`
- `Inputs/README.md`

## Notes on older scripts

The repository may contain older scripts and outputs used during development. The **intended reproducing path for readers/reviewers** is the `Scripts/Nature_Revision_2/` pipeline described above.