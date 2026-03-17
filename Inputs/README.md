# Inputs (redistributable)

This folder contains **all non-raw inputs** required to reproduce the Nature Revision 2 pipeline in `Scripts/Nature_Revision_2/`.

## Required files

- `Venmans_Value of an offset 25_1evening.csv`
  - Used by `Scripts/Nature_Revision_2/03_calculate_social_discount_rates.R`
  - Provides temperature trajectories by RCP (used to construct SCC discount-rate curves).

- `FixedScenarioParmams.R`
  - Used by `Scripts/Nature_Revision_2/04_propage_carbon_thru_scenarios.R`
  - Defines the starting landscape composition objects (e.g. `all_start_landscape`, and related scenario parameters).

- `MasterAllScenarios.rds`
  - Used by `Scripts/Nature_Revision_2/04_propage_carbon_thru_scenarios.R`
  - Scenario composition (used for plotting/summary joins).

- `MasterAllScenarios_withHarvestDelays.rds`
  - Used by `Scripts/Nature_Revision_2/04_propage_carbon_thru_scenarios.R`
  - Full scenario schedules (including harvest-delay variants).

If you rename or relocate any of these, update the corresponding scripts.

