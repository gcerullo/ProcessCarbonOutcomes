# ----------------------------------------------------------------------------
# Nature Revision 2 — carbon pipeline runner
#
# I source this when I want a clean-room rebuild of the NR2 carbon chain from model fit through discounting and propagation.
# Inputs: ProcessCarbonOutcomes project layout; Scripts/Nature_Revision_2/_config.R and the numbered step scripts it calls.
# Outputs: everything each step writes under nr2_out_root (see _config.R).
# ----------------------------------------------------------------------------

stopifnot(dir.exists("Scripts"), file.exists(file.path("ProcessCarbonOutcomes.Rproj")))

source(file.path("Scripts", "Nature_Revision_2", "_config.R"))
nr2_write_run_info()

scripts <- c(
  "01_one_model_to_rule_them_all.R",
  "02_prepare_draws.R",
  "03_calculate_social_discount_rates.R",
  "04_propage_carbon_thru_scenarios.R"
)

for (s in scripts) {
  message("\n==============================")
  message("Running: ", s)
  message("==============================")
  sys.source(file.path("Scripts", "Nature_Revision_2", s), envir = new.env(parent = globalenv()))
}

message("\nAll done.")
message("Outputs written to: ", normalizePath(nr2_out_root, winslash = "/", mustWork = FALSE))
