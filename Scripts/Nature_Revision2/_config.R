# ----------------------------------------------------------------------------
# Nature Revision 2 — carbon pipeline path helpers
#
# I keep NR2 output roots and small label-normalisation helpers here so the long carbon scripts stay readable.
# Inputs: optional environment variable NR2_OUT_ROOT; otherwise defaults under Outputs/Nature_Revision_Outputs/NR2/current.
# Outputs: directories per step via nr2_ensure_dirs; RUN_INFO.txt from nr2_write_run_info() when the runner calls it.
# ----------------------------------------------------------------------------

nr2_out_root <- Sys.getenv("NR2_OUT_ROOT")
if (identical(nr2_out_root, "")) {
  # Deterministic, simple default for publication/review
  nr2_out_root <- file.path("Outputs", "Nature_Revision_Outputs", "NR2", "current")
}

nr2_step_paths <- function(step, include_models = FALSE) {
  step_root <- file.path(nr2_out_root, step)
  paths <- list(
    step = step,
    root = step_root,
    figures = file.path(step_root, "figures"),
    tables = file.path(step_root, "tables"),
    rds = file.path(step_root, "rds")
  )

  if (isTRUE(include_models)) {
    paths$models <- file.path(step_root, "models")
  }

  paths
}

# MasterAllScenarios*.rds often use "once_logged" / "twice_logged" (underscores) while draws
# and routing logic use "once-logged" / "twice-logged". Normalize before schedule rules run.
nr2_normalize_habitat_labels <- function(x) {
  x <- trimws(as.character(x))
  x <- ifelse(x == "once_logged", "once-logged", x)
  ifelse(x == "twice_logged", "twice-logged", x)
}

nr2_ensure_dirs <- function(paths) {
  dir.create(paths$root, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$figures, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$tables, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$rds, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(paths$models)) dir.create(paths$models, recursive = TRUE, showWarnings = FALSE)
  invisible(paths)
}

nr2_write_run_info <- function() {
  dir.create(nr2_out_root, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      paste0("timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0("wd: ", getwd()),
      capture.output(sessionInfo())
    ),
    con = file.path(nr2_out_root, "RUN_INFO.txt")
  )
}
