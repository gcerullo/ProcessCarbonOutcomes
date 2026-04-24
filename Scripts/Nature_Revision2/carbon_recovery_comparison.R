#!/usr/bin/env Rscript

# =============================================================================
# Nature Revision 2 — Carbon recovery comparison table
# =============================================================================
#
# Purpose
# - Build a compact, uncertainty-aware table of average annual ACD increase
#   suitable for comparison against values reported in the broader literature.
# - Uses posterior draws from the final carbon model output (02_draws).
#
# Inputs
# - Posterior / scenario carbon RDS from the NR2 carbon draws pipeline (see sourcing below).
#
# Outputs
# - Tables + RDS under: Outputs/Nature_Revision_Outputs/NR2/current/carbon_recovery_comparison

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

source(file.path("Scripts", "Nature_Revision_2", "_config.R"))
paths <- nr2_step_paths("carbon_recovery_comparison")
nr2_ensure_dirs(paths)

out_root <- paths$root
tab_dir <- paths$tables
rds_dir <- paths$rds

log_path <- file.path(out_root, "run_log.txt")
log_line <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(txt, "\n", file = log_path, append = TRUE)
  message(txt)
}

warn_if_windows_path_long <- function(path, warn_limit = 240, hard_limit = 259) {
  if (!identical(.Platform$OS.type, "windows")) return(invisible(path))
  full_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  n <- nchar(full_path)
  if (n >= hard_limit) {
    log_line("WARNING: very long output path (", n, " chars): ", full_path)
  } else if (n >= warn_limit) {
    log_line("Warning: long output path (", n, " chars): ", full_path)
  }
  invisible(path)
}

write_csv_checked <- function(df, path) {
  warn_if_windows_path_long(path)
  write.csv(as.data.frame(df), path, row.names = FALSE)
}

save_rds_checked <- function(object, path) {
  warn_if_windows_path_long(path)
  saveRDS(object, path)
}

# -----------------------------
# Run settings
# -----------------------------
# We report rates from time 0 to each horizon endpoint so literature comparisons
# can be done at short/medium/long horizons.
AGE_START <- 0
WINDOW_ENDS <- c(20, 40, 60)
TWO_L_SLOPES <- c("0.8", "1", "1.2")
PHILIPSON_WINDOWS <- tibble::tibble(window_start = c(30), window_end = c(35))

log_line("nr2_out_root = ", normalizePath(nr2_out_root, winslash = "/", mustWork = FALSE))
log_line("AGE_START = ", AGE_START)
log_line("WINDOW_ENDS = ", paste(WINDOW_ENDS, collapse = ", "))
log_line("TWO_L_SLOPES = ", paste(TWO_L_SLOPES, collapse = ", "))
log_line(
  "PHILIPSON_WINDOWS = ",
  paste0(PHILIPSON_WINDOWS$window_start, "-", PHILIPSON_WINDOWS$window_end, collapse = ", ")
)

if (!dir.exists(tab_dir) || !dir.exists(rds_dir)) {
  stop(
    "Output directories not available. Check NR2_OUT_ROOT path length and permissions.\n",
    "tab_dir: ", normalizePath(tab_dir, winslash = "/", mustWork = FALSE), "\n",
    "rds_dir: ", normalizePath(rds_dir, winslash = "/", mustWork = FALSE)
  )
}

draws_path <- file.path(nr2_out_root, "02_draws", "rds", "acdraws_aboveground.rds")
stopifnot(file.exists(draws_path))
acd_draws <- readRDS(draws_path)

# Harmonize habitat labels to avoid underscore/hyphen mismatches.
acd_draws <- acd_draws %>%
  mutate(
    habitat = as.character(habitat),
    habitat = if_else(habitat == "once_logged", "once-logged", habitat),
    habitat = if_else(habitat == "twice_logged", "twice-logged", habitat),
    slope_factor = as.character(slope_factor)
  )

# Recovery trajectories used for literature comparison.
# Use non-"*_start" habitat labels only, per request, so estimates map directly
# to the named habitat states used in figures/tables.
recovery_habitats <- c("primary", "once-logged", "restored", "twice-logged")

missing_recovery_habitats <- setdiff(recovery_habitats, unique(acd_draws$habitat))
if (length(missing_recovery_habitats) > 0) {
  stop(
    "Missing expected recovery habitat(s) in draws: ",
    paste(missing_recovery_habitats, collapse = ", ")
  )
}

# Keep slope-specific trajectories for twice-logged; keep baseline slope
# for primary/once-logged/restored.
acd_recovery <- acd_draws %>%
  filter(habitat %in% recovery_habitats) %>%
  filter(
    (habitat == "twice-logged" & slope_factor %in% TWO_L_SLOPES) |
      (habitat %in% c("primary", "once-logged", "restored") & slope_factor == "1")
  )

# Compute draw-level annualized ACD change between AGE_START and each window endpoint.
# For each habitat/slope/draw/window:
#   1) keep age range [AGE_START, window_end]
#   2) take mean ACD at the first and last available age within that range
#   3) annualize endpoint difference by years spanned
window_specs <- tibble(window_start = AGE_START, window_end = WINDOW_ENDS)

per_draw_rates <- bind_rows(
  lapply(seq_len(nrow(window_specs)), function(i) {
    w0 <- window_specs$window_start[i]
    w1 <- window_specs$window_end[i]

    acd_recovery %>%
      filter(functionalhabAge >= w0, functionalhabAge <= w1) %>%
      group_by(habitat, slope_factor, draw) %>%
      group_modify(~ {
        d <- .x %>% arrange(functionalhabAge)
        ages <- sort(unique(d$functionalhabAge))
        if (length(ages) < 2) {
          return(tibble(
            acd_start = NA_real_,
            acd_end = NA_real_,
            age_start = NA_real_,
            age_end = NA_real_,
            years_span = NA_real_,
            annual_acd_increase = NA_real_
          ))
        }

        age_start <- min(ages, na.rm = TRUE)
        age_end <- max(ages, na.rm = TRUE)
        acd_start <- d %>% filter(functionalhabAge == age_start) %>% summarise(v = mean(ACD, na.rm = TRUE), .groups = "drop") %>% pull(v)
        acd_end <- d %>% filter(functionalhabAge == age_end) %>% summarise(v = mean(ACD, na.rm = TRUE), .groups = "drop") %>% pull(v)
        years_span <- age_end - age_start

        tibble(
          acd_start = acd_start,
          acd_end = acd_end,
          age_start = age_start,
          age_end = age_end,
          years_span = years_span,
          annual_acd_increase = ifelse(years_span > 0, (acd_end - acd_start) / years_span, NA_real_)
        )
      }) %>%
      ungroup() %>%
      mutate(window_start = w0, window_end = w1)
  })
)

summary_table <- per_draw_rates %>%
  group_by(habitat, slope_factor, window_start, window_end) %>%
  summarise(
    n_draws = sum(!is.na(annual_acd_increase)),
    age_start = first(age_start),
    age_end = first(age_end),
    years_span = first(years_span),
    annual_acd_increase_mean = mean(annual_acd_increase, na.rm = TRUE),
    annual_acd_increase_median = median(annual_acd_increase, na.rm = TRUE),
    annual_acd_increase_lwr95 = quantile(annual_acd_increase, 0.025, na.rm = TRUE),
    annual_acd_increase_upr95 = quantile(annual_acd_increase, 0.975, na.rm = TRUE),
    annual_acd_increase_lwr80 = quantile(annual_acd_increase, 0.10, na.rm = TRUE),
    annual_acd_increase_upr80 = quantile(annual_acd_increase, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    habitat_for_comparison = case_when(
      habitat == "primary" ~ "primary",
      habitat == "once-logged" ~ "once-logged",
      habitat == "restored" ~ "restored",
      habitat == "twice-logged" ~ "twice-logged",
      TRUE ~ habitat
    ),
    trajectory_assumption = case_when(
      habitat == "twice-logged" & slope_factor == "0.8" ~ "twice-logged slope = 0.8",
      habitat == "twice-logged" & slope_factor == "1" ~ "twice-logged slope = 1.0",
      habitat == "twice-logged" & slope_factor == "1.2" ~ "twice-logged slope = 1.2",
      habitat == "primary" ~ "reference/background trajectory",
      habitat == "once-logged" ~ "baseline model trajectory",
      habitat == "restored" ~ "baseline model trajectory",
      TRUE ~ "baseline model trajectory"
    )
  ) %>%
  select(
    habitat_for_comparison, trajectory_assumption, slope_factor,
    window_start, window_end, n_draws, age_start, age_end, years_span,
    annual_acd_increase_mean, annual_acd_increase_median,
    annual_acd_increase_lwr80, annual_acd_increase_upr80,
    annual_acd_increase_lwr95, annual_acd_increase_upr95
  ) %>%
  arrange(habitat_for_comparison, window_end, slope_factor)

summary_csv <- file.path(tab_dir, "carbon_recovery__annual_acd_inc__windowed_summary.csv")
per_draw_csv <- file.path(tab_dir, "carbon_recovery__annual_acd_inc__windowed_per_draw.csv")
summary_rds <- file.path(rds_dir, "carbon_recovery__annual_acd_inc__windowed_summary.rds")
per_draw_rds <- file.path(rds_dir, "carbon_recovery__annual_acd_inc__windowed_per_draw.rds")

write_csv_checked(summary_table, summary_csv)
write_csv_checked(per_draw_rates, per_draw_csv)

save_rds_checked(summary_table, summary_rds)
save_rds_checked(per_draw_rates, per_draw_rds)

# -----------------------------------------------------------------------------
# Philipson-comparable recovery rates (DIRECT FROM MODEL SLOPES; no *_start)
# -----------------------------------------------------------------------------
# The Philipson comparison is fundamentally about slope (annual recovery rate)
# for natural regeneration vs active restoration. In this model:
#   ACD ~ 0 + state + 0 + state:time
# the annual rate for each state is directly the corresponding state:time
# coefficient. We therefore extract posterior draws of:
#   - once-logged slope: b_stateonce_logged:time
#   - restored slope:    b_staterestored:time
# and their contrast (restored - once-logged), without using any *_start labels.
joint_model_path <- file.path(nr2_out_root, "01_one_model", "models", "unified_linear_carbon_model.rds")
if (!file.exists(joint_model_path)) {
  stop("Could not find joint model at: ", joint_model_path)
}
if (!requireNamespace("posterior", quietly = TRUE)) {
  stop("Package `posterior` is required for draw extraction. Please install.packages('posterior').")
}

joint_model <- readRDS(joint_model_path)
draws_df <- as.data.frame(posterior::as_draws_df(joint_model))

coef_once <- "b_stateonce_logged:time"
coef_rest <- "b_staterestored:time"
if (!(coef_once %in% names(draws_df) && coef_rest %in% names(draws_df))) {
  stop(
    "Expected slope coefficients not found in model draws. Missing one of: ",
    coef_once, ", ", coef_rest
  )
}

philipson_per_draw <- tibble(
  draw = seq_len(nrow(draws_df)),
  once_logged_rate = as.numeric(draws_df[[coef_once]]),
  restored_rate = as.numeric(draws_df[[coef_rest]])
) %>%
  mutate(
    restored_minus_once_logged = restored_rate - once_logged_rate,
    window_start = PHILIPSON_WINDOWS$window_start[1],
    window_end = PHILIPSON_WINDOWS$window_end[1],
    comparison_mode = "philipson_like_model_direct"
  )

summ_tbl <- function(x, habitat_name, slope_name) {
  tibble(
    habitat_for_comparison = habitat_name,
    slope_factor = slope_name,
    n_draws = sum(!is.na(x)),
    annual_acd_increase_mean = mean(x, na.rm = TRUE),
    annual_acd_increase_median = median(x, na.rm = TRUE),
    annual_acd_increase_lwr80 = as.numeric(quantile(x, 0.10, na.rm = TRUE)),
    annual_acd_increase_upr80 = as.numeric(quantile(x, 0.90, na.rm = TRUE)),
    annual_acd_increase_lwr95 = as.numeric(quantile(x, 0.025, na.rm = TRUE)),
    annual_acd_increase_upr95 = as.numeric(quantile(x, 0.975, na.rm = TRUE))
  )
}

philipson_summary <- bind_rows(
  summ_tbl(philipson_per_draw$once_logged_rate, "once-logged", "1"),
  summ_tbl(philipson_per_draw$restored_rate, "restored", "1"),
  summ_tbl(philipson_per_draw$restored_minus_once_logged, "restored-minus-once-logged", "contrast") %>%
    mutate(p_restored_gt_once_logged = mean(philipson_per_draw$restored_minus_once_logged > 0, na.rm = TRUE))
) %>%
  mutate(
    comparison_mode = "philipson_like_model_direct",
    window_start = PHILIPSON_WINDOWS$window_start[1],
    window_end = PHILIPSON_WINDOWS$window_end[1]
  ) %>%
  select(
    comparison_mode, habitat_for_comparison, slope_factor,
    window_start, window_end, n_draws,
    annual_acd_increase_mean, annual_acd_increase_median,
    annual_acd_increase_lwr80, annual_acd_increase_upr80,
    annual_acd_increase_lwr95, annual_acd_increase_upr95,
    everything()
  )

philipson_summary_csv <- file.path(tab_dir, "carbon_recovery__philipson_like__summary.csv")
philipson_per_draw_csv <- file.path(tab_dir, "carbon_recovery__philipson_like__per_draw.csv")
philipson_summary_rds <- file.path(rds_dir, "carbon_recovery__philipson_like__summary.rds")
philipson_per_draw_rds <- file.path(rds_dir, "carbon_recovery__philipson_like__per_draw.rds")

write_csv_checked(philipson_summary, philipson_summary_csv)
write_csv_checked(philipson_per_draw, philipson_per_draw_csv)
save_rds_checked(philipson_summary, philipson_summary_rds)
save_rds_checked(philipson_per_draw, philipson_per_draw_rds)

# -----------------------------------------------------------------------------
# Quick visual comparison: model-direct slopes vs Philipson reference values
# -----------------------------------------------------------------------------
# This plot is for at-a-glance comparison of annual ACD increase estimates.
# - "your_model_direct" uses posterior fixed-effect slopes from the fitted
#   joint model.
# - "philipson_reference" uses headline values reported in Philipson et al. 2020
#   for naturally regenerating and restored forest during years 30-35 post-logging.
# Primary is shown as 0 in both series (no recovery slope expected for unlogged).
fx <- summary(joint_model)[["fixed"]]

model_direct_tbl <- tibble(
  habitat = c("primary", "once-logged", "restored"),
  source = "your_model_direct",
  estimate = c(
    fx["stateprimary:time", "Estimate"],
    fx["stateonce_logged:time", "Estimate"],
    fx["staterestored:time", "Estimate"]
  ),
  lwr95 = c(
    fx["stateprimary:time", "l-95% CI"],
    fx["stateonce_logged:time", "l-95% CI"],
    fx["staterestored:time", "l-95% CI"]
  ),
  upr95 = c(
    fx["stateprimary:time", "u-95% CI"],
    fx["stateonce_logged:time", "u-95% CI"],
    fx["staterestored:time", "u-95% CI"]
  )
)

philipson_reference_tbl <- tibble(
  habitat = c("primary", "once-logged", "restored"),
  source = "philipson_reference",
  estimate = c(0.0, 2.9, 4.4),
  lwr95 = c(NA_real_, 2.1, 3.6),
  upr95 = c(NA_real_, 3.7, 5.2)
)

comparison_plot_tbl <- bind_rows(model_direct_tbl, philipson_reference_tbl) %>%
  mutate(
    habitat = factor(habitat, levels = c("primary", "once-logged", "restored")),
    source = factor(source, levels = c("your_model_direct", "philipson_reference"))
  )

comparison_csv <- file.path(tab_dir, "carbon_recovery__model_vs_philipson__annual_slope.csv")
comparison_rds <- file.path(rds_dir, "carbon_recovery__model_vs_philipson__annual_slope.rds")
comparison_png <- file.path(paths$figures, "carbon_recovery__model_vs_philipson__annual_slope.png")
comparison_pdf <- file.path(paths$figures, "carbon_recovery__model_vs_philipson__annual_slope.pdf")

write_csv_checked(comparison_plot_tbl, comparison_csv)
save_rds_checked(comparison_plot_tbl, comparison_rds)

p_compare <- ggplot(comparison_plot_tbl, aes(x = habitat, y = estimate, color = source)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(
    aes(ymin = lwr95, ymax = upr95),
    width = 0.12,
    position = position_dodge(width = 0.45),
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  geom_point(
    position = position_dodge(width = 0.45),
    size = 2.8
  ) +
  scale_color_manual(
    values = c("your_model_direct" = "#1f78b4", "philipson_reference" = "#d95f02"),
    labels = c("your_model_direct" = "Your model (posterior slope)", "philipson_reference" = "Philipson 2020 reference")
  ) +
  labs(
    x = NULL,
    y = "Annual ACD increase (Mg C ha^-1 yr^-1)",
    color = NULL,
    title = "Annual ACD recovery: your model vs Philipson reference"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

warn_if_windows_path_long(comparison_png)
warn_if_windows_path_long(comparison_pdf)
ggsave(comparison_png, p_compare, width = 8, height = 5.2, dpi = 220)
ggsave(comparison_pdf, p_compare, width = 8, height = 5.2)

writeLines(capture.output(sessionInfo()), con = file.path(out_root, "sessionInfo.txt"))

log_line("Saved summary table: ", summary_csv)
log_line("Saved per-draw table: ", per_draw_csv)
log_line("Saved Philipson-like summary: ", philipson_summary_csv)
log_line("Saved Philipson-like per-draw: ", philipson_per_draw_csv)
log_line("Saved quick comparison table: ", comparison_csv)
log_line("Saved quick comparison figure: ", comparison_png)
log_line("Done.")

