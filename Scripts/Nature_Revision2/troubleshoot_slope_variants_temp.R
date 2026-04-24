# ----------------------------------------------------------------------------
# Nature Revision 2 — temporary slope-variant troubleshooting
#
# I keep this around for ad hoc checks when I experiment with twice-logged slope trajectories; it is not part of the main run_all sequence.
# Inputs: whichever partial RDS I point it at inside the file.
# Outputs: console plots or temporary RDS — delete or archive when the issue is closed.
# ----------------------------------------------------------------------------


# Temporary diagnostic: for mostly_2L starting scenarios, does full-horizon
# cumulative carbon stock (stock-years: sum of annual landscape ACD) differ
# between twice-logged slope trajectory 1 vs 1.2?
#
# Reads propagated outputs from 04_propagate tables (per-trajectory CSVs or
# all_trajectories combine).

#tldr; there is a differnce between 1 and 1.2 slopes; but its in the mean_abs_difference of ~40 mullion MgChayr; whereas plots are in billions 
#so not readily visible. 

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

source(file.path("Scripts", "Nature_Revision_2", "_config.R"))
paths <- nr2_step_paths("04_propagate")
tab_dir <- paths$tables

trajectory_file_suffix <- function(x) {
  y <- as.character(x)
  y <- gsub("[^A-Za-z0-9]+", "_", y)
  y <- gsub("_+", "_", y)
  y <- gsub("^_|_$", "", y)
  paste0("2Ltraj_", y)
}

SLOPE_A <- "1"
SLOPE_B <- "1.2"

# --- Load stock-years summaries for two slope assumptions --------------------
all_path <- file.path(tab_dir, "final_perf__carbon_stock_years__unc__all_trajectories.csv")
use_all <- file.exists(all_path)

if (use_all) {
  raw <- read.csv(all_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"twice_logged_slope_trajectory" %in% names(raw)) {
    stop("Expected column `twice_logged_slope_trajectory` in: ", all_path)
  }
  raw <- raw %>%
    mutate(twice_logged_slope_trajectory = as.character(twice_logged_slope_trajectory))
  d1 <- raw %>% filter(twice_logged_slope_trajectory == SLOPE_A)
  d2 <- raw %>% filter(twice_logged_slope_trajectory == SLOPE_B)
} else {
  f1 <- file.path(
    tab_dir,
    paste0("final_perf__carbon_stock_years__unc__", trajectory_file_suffix(SLOPE_A), ".csv")
  )
  f2 <- file.path(
    tab_dir,
    paste0("final_perf__carbon_stock_years__unc__", trajectory_file_suffix(SLOPE_B), ".csv")
  )
  if (!file.exists(f1) || !file.exists(f2)) {
    stop(
      "Missing per-trajectory stock-years tables.\n",
      "Tried:\n  ", f1, "\n  ", f2, "\n",
      "Or place `final_perf__carbon_stock_years__unc__all_trajectories.csv` in:\n  ",
      tab_dir
    )
  }
  d1 <- read.csv(f1, stringsAsFactors = FALSE, check.names = FALSE)
  d2 <- read.csv(f2, stringsAsFactors = FALSE, check.names = FALSE)
  d1$twice_logged_slope_trajectory <- SLOPE_A
  d2$twice_logged_slope_trajectory <- SLOPE_B
}

key_cols <- c("index", "production_target", "scenarioName", "scenarioStart")
stopifnot(all(key_cols %in% names(d1)), all(key_cols %in% names(d2)))
stopifnot("mean_cum_stock_year" %in% names(d1))

# --- Restrict to mostly twice-logged starting scenarios ------------------------
m2l <- function(x) {
  grepl("mostly_2L", as.character(x), ignore.case = TRUE)
}

d1m <- d1 %>% filter(m2l(scenarioStart))
d2m <- d2 %>% filter(m2l(scenarioStart))

message(
  "Rows after mostly_2L filter: slope ", SLOPE_A, " -> ", nrow(d1m),
  "; slope ", SLOPE_B, " -> ", nrow(d2m)
)

# --- Pairwise compare same scenario keys -------------------------------------
cmp <- inner_join(
  d1m %>%
    select(all_of(key_cols), mean_cum_stock_year) %>%
    rename(mean_cum_stock_year_slope_a = mean_cum_stock_year),
  d2m %>%
    select(all_of(key_cols), mean_cum_stock_year) %>%
    rename(mean_cum_stock_year_slope_b = mean_cum_stock_year),
  by = key_cols
) %>%
  mutate(
    diff_B_minus_A = mean_cum_stock_year_slope_b - mean_cum_stock_year_slope_a,
    rel_diff_pct = if_else(
      abs(mean_cum_stock_year_slope_a) > 1e-9,
      100 * diff_B_minus_A / abs(mean_cum_stock_year_slope_a),
      NA_real_
    )
  )

if (nrow(cmp) == 0) {
  stop(
    "No paired rows for mostly_2L after joining slopes ",
    SLOPE_A,
    " vs ",
    SLOPE_B,
    ". Check scenarioStart spelling and that both runs exist."
  )
}

# --- Summaries ----------------------------------------------------------------
summ <- tibble(
  n_pairs = nrow(cmp),
  mean_abs_diff = mean(abs(cmp$diff_B_minus_A), na.rm = TRUE),
  median_abs_diff = median(abs(cmp$diff_B_minus_A), na.rm = TRUE),
  max_abs_diff = max(abs(cmp$diff_B_minus_A), na.rm = TRUE),
  n_identical = sum(abs(cmp$diff_B_minus_A) < 1e-12, na.rm = TRUE),
  n_rel_gt_0p1_pct = sum(!is.na(cmp$rel_diff_pct) & abs(cmp$rel_diff_pct) > 0.1, na.rm = TRUE)
)

message("\n=== Paired comparison: cumulative stock-years (mean_cum_stock_year) ===")
message("Slope B = ", SLOPE_B, " minus slope A = ", SLOPE_A)
print(summ)

# Paired tests (nonparametric + parametric)
wt <- tryCatch(
  stats::wilcox.test(
    cmp$mean_cum_stock_year_slope_b,
    cmp$mean_cum_stock_year_slope_a,
    paired = TRUE,
    exact = FALSE
  ),
  error = function(e) NULL
)
tt <- tryCatch(
  stats::t.test(
    cmp$mean_cum_stock_year_slope_b,
    cmp$mean_cum_stock_year_slope_a,
    paired = TRUE
  ),
  error = function(e) NULL
)

if (!is.null(wt)) {
  message("\nWilcoxon signed-rank (paired): V = ", signif(wt$statistic, 5), ", p = ", format.pval(wt$p.value, digits = 4))
}
if (!is.null(tt)) {
  message(
    "Paired t-test: diff mean = ", signif(mean(cmp$diff_B_minus_A), 6),
    ", p = ", format.pval(tt$p.value, digits = 4)
  )
}

# Largest absolute differences (top 15)
top <- cmp %>%
  arrange(desc(abs(diff_B_minus_A))) %>%
  select(
    all_of(key_cols),
    mean_cum_stock_year_slope_a,
    mean_cum_stock_year_slope_b,
    diff_B_minus_A,
    rel_diff_pct
  ) %>%
  head(15)

message("\nTop rows by |diff| (B minus A):")
print(as.data.frame(top), row.names = FALSE)

out_csv <- file.path(tab_dir, "troubleshoot_slope_variants_temp__paired_m2l.csv")
write.csv(cmp, out_csv, row.names = FALSE)
message("\nWrote full paired table: ", normalizePath(out_csv, winslash = "/", mustWork = FALSE))
