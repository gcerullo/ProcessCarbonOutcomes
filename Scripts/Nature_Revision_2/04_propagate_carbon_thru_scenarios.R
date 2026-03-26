# =============================================================================
# Nature Revision 2 — Step 04: propagate carbon through scenarios
# =============================================================================
#
# Purpose
# - Main script for propagating carbon through scenarios with posterior draws.
# - Keeps the modular structure so each stage is easier to verify and maintain.
#
# NOTE ON COMMENTS
# - You asked to keep the original narrative comments because they capture the rationale
#   and make it easier to communicate what each stage is doing.
# - Where helpful, the blocks below reuse your original wording (sometimes adapted to the
#   refactored function names).

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyverse)
  library(data.table)
})

if (!requireNamespace("biscale", quietly = TRUE)) {
  stop("Package `biscale` is required for bivariate scenario colouring. Please run: install.packages('biscale')")
}

# =============================================================================
# 0) OUTPUTS + RUN SETTINGS
# =============================================================================

source(file.path("Scripts", "Nature_Revision_2", "_config.R"))
paths <- nr2_step_paths("04_propagate")
nr2_ensure_dirs(paths)

# Base output folder for this step
out_root <- paths$root
fig_dir <- paths$figures
tab_dir <- paths$tables
rds_dir <- paths$rds

log_path <- file.path(out_root, "run_log.txt")
# Lightweight logger: mirror status to console and persistent run log.
log_line <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(txt, "\n", file = log_path, append = TRUE)
  message(txt)
}

# -----------------------------
# Main knobs for fast iteration
# -----------------------------
TEST_MODE <- FALSE

# Cache the slowest step (scenario schedules prep) under the temp output folder.
# This makes repeated iterations much faster.
USE_CACHE <- FALSE

# Scenario-set selection
# - A "scenario set" = one element of `Inputs/MasterAllScenarios_withHarvestDelays.rds` (a list).
# - In full runs, you typically process ALL sets.
scenario_set_idxs <- seq(1,12)     # e.g. c(1, 5, 9). If NULL: auto (TEST_MODE) or all (full).

# Within-set filtering (optional): keep only a few scenario `index` values per set
limit_indices_per_set <- NULL    # e.g. 5 (keep first 5 unique index values)

# Posterior draw subsetting (optional): use fewer draws for smoke tests
max_draws <- 200               # e.g. 50 (keep first 50 draws)

# Twice-logged recovery trajectory choice(s)
# - Set to one value (e.g., 1) for a single run.
# - Set multiple values (e.g., c(0.75, 1, 1.25)) for an outer loop.
# - Set NULL to auto-run all slope_factor values present in ACD draws.
twice_logged_slope_trajectories <- c(0.8, 1, 1.2)

set.seed(123)

if (isTRUE(TEST_MODE)) {
  # Very small default smoke-test (fast iteration).
  # Increase these gradually as you gain confidence.
  if (is.null(scenario_set_idxs)) scenario_set_idxs <- 1
 # if (is.null(limit_indices_per_set)) limit_indices_per_set <- 2
  if (is.null(max_draws)) max_draws <- 10
}

log_line("TEST_MODE = ", TEST_MODE)
log_line("scenario_set_idxs = ", paste(scenario_set_idxs %||% "ALL", collapse = ", "))
log_line("limit_indices_per_set = ", limit_indices_per_set %||% "NULL")
log_line("max_draws = ", max_draws %||% "NULL")
log_line("twice_logged_slope_trajectories = ", paste(twice_logged_slope_trajectories %||% "AUTO", collapse = ", "))
log_line("nr2_out_root = ", normalizePath(nr2_out_root, winslash = "/", mustWork = FALSE))
log_line("USE_CACHE = ", USE_CACHE)

cache_dir <- file.path(out_root, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
cache_scenarios_prepped <- file.path(cache_dir, "scenarios_prepped.rds")
cache_version <- "04_propagate_scenarios_v1"

# Utility timer/logger so long stages report start/finish and runtime.
time_it <- function(label, expr) {
  t0 <- Sys.time()
  log_line(label, " ...")
  out <- force(expr)
  dt <- difftime(Sys.time(), t0, units = "secs")
  log_line(label, " done (", sprintf("%.1f", as.numeric(dt)), "s)")
  out
}

# Filename-safe suffix to label outputs by 2L trajectory value.
trajectory_file_suffix <- function(x) {
  y <- as.character(x)
  y <- gsub("[^A-Za-z0-9]+", "_", y)
  y <- gsub("_+", "_", y)
  y <- gsub("^_|_$", "", y)
  paste0("2Ltraj_", y)
}

# =============================================================================
# 1) INPUTS
# =============================================================================
# Standard ordering helper so schedules are easier to inspect/debug.
scenario_row_order_visualisation_fun <- function(x){
  x %>%  
    group_by(index, production_target, habitat, original_habitat, harvest_delay) %>%  
    arrange(true_year, .by_group = TRUE) %>%
    ungroup()
}

# Load static scenario metadata and starting landscape composition.
load_fixed_scenario_params <- function(path = file.path("Inputs", "FixedScenarioParmams.R")) {
  stopifnot(file.exists(path))
  # Use an isolated environment but inherit from the current session so that
  # base functions + attached packages are available to the sourced file.
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)

  if (!exists("all_start_landscape", envir = env, inherits = FALSE)) {
    stop("Expected `all_start_landscape` to be defined by Inputs/FixedScenarioParmams.R")
  }

  list(all_start_landscape = get("all_start_landscape", envir = env, inherits = FALSE))
}
# Expand baseline starting landscape across years for join-ready trajectories.
expand_start_landscape_through_time <- function(all_start_landscape, years = 0:60) {
  all_start_landscape %>%
    mutate(functional_habitat = habitat) %>%
    mutate(
      functional_habitat = case_when(
        habitat == "once-logged" ~ "once-logged_start",
        habitat == "twice-logged" ~ "twice-logged_start",
        TRUE ~ habitat
      )
    ) %>%
    tidyr::uncount(weights = length(years)) %>%
    group_by(across(everything())) %>%
    mutate(functionalhabAge = years) %>%
    ungroup()
}
# Load SCC discount-rate table and derive normalized SCC ratios by year.
load_social_discount_rates <- function(path = file.path(nr2_out_root, "03_scc", "tables", "scc_dr_2_4_6.csv")) {
  stopifnot(file.exists(path))
  read.csv(path) %>%
    group_by(discount_rate) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(scc_discounted_ratio = scc_discounted / first(scc_discounted)) %>%
    ungroup()
}
# Build scenario composition descriptors used in summaries/plots.
load_scenario_composition <- function(path = file.path("Inputs", "MasterAllScenarios.rds")) {
  stopifnot(file.exists(path))
  scenarios2 <- readRDS(path)
  scen_dt <- rbindlist(scenarios2, use.names = TRUE) %>%
    # match your older plotting script: infer starting landscape from scenarioName
    mutate(
      scenarioStart = scenarioName,
      scenarioStart = stringr::str_remove(scenarioStart, "_IY_ND\\.csv$"),
      scenarioStart = stringr::str_remove(scenarioStart, "_CY_ND\\.csv$"),
      scenarioStart = stringr::str_remove(scenarioStart, "_IY_D\\.csv$"),
      scenarioStart = stringr::str_remove(scenarioStart, "_CY_D\\.csv$")
    )

  # Full habitat composition: one row per (scenarioName, index, production_target)
  hab_props <- scen_dt %>%
    group_by(scenarioName, scenarioStart, index, production_target, habitat) %>%
    summarise(prop = sum(num_parcels, na.rm = TRUE) / 1000, .groups = "drop")

  # Make the composition display stable/readable in hovers
  hab_props_ord <- hab_props %>%
    mutate(
      hab_order = match(
        habitat,
        c(
          "primary",
          "once-logged",
          "twice-logged",
          "restored",
          "deforested",
          "albizia_current",
          "albizia_future",
          "eucalyptus_current",
          "eucalyptus_future"
        )
      )
    ) %>%
    arrange(is.na(hab_order), hab_order, habitat)

  comp_tbl <- hab_props_ord %>%
    group_by(scenarioName, scenarioStart, index, production_target) %>%
    summarise(
      comp_sum = sum(prop, na.rm = TRUE),
      comp_html = paste0(
        "<b>composition</b><br>",
        paste0(habitat, ": ", signif(prop, 3), collapse = "<br>")
      ),
      .groups = "drop"
    )

  out <- hab_props %>%
    group_by(scenarioName, scenarioStart, index, production_target) %>%
    summarise(
      # proportions of TOTAL landscape [1000 parcels] in different habitat type
      propOG = sum(prop[habitat == "primary"], na.rm = TRUE),
      prop1L = sum(prop[habitat == "once-logged"], na.rm = TRUE),
      prop2L = sum(prop[habitat == "twice-logged"], na.rm = TRUE),
      propDeforested = sum(prop[habitat == "deforested"], na.rm = TRUE),
      propRestored = sum(prop[habitat == "restored"], na.rm = TRUE),
      propEucalyptus = sum(prop[grepl("^eucalyptus", habitat)], na.rm = TRUE),
      propAlbizia = sum(prop[grepl("^albizia", habitat)], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(propPlant = propEucalyptus + propAlbizia) %>%
    left_join(comp_tbl, by = c("scenarioName", "scenarioStart", "index", "production_target"))

  # starting landscape habitat proportions (same mapping you used previously)
  habInStart <- tibble::tibble(
    scenarioStart = c(
      "all_primary", "mostly_1L", "mostly_2L",
      "primary_deforested", "mostly_1L_deforested", "mostly_2L_deforested"
    ),
    originalOG = c(1, 0.2, 0.2, 0.8, 0.2, 0.2),
    original1L = c(0, 0.8, 0, 0, 0.6, 0),
    original2L = c(0, 0, 0.8, 0, 0, 0.6)
  )

  out <- out %>%
    left_join(habInStart, by = "scenarioStart") %>%
    mutate(
      remainingOG = if_else(!is.na(originalOG) & originalOG > 0, propOG / originalOG, 0),
      remaining1L = if_else(!is.na(original1L) & original1L > 0, prop1L / original1L, 0),
      remaining2L = if_else(!is.na(original2L) & original2L > 0, prop2L / original2L, 0)
    )

  # bivariate palette (BlueOr; dim=4) + hex colours
  hex_vec <- biscale::bi_pal("BlueOr", dim = 4, preview = FALSE)
  cols <- tibble::tibble(bi_class = names(hex_vec), hex = as.character(hex_vec))

  out_prim <- out %>%
    filter(scenarioStart %in% c("all_primary", "primary_deforested")) %>%
    biscale::bi_class(x = propOG, y = prop1L, dim = 4, style = "equal") %>%
    left_join(cols, by = "bi_class") %>%
    rename(hexP = hex)

  out_1l <- out %>%
    filter(scenarioStart %in% c("mostly_1L", "mostly_1L_deforested")) %>%
    biscale::bi_class(x = remainingOG, y = remaining1L, dim = 4, style = "equal") %>%
    left_join(cols, by = "bi_class") %>%
    rename(hex1L = hex)

  out_2l <- out %>%
    filter(scenarioStart %in% c("mostly_2L", "mostly_2L_deforested")) %>%
    biscale::bi_class(x = remainingOG, y = remaining2L, dim = 4, style = "equal") %>%
    left_join(cols, by = "bi_class") %>%
    rename(hex2L = hex)

  out_other <- out %>%
    filter(!scenarioStart %in% c(
      "all_primary", "primary_deforested",
      "mostly_1L", "mostly_1L_deforested",
      "mostly_2L", "mostly_2L_deforested"
    ))

  out <- bind_rows(out_prim, out_1l, out_2l, out_other)
  rm(scenarios2)
  out
}

# Load scenario schedules (with harvest-delay variants) from master input RDS.
load_scenario_schedules <- function(path = file.path("Inputs", "MasterAllScenarios_withHarvestDelays.rds")) {
  stopifnot(file.exists(path))
  readRDS(path)
}

# Load ACD posterior draws prepared in step 02.
load_acd_draws <- function(path = file.path(nr2_out_root, "02_draws", "rds", "acdraws_aboveground.rds")) {
  stopifnot(file.exists(path))
  readRDS(path)
}

# =============================================================================
# 2) SCENARIO PREP (harvest-delay rules etc)
# =============================================================================

# ---------------------------------------------------------------------------
# transition rules  (original notes)
# ---------------------------------------------------------------------------
# primary - no rules
#
# once logged
# a. P -> 1L  If parcel goes to from primary to once-logged, we start at 1L yr 0
# b. 1L -> 1L if parcel starts as once-logged and stays 1L, keep harvest_delay == 0 (no transition)
# c. 1L -> 2L if parcel starts as once logged before re-harvest, we assume 1L recovery
#    (already explicit in scenarios) that resets to yr 0 once-logged values.
#
# 1L -> restored: route to restored_start curves; keep harvest_delay == 0 at scenario start
#
# if parcel starts as twice-logged, use the 15yr after 1L curves - with different slope factors,
# tracking the different slope assumption rates for twice-logged forest
# if parcel starts as primary and goes to twice-logged, use the 30yrafter 1L, with different slope factors
# ---------------------------------------------------------------------------

# For no-transition rows (original_habitat == habitat), keep only delay 0.
filter_original_eq_habitat <- function(df_list) {
  # Other original==habitat pairs (e.g. primary-primary); 1L->1L delay==0 is enforced in update_1L_functional_habitat
  df_list %>% map(~ .x %>% filter(!(original_habitat == habitat & (is.na(harvest_delay) | harvest_delay != 0))))
}

# update_1L_functional_habitat: see prepare_scenarios() — habitat labels are normalized *before* this runs.
# Rules only:
#   (1) original once-logged + habitat restored  -> functional_habitat = restored_start, functionalhabAge = true_year;
#       drop rows unless harvest_delay == 0.
#   (2) original once-logged + habitat once-logged -> functional_habitat = once-logged_start;
#       drop rows unless harvest_delay == 0.

# Apply once-logged-specific routing rules so schedule labels match draw habitats.
update_1L_functional_habitat <- function(df_list) {
  df_list %>%
    map(~ .x %>%
          mutate(
            functional_habitat = as.character(functional_habitat),
            functional_habitat = dplyr::case_when(
              original_habitat == "once-logged" & habitat == "restored" ~ "restored_start",
              original_habitat == "once-logged" & habitat == "once-logged" ~ "once-logged_start",
              TRUE ~ functional_habitat
            ),
            functionalhabAge = dplyr::if_else(
              original_habitat == "once-logged" & habitat == "restored",
              as.numeric(true_year),
              as.numeric(functionalhabAge)
            )
          ) %>%
          filter(
            !((original_habitat == "once-logged" & habitat == "restored") &
              (is.na(harvest_delay) | as.numeric(harvest_delay) != 0)),
            !((original_habitat == "once-logged" & habitat == "once-logged") &
              (is.na(harvest_delay) | as.numeric(harvest_delay) != 0))
          ))
}

# Apply twice-logged-specific routing and minimum-delay rule.
update_2L_functional_habitat <- function(df_list) {
  # scenarios that start as twice-logged are treated as re-harvested just before scenario start
  df_list %>%
    map(~ .x %>%
          mutate(
            functional_habitat = if_else(
              original_habitat == "twice-logged" & functional_habitat == "twice-logged",
              "twice_logged_start",
              functional_habitat
            )
          ) %>%
          filter(!(original_habitat == "twice-logged" & harvest_delay < 15))
    )
}

# Drop scenarios with production_target <= 0 (not relevant for timber pathways).
remove_no_timber_scenarios <- function(df_list) {
  df_list %>% map(~ .x %>% filter(production_target > 0))
}

# Prepare schedules once: normalize labels, coerce delays, apply transition rules.
prepare_scenarios <- function(scenarios) {
  scenarios <- scenarios %>%
    map(~ .x %>%
          mutate(harvest_delay = as.numeric(stringr::str_extract(as.character(harvest_delay), "\\d+"))) %>%
          mutate(
            original_habitat = nr2_normalize_habitat_labels(original_habitat),
            habitat = nr2_normalize_habitat_labels(habitat)
          ) %>%
          mutate(across(any_of("functional_habitat"), nr2_normalize_habitat_labels)))

  scenarios <- scenarios %>%
    update_1L_functional_habitat() %>%
    update_2L_functional_habitat() %>%
    filter_original_eq_habitat() %>%
    remove_no_timber_scenarios()

  lapply(scenarios, as.data.table)
}

# =============================================================================
# 3) DRAWS PREP
# =============================================================================

# Build per-draw lookup tables, filtered to the chosen twice-logged slope trajectory.
prepare_draw_list <- function(acd_draws, twice_logged_slope_trajectory, max_draws = NULL) {
  # Match legacy naming and select the desired twice-logged slope scenario
  acd_draws <- acd_draws %>%
    mutate(habitat = if_else(habitat == "once_logged", "once-logged", habitat)) %>%
    filter(as.character(slope_factor) == as.character(twice_logged_slope_trajectory)) %>%
    rename(functional_habitat = habitat) %>%
    select(-start_age)

  if (!is.null(max_draws)) {
    keep <- sort(unique(acd_draws$draw))[seq_len(min(max_draws, length(unique(acd_draws$draw))))]
    acd_draws <- acd_draws %>% filter(draw %in% keep)
  }

  acd_draws %>%
    group_split(draw) %>%
    lapply(as.data.table)
}

# =============================================================================
# 4) CORE COMPUTATIONS (per scenario-set, across draws)
# =============================================================================

# Join one schedule table to one draw table to attach ACD by habitat/age.
add_carbon_to_schedule <- function(schedule_dt, draw_dt) {
  setDT(schedule_dt)
  setDT(draw_dt)

  x_primary <- schedule_dt[functional_habitat == "primary"]
  x_other <- schedule_dt[functional_habitat != "primary"]

  joined_primary <- x_primary[
    draw_dt[functional_habitat == "primary"],
    on = .(functional_habitat),
    nomatch = 0
  ]
  if ("i.functionalhabAge" %in% names(joined_primary)) joined_primary[, i.functionalhabAge := NULL]

  joined_other <- x_other[
    draw_dt[functional_habitat != "primary"],
    on = .(functional_habitat, functionalhabAge),
    nomatch = 0
  ]

  rbindlist(list(joined_primary, joined_other), use.names = TRUE)
}

# Convert parcel counts into per-delay parcel allocation within each transition/year.
calculate_delays_per_transition <- function(x) {
  setDT(x)
  transitions <- x[, .(num_transition_delays = uniqueN(harvest_delay)),
                   by = .(index, production_target, original_habitat, habitat)]

  x[transitions, on = .(index, production_target, original_habitat, habitat),
    num_transition_delays := i.num_transition_delays]

  x[, parcels_per_delay := num_parcels / pmax(1, num_transition_delays)]
  x
}

# Aggregate joined parcel-level data to scenario-year ACD/carbon summaries.
scenario_acd_by_year <- function(x) {
  #_______________________________________________
  # Function: scenario_ACD_fun  (original notes; adapted)
  # Purpose:  Calculate aboveground carbon density (ACD) at multiple scales:
  #            (1) Parcel-level (10 km² units)
  #            (2) Habitat-transition by year
  #            (3) Scenario-level by year
  # Inputs:
  #    x — data.table containing columns:
  #         index, production_target, original_habitat, habitat,
  #         true_year, ACD, num_parcels, num_transition_delays,
  #         scenarioName, scenarioStart
  # Outputs:
  #    data.table with total scenario ACD per (index, production_target, true_year)
  #_______________________________________________
  setDT(x)

  # ACD (Mg/ha) scaled to carbon stock per 10 km^2 parcel (1000 ha)
  x[, ACD_per_parcel := ACD * 1000]
  x[, ACD_10km2_stag := ACD_per_parcel * parcels_per_delay]

  habitat_year <- x[, .(
    hab_ACD_year = sum(ACD_10km2_stag, na.rm = TRUE),
    scenarioName = first(scenarioName),
    scenarioStart = first(scenarioStart)
  ), by = .(index, production_target, original_habitat, habitat, true_year)]

  habitat_year[, .(
    scen_ACD_year = sum(hab_ACD_year, na.rm = TRUE),
    scenarioName = first(scenarioName),
    scenarioStart = first(scenarioStart)
  ), by = .(index, production_target, true_year)]
}

# Compute cumulative carbon stock-years per draw across the scenario horizon.
carbon_stock_years_by_draw <- function(scenario_year_dt) {
  #____________________________________________________________
  # Function: carbon_stock_years_fun  (original notes; adapted)
  # Purpose:
  #   Calculate total "carbon stock-years" — the sum of annual aboveground carbon
  #   density (ACD) values across the entire duration of each scenario.
  #
  # Inputs:
  #   x — data.table containing at least:
  #        - production_target : numeric or factor indicating scenario target
  #        - index             : scenario or simulation identifier
  #        - scen_ACD_year     : total ACD (Mg) in that scenario for each year
  #
  # Outputs:
  #   data.table with total summed carbon stock (stock-years) per scenario and
  #   production target.
  #____________________________________________________________
  setDT(scenario_year_dt)
  scenario_year_dt[, .(
    cumaltive_stock_year = sum(scen_ACD_year, na.rm = TRUE),
    scenarioName = first(scenarioName),
    scenarioStart = first(scenarioStart)
  ), by = .(production_target, index)]
}

# Summarise stock-years uncertainty (mean, median, 80% interval) across draws.
summarise_stock_years <- function(carbon_yrs_dt) {
  setDT(carbon_yrs_dt)
  carbon_yrs_dt[, .(
    mean_cum_stock_year = mean(cumaltive_stock_year, na.rm = TRUE),
    lwr_cum_stock_year_95 = quantile(cumaltive_stock_year, probs = 0.025, na.rm = TRUE),
    upr_cum_stock_year_95 = quantile(cumaltive_stock_year, probs = 0.975, na.rm = TRUE),
    lwr_cum_stock_year_80 = quantile(cumaltive_stock_year, probs = 0.10, na.rm = TRUE),
    upr_cum_stock_year_80 = quantile(cumaltive_stock_year, probs = 0.90, na.rm = TRUE),
    lwr_cum_stock_year_50 = quantile(cumaltive_stock_year, probs = 0.25, na.rm = TRUE),
    upr_cum_stock_year_50 = quantile(cumaltive_stock_year, probs = 0.75, na.rm = TRUE),
    scenarioName = first(scenarioName),
    scenarioStart = first(scenarioStart)
  ), by = .(index, production_target)]
}

# Build baseline starting-landscape trajectory for each draw/year.
starting_landscape_acd_by_year <- function(draw_dt, all_start_landscape) {
  setDT(draw_dt)
  sl <- copy(as.data.table(all_start_landscape))

  primary_dt <- draw_dt[functional_habitat == "primary"]
  non_primary_dt <- draw_dt[functional_habitat != "primary"]

  primary_joined <- sl[primary_dt, on = .(functional_habitat), nomatch = 0L]
  non_primary_joined <- sl[non_primary_dt, on = .(functional_habitat, functionalhabAge), nomatch = 0L]

  merged <- rbindlist(list(primary_joined, non_primary_joined), use.names = TRUE, fill = TRUE)
  if ("i.functionalhabAge" %in% names(merged)) merged[, `i.functionalhabAge` := NULL]

  merged[, SL_ACD_year := ACD * 1000 * num_parcels]
  merged[, .(SL_ACD_year = sum(SL_ACD_year, na.rm = TRUE)),
         by = .(scenarioStart, true_year = functionalhabAge)]
}

# Align scenario trajectory with starting-landscape baseline for change calculations.
join_scenario_with_start_acd <- function(scenario_year_dt, sl_year_dt) {
  #____________________________________________________
  # Function: join_scenario_with_start_ACD  (original notes; adapted)
  # Purpose:
  #   Join each scenario data.table with its corresponding starting landscape ACD
  #   (both lists must have the same length and order)
  #
  # Inputs:
  #   scenario_list     : list of data.tables, each containing scenario ACD
  #   start_ACD_list    : list of data.tables, each containing starting landscape ACD
  #   join_cols         : character vector of columns to join by (default: scenarioStart and true_year)
  #
  # Output:
  #   List of data.tables with starting landscape ACD added to each scenario
  #____________________________________________________
  merge(
    scenario_year_dt,
    sl_year_dt,
    by = c("scenarioStart", "true_year"),
    all.x = TRUE,
    sort = FALSE
  ) %>% unique()
}

# Convert level trajectories to annual ACD change relative to baseline.
acd_change_dt <- function(dt) {
  # --------------------------------------------------------------------------------
  # Function: ACD_change_function_dt  (original notes; adapted)
  # Purpose:  Calculate annual changes in Aboveground Carbon Density (ACD) and total
  #           carbon for scenario landscapes and the starting landscape.
  # Reason:   Clean, efficient version for posterior draws where error columns are
  #           no longer needed. Works on a data.table and can be applied to lists of
  #           scenarios/draws.
  # --------------------------------------------------------------------------------
  setDT(dt)
  setorder(dt, index, production_target, true_year)

  dt[, `:=`(
    scen_ACD_change = scen_ACD_year - shift(scen_ACD_year),
    SL_ACD_change = SL_ACD_year - shift(SL_ACD_year)
  ), by = .(index, production_target)]

  dt
}

# Convert annual ACD change to CO2-equivalent flux units for SCC monetization.
flux_conversion <- function(dt) {
  #________________________________________________________
  # Function: flux_conversion_function  (original notes; adapted)
  # Purpose: Convert annual ACD (Aboveground Carbon Density) changes
  #          into corresponding carbon fluxes (CO2-equivalent)
  #________________________________________________________
  mw <- 44.01 / 12.01
  dt[, `:=`(
    scen_flux_ACD = scen_ACD_change * mw,
    SL_flux_ACD = SL_ACD_change * mw
  )]
  dt
}

# Prepare SCC lookup tables by discount rate for fast joins.
prepare_scc_tables <- function(socialDR) {
  # Prepare Social Cost of Carbon (SCC) tables
  # SCC datasets contain discounted values of the social cost of carbon (scc_discounted)
  # for different discount rates (2%, 4%, 6%) and years.
  # We'll rename 'year' -> 'true_year' for clean joins later.
  list(
    `2%` = socialDR %>% filter(discount_rate == "2%") %>% rename(true_year = year),
    `4%` = socialDR %>% filter(discount_rate == "4%") %>% rename(true_year = year),
    `6%` = socialDR %>% filter(discount_rate == "6%") %>% rename(true_year = year)
  )
}

# Apply SCC lookup to annual fluxes and compute total discounted SCC impact.
total_scc_impact <- function(dt, scc_dt) {
  #_______________________________________
  # Define core function: socialDR_fun()  (original notes; adapted)
  #_______________________________________
  # This function merges scenario fluxes with the chosen SCC table (SDR),
  # then applies the SCC discounting to compute annual and total monetized
  # carbon impacts (for aboveground carbon).
  # Note: kept structurally aligned with current script
  dt %>%
    left_join(scc_dt, by = "true_year") %>%
    mutate(annual_carbon_impact_ACD = (scen_flux_ACD - SL_flux_ACD) * scc_discounted) %>%
    group_by(index, production_target, scenarioStart, scenarioName) %>%
    summarise(TOTcarbon_ACD_impact = sum(annual_carbon_impact_ACD, na.rm = TRUE), .groups = "drop")
}

# Summarise SCC impacts across draws (mean/median and uncertainty intervals).
summarise_totcarbon_draws <- function(df) {
  #----------------------------------------------
  # Function: summarize_totcarbon_draws  (original notes; adapted)
  # Purpose: Summarize TOTcarbon_ACD_impact across posterior draws
  #----------------------------------------------
  df %>%
    group_by(index, production_target, scenarioStart, scenarioName) %>%
    summarise(
      TOTcarbon_ACD_mean = mean(TOTcarbon_ACD_impact, na.rm = TRUE),
      # 95% and 80% credible intervals
      TOTcarbon_ACD_lwr95 = quantile(TOTcarbon_ACD_impact, 0.025, na.rm = TRUE),
      TOTcarbon_ACD_upr95 = quantile(TOTcarbon_ACD_impact, 0.975, na.rm = TRUE),
      TOTcarbon_ACD_lwr80 = quantile(TOTcarbon_ACD_impact, 0.10, na.rm = TRUE),
      TOTcarbon_ACD_upr80 = quantile(TOTcarbon_ACD_impact, 0.90, na.rm = TRUE),
      .groups = "drop"
    )
}

# End-to-end runner for one scenario set across all selected draws.
run_one_scenario_set <- function(
    scenario_set_dt,
    draw_list,
    all_start_landscape,
    scc_tables,
    scenario_composition,
    limit_indices_per_set = NULL
) {
  setDT(scenario_set_dt)

  if (!is.null(limit_indices_per_set)) {
    keep_idx <- unique(scenario_set_dt$index)[seq_len(min(limit_indices_per_set, length(unique(scenario_set_dt$index))))]
    scenario_set_dt <- scenario_set_dt[index %in% keep_idx]
  }

  # 1) Join carbon to schedule and compute scenario ACD by year for each draw
  scenario_by_draw <- lapply(draw_list, function(d) add_carbon_to_schedule(scenario_set_dt, d))
  scenario_by_draw <- lapply(scenario_by_draw, calculate_delays_per_transition)
  scenario_year_list <- lapply(scenario_by_draw, scenario_acd_by_year)

  # 2) Carbon stock-years (scenario-year ACD summed across years)
  carbon_yrs_list <- lapply(scenario_year_list, carbon_stock_years_by_draw)
  carbon_yrs_dt <- rbindlist(carbon_yrs_list, use.names = TRUE, fill = TRUE, idcol = "draw")
  carbon_stock_years_summary <- summarise_stock_years(carbon_yrs_dt)

  # 3) Starting landscape ACD trajectories per draw
  sl_year_list <- lapply(draw_list, starting_landscape_acd_by_year, all_start_landscape = all_start_landscape)

  # 4) Join scenario and starting landscape, compute fluxes
  joined_list <- Map(join_scenario_with_start_acd, scenario_year_list, sl_year_list)
  joined_list <- lapply(joined_list, acd_change_dt)
  joined_list <- lapply(joined_list, flux_conversion)

  # 5) Monetise fluxes under each discount rate
  impacts <- lapply(names(scc_tables), function(dr) {
    scc_dt <- scc_tables[[dr]]
    impact_by_draw <- lapply(joined_list, total_scc_impact, scc_dt = scc_dt)
    impact_dt <- rbindlist(impact_by_draw, use.names = TRUE, fill = TRUE, idcol = "draw")
    summary_dt <- summarise_totcarbon_draws(impact_dt) %>% mutate(discount_rate = dr)
    summary_dt
  })

  posterior_summary_comb <- bind_rows(impacts) %>%
    left_join(carbon_stock_years_summary, by = c("index", "production_target", "scenarioStart", "scenarioName"))

  # optional: attach scenario composition for downstream plotting convenience
  posterior_summary_comb <- posterior_summary_comb %>%
    left_join(
      scenario_composition %>%
        select(
          scenarioName, index, production_target,
          comp_sum, comp_html,
          propOG, propPlant, prop1L, prop2L, propDeforested, propRestored,
          propEucalyptus, propAlbizia,
          remainingOG, remaining1L, remaining2L,
          hexP, hex1L, hex2L
        ) %>%
        distinct(),
      by = c("scenarioName", "index", "production_target")
    )

  list(
    posterior_summary = posterior_summary_comb,
    carbon_stock_years_summary = carbon_stock_years_summary
  )
}

# =============================================================================
# 5) RUN
# =============================================================================

fixed_params <- time_it("Load fixed scenario params", load_fixed_scenario_params())
all_start_landscape <- time_it(
  "Expand starting landscape (0–60y)",
  expand_start_landscape_through_time(fixed_params$all_start_landscape, years = 0:60)
)

socialDR <- time_it("Load SCC discount rates (03_scc output)", load_social_discount_rates())
scc_tables <- time_it("Prepare SCC tables", prepare_scc_tables(socialDR))

scenario_composition <- time_it("Load scenario composition (MasterAllScenarios.rds)", load_scenario_composition())

scenarios <- NULL
if (isTRUE(USE_CACHE) && file.exists(cache_scenarios_prepped)) {
  cached <- readRDS(cache_scenarios_prepped)
  if (is.list(cached) && identical(cached$version, cache_version) && !is.null(cached$scenarios)) {
    scenarios <- cached$scenarios
    log_line("Loaded cached prepared scenarios: ", cache_scenarios_prepped)
  }
}

if (is.null(scenarios)) {
  scenarios_raw <- time_it("Load scenario schedules (with harvest delays)", load_scenario_schedules())
  scenarios <- time_it("Prepare scenario schedules (apply harvest-delay rules)", prepare_scenarios(scenarios_raw))
  rm(scenarios_raw)

  if (isTRUE(USE_CACHE)) {
    saveRDS(
      list(
        version = cache_version,
        created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        scenarios = scenarios
      ),
      cache_scenarios_prepped
    )
    log_line("Wrote cached prepared scenarios: ", cache_scenarios_prepped)
  }
}

# -----------------------------------------------------------------------------
# CHECKPOINT: inspect prepared schedules BEFORE any carbon joins
# -----------------------------------------------------------------------------
# These are created explicitly so you can inspect rule application in the
# environment before draw joins / carbon propagation starts.
scenarios_prepped_for_review <- scenarios

# Combined long table with scenario-set id for quick filtering / View().
scenarios_prepped_combined <- rbindlist(
  Map(
    function(dt, id) {
      x <- as.data.table(copy(dt))
      x[, scenario_set_id := id]
      x
    },
    scenarios_prepped_for_review,
    seq_along(scenarios_prepped_for_review)
  ),
  use.names = TRUE,
  fill = TRUE
)

# Helpful console guidance when running interactively.
if (interactive()) {
  message("Prepared schedules available for review:")
  message(" - scenarios_prepped_for_review (list by scenario set)")
  message(" - scenarios_prepped_combined (single table; includes scenario_set_id)")
  message("Example checks:")
  message("  View(scenarios_prepped_combined)")
  message("  scenarios_prepped_combined[original_habitat == habitat & harvest_delay != 0]")
  message("  scenarios_prepped_combined[original_habitat == 'once-logged' & habitat == 'restored']")
}

acd_draws <- time_it("Load ACD draws (02_draws output)", load_acd_draws())
trajectory_values <- twice_logged_slope_trajectories
if (is.null(trajectory_values) || length(trajectory_values) == 0) {
  trajectory_values <- sort(unique(as.character(acd_draws$slope_factor)))
}
trajectory_values <- unique(as.character(trajectory_values))
log_line("Running trajectories: ", paste(trajectory_values, collapse = ", "))

if (!is.null(scenario_set_idxs)) {
  scenario_set_idxs <- scenario_set_idxs[scenario_set_idxs >= 1 & scenario_set_idxs <= length(scenarios)]
  scenarios <- scenarios[scenario_set_idxs]
}

posterior_summary_all_by_trajectory <- list()
posterior_summary_combined_by_trajectory <- list()
final_performance_carbon_by_trajectory <- list()
carbon_stock_years_out_by_trajectory <- list()

for (traj in trajectory_values) {
  traj_suffix <- trajectory_file_suffix(traj)
  log_line("----- 2L trajectory: ", traj, " -----")

  draw_list <- time_it(
    paste0("Prepare per-draw lookup tables (", traj_suffix, ")"),
    prepare_draw_list(acd_draws, traj, max_draws = max_draws)
  )

  posterior_summary_all <- vector("list", length(scenarios))
  names(posterior_summary_all) <- paste0("scenario_set_", seq_along(scenarios))
  log_line("Running ", length(scenarios), " scenario set(s) for ", traj_suffix, " ...")

  for (i in seq_along(scenarios)) {
    log_line("Scenario set ", i, " / ", length(scenarios), " [", traj_suffix, "]")
    res <- run_one_scenario_set(
      scenario_set_dt = scenarios[[i]],
      draw_list = draw_list,
      all_start_landscape = all_start_landscape,
      scc_tables = scc_tables,
      scenario_composition = scenario_composition,
      limit_indices_per_set = limit_indices_per_set
    )
    posterior_summary_all[[i]] <- res$posterior_summary
  }

  # Combine scenario-set outputs into one table for export + plotting
  posterior_summary_combined <- rbindlist(posterior_summary_all, use.names = TRUE, fill = TRUE) %>%
    mutate(twice_logged_slope_trajectory = as.character(traj))

  # Final performance exports (for consolidated cross-outcome figures)
  final_performance_carbon <- posterior_summary_combined %>%
    select(
      index, production_target, scenarioName, scenarioStart, discount_rate,
      TOTcarbon_ACD_mean, TOTcarbon_ACD_lwr80, TOTcarbon_ACD_upr80, TOTcarbon_ACD_lwr95, TOTcarbon_ACD_upr95,
      twice_logged_slope_trajectory
    ) %>%
    distinct() %>%
    mutate(outcome = "carbon")

  carbon_stock_years_out <- posterior_summary_combined %>%
    select(
      index, production_target, scenarioName, scenarioStart,
      mean_cum_stock_year, lwr_cum_stock_year_80, upr_cum_stock_year_80, lwr_cum_stock_year_95, upr_cum_stock_year_95,
      twice_logged_slope_trajectory
    ) %>%
    distinct() %>%
    mutate(outcome = "carbon")

  # Per-trajectory saves (explicitly labelled by trajectory)
  saveRDS(posterior_summary_all, file.path(rds_dir, paste0("carbon_outcomes__", traj_suffix, ".rds")))
  write.csv(
    as.data.frame(posterior_summary_combined),
    file.path(tab_dir, paste0("carbon_outcomes_combined__", traj_suffix, ".csv")),
    row.names = FALSE
  )
  saveRDS(
    final_performance_carbon,
    file.path(rds_dir, paste0("final_performance__carbon__with_uncertainty__", traj_suffix, ".rds"))
  )
  write.csv(
    as.data.frame(final_performance_carbon),
    file.path(tab_dir, paste0("final_performance__carbon__with_uncertainty__", traj_suffix, ".csv")),
    row.names = FALSE
  )
  saveRDS(
    carbon_stock_years_out,
    file.path(rds_dir, paste0("final_performance__carbon_stock_years__with_uncertainty__", traj_suffix, ".rds"))
  )
  write.csv(
    as.data.frame(carbon_stock_years_out),
    file.path(tab_dir, paste0("final_performance__carbon_stock_years__with_uncertainty__", traj_suffix, ".csv")),
    row.names = FALSE
  )

  posterior_summary_all_by_trajectory[[traj_suffix]] <- posterior_summary_all
  posterior_summary_combined_by_trajectory[[traj_suffix]] <- posterior_summary_combined
  final_performance_carbon_by_trajectory[[traj_suffix]] <- final_performance_carbon
  carbon_stock_years_out_by_trajectory[[traj_suffix]] <- carbon_stock_years_out
}

rm(acd_draws)

# Keep backwards-compatible objects for downstream plotting (first trajectory).
first_traj <- names(posterior_summary_all_by_trajectory)[1]
posterior_summary_all <- posterior_summary_all_by_trajectory[[first_traj]]
posterior_summary_combined <- posterior_summary_combined_by_trajectory[[first_traj]]
final_performance_carbon <- final_performance_carbon_by_trajectory[[first_traj]]
carbon_stock_years_out <- carbon_stock_years_out_by_trajectory[[first_traj]]

# =============================================================================
# 6) SAVE OUTPUTS
# =============================================================================

# Combined outputs across all trajectories
posterior_summary_combined_all_trajectories <- bind_rows(posterior_summary_combined_by_trajectory)
final_performance_carbon_all_trajectories <- bind_rows(final_performance_carbon_by_trajectory)
carbon_stock_years_out_all_trajectories <- bind_rows(carbon_stock_years_out_by_trajectory)

saveRDS(posterior_summary_all_by_trajectory, file.path(rds_dir, "carbon_outcomes__all_trajectories.rds"))
writeLines(capture.output(sessionInfo()), con = file.path(out_root, "sessionInfo.txt"))

# Save a single combined CSV for easy inspection / plotting elsewhere
write.csv(
  as.data.frame(posterior_summary_combined_all_trajectories),
  file.path(tab_dir, "carbon_outcomes_combined__all_trajectories.csv"),
  row.names = FALSE
)

saveRDS(
  final_performance_carbon_all_trajectories,
  file.path(rds_dir, "final_performance__carbon__with_uncertainty__all_trajectories.rds")
)
write.csv(
  as.data.frame(final_performance_carbon_all_trajectories),
  file.path(tab_dir, "final_performance__carbon__with_uncertainty__all_trajectories.csv"),
  row.names = FALSE
)

saveRDS(
  carbon_stock_years_out_all_trajectories,
  file.path(rds_dir, "final_performance__carbon_stock_years__with_uncertainty__all_trajectories.rds")
)
write.csv(
  as.data.frame(carbon_stock_years_out_all_trajectories),
  file.path(tab_dir, "final_performance__carbon_stock_years__with_uncertainty__all_trajectories.csv"),
  row.names = FALSE
)

# =============================================================================
# 7) PLOTS (clean summary figures)
# =============================================================================

# Master-style plots (what you used in the manuscript figure code):
# - colour by amount of old-growth (propOG)
# - shape indicates whether there is any plantation in the scenario (propPlant > 0)
# - points + intervals share the SAME jitter so they stay connected
#
# Note: `posterior_summary_*` includes `propOG` + `propPlant` from `scenario_composition`.

# Static SCC figure for manuscript-style comparison across scenarios.
plot_scc_master <- function(df) {
  x <- df %>%
    mutate(
      has_plantation = if_else(!is.na(propPlant) & propPlant > 0, "Plantation", "No plantation"),
      scen_col = case_when(
        scenarioStart %in% c("all_primary", "primary_deforested") ~ hexP,
        scenarioStart %in% c("mostly_1L", "mostly_1L_deforested") ~ hex1L,
        scenarioStart %in% c("mostly_2L", "mostly_2L_deforested") ~ hex2L,
        TRUE ~ "#9E9E9E"
      )
    )

  # Use the SAME jitter for points and error bars so they remain connected
  pos <- position_jitter(width = 0.05, height = 0)

  x %>%
    ggplot(aes(
      x = production_target,
      y = TOTcarbon_ACD_mean / 1e9,
      colour = scen_col,
      shape = has_plantation
    )) +
    geom_errorbar(
      aes(ymin = TOTcarbon_ACD_lwr80 / 1e9, ymax = TOTcarbon_ACD_upr80 / 1e9),
      position = pos,
      width = 0.00,
      alpha = 0.45,
      linewidth = 0.6
    ) +
    geom_point(
      position = pos,
      alpha = 0.55,
      stroke = 0.7,
      size = 2.0
    ) +
    scale_colour_identity() +
    scale_shape_manual(values = c("No plantation" = 19, "Plantation" = 2)) +
    labs(
      x = "Production target",
      y = "Total carbon impact (billion; ACD only)",
      title = "SCC-discounted carbon impacts across scenarios (80% intervals)"
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    facet_grid(scenarioName ~ discount_rate) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      strip.text.y = element_blank()
    )
}

# Static stock-years figure (uncertainty-aware) across scenarios.
plot_stock_years_master <- function(df) {
  stock_df <- df %>%
    distinct(
      scenarioName, index, production_target, scenarioStart,
      propOG, propPlant, hexP, hex1L, hex2L,
      mean_cum_stock_year, lwr_cum_stock_year_80, upr_cum_stock_year_80
    ) %>%
    mutate(
      has_plantation = if_else(!is.na(propPlant) & propPlant > 0, "Plantation", "No plantation"),
      scen_col = case_when(
        scenarioStart %in% c("all_primary", "primary_deforested") ~ hexP,
        scenarioStart %in% c("mostly_1L", "mostly_1L_deforested") ~ hex1L,
        scenarioStart %in% c("mostly_2L", "mostly_2L_deforested") ~ hex2L,
        TRUE ~ "#9E9E9E"
      )
    )

  # Use the SAME jitter for points and error bars so they remain connected
  pos <- position_jitter(width = 0.015, height = 0)

  stock_df %>%
    ggplot(aes(
      x = production_target,
      y = mean_cum_stock_year,
      colour = scen_col,
      shape = has_plantation
    )) +
    geom_errorbar(
      aes(ymin = lwr_cum_stock_year_80, ymax = upr_cum_stock_year_80),
      position = pos,
      width = 0.015,
      linewidth = 0.6,
      alpha = 0.5
    ) +
    geom_point(
      position = pos,
      alpha = 0.55,
      stroke = 0.7,
      size = 2.0
    ) +
    scale_colour_identity() +
    scale_shape_manual(values = c("No plantation" = 19, "Plantation" = 2)) +
    labs(
      x = "Production target",
      y = "Carbon stock-years (Mg C)",
      title = "Carbon stock-years across scenarios (80% intervals)"
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    facet_wrap(~scenarioName, ncol = 4) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      strip.text = element_blank()
    )
}

# Sanitize labels so filenames are valid/stable across OSs.
safe_filename <- function(x, max_len = 120) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (nchar(x) > max_len) x <- substr(x, 1, max_len)
  if (!nzchar(x)) x <- "unnamed"
  x
}

# -----------------------------------------------------------------------------
# Interactive (hover + lasso/box select) HTML versions of the master plots
# -----------------------------------------------------------------------------
if (!requireNamespace("plotly", quietly = TRUE)) {
  stop("Package `plotly` is required for interactive exploration. Please run: install.packages('plotly')")
}
if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
  stop("Package `htmlwidgets` is required for saving interactive HTML. Please run: install.packages('htmlwidgets')")
}

# Save plotly/htmlwidgets output robustly (short lib path for Windows/OneDrive).
save_widget_in_dir <- function(widget, file_path, selfcontained = FALSE) {
  out_dir <- dirname(file_path)
  out_file <- basename(file_path)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(out_dir)

  htmlwidgets::saveWidget(
    widget,
    file = out_file,
    # On Windows/OneDrive we can hit path-length limits if saveWidget creates a
    # long "<htmlname>_files/..." dependency folder. Using a fixed short libdir
    # is much more robust.
    selfcontained = isTRUE(selfcontained),
    libdir = "lib"
  )
  invisible(file_path)
}

# Interactive SCC view (hover + lasso) to inspect individual scenarios.
plotly_scc_master <- function(df) {
  x <- df %>%
    mutate(
      has_plantation = if_else(!is.na(propPlant) & propPlant > 0, "Plantation", "No plantation"),
      scen_col = case_when(
        scenarioStart %in% c("all_primary", "primary_deforested") ~ hexP,
        scenarioStart %in% c("mostly_1L", "mostly_1L_deforested") ~ hex1L,
        scenarioStart %in% c("mostly_2L", "mostly_2L_deforested") ~ hex2L,
        TRUE ~ "#9E9E9E"
      ),
      hover = paste0(
        "<b>scenarioName</b>: ", scenarioName,
        "<br><b>discount_rate</b>: ", discount_rate,
        "<br><b>index</b>: ", index,
        "<br><b>production_target</b>: ", production_target,
        "<br><b>scenarioStart</b>: ", scenarioStart,
        "<br><b>composition sum</b>: ", signif(comp_sum, 4),
        "<br>", comp_html,
        "<br><b>TOTcarbon_ACD_mean (bn)</b>: ", signif(TOTcarbon_ACD_mean / 1e9, 4),
        "<br><b>80% CI (bn)</b>: [", signif(TOTcarbon_ACD_lwr80 / 1e9, 4), ", ", signif(TOTcarbon_ACD_upr80 / 1e9, 4), "]"
      )
    )

  drs <- sort(unique(x$discount_rate))
  plots <- lapply(drs, function(dr) {
    xd <- x %>% filter(discount_rate == dr)
    plotly::plot_ly(
      data = xd,
      x = ~production_target,
      y = ~(TOTcarbon_ACD_mean / 1e9),
      type = "scatter",
      mode = "markers",
      text = ~hover,
      hoverinfo = "text",
      marker = list(
        color = ~scen_col,
        symbol = ~if_else(has_plantation == "Plantation", "x-open", "circle"),
        size = 7,
        line = list(width = 0.6)
      ),
      error_y = list(
        type = "data",
        symmetric = FALSE,
        array = ~(TOTcarbon_ACD_upr80 - TOTcarbon_ACD_mean) / 1e9,
        arrayminus = ~(TOTcarbon_ACD_mean - TOTcarbon_ACD_lwr80) / 1e9,
        thickness = 1
      )
    ) %>%
      plotly::layout(
        title = list(text = paste0("Discount rate: ", dr)),
        xaxis = list(title = "Production target", range = c(0, 1)),
        yaxis = list(title = "Total carbon impact (billion; ACD only)")
      )
  })

  plotly::subplot(plots, nrows = 1, shareX = TRUE, titleX = TRUE, titleY = TRUE) %>%
    plotly::layout(dragmode = "lasso", showlegend = FALSE)
}

# Interactive stock-years view (hover + lasso) for exploratory QA.
plotly_stock_years_master <- function(df) {
  x <- df %>%
    distinct(
      scenarioName, index, production_target, scenarioStart,
      comp_sum, comp_html,
      propOG, propPlant, prop1L, prop2L, propDeforested, propRestored, hexP, hex1L, hex2L,
      mean_cum_stock_year, lwr_cum_stock_year_80, upr_cum_stock_year_80
    ) %>%
    mutate(
      has_plantation = if_else(!is.na(propPlant) & propPlant > 0, "Plantation", "No plantation"),
      scen_col = case_when(
        scenarioStart %in% c("all_primary", "primary_deforested") ~ hexP,
        scenarioStart %in% c("mostly_1L", "mostly_1L_deforested") ~ hex1L,
        scenarioStart %in% c("mostly_2L", "mostly_2L_deforested") ~ hex2L,
        TRUE ~ "#9E9E9E"
      ),
      hover = paste0(
        "<b>scenarioName</b>: ", scenarioName,
        "<br><b>index</b>: ", index,
        "<br><b>production_target</b>: ", production_target,
        "<br><b>scenarioStart</b>: ", scenarioStart,
        "<br><b>composition sum</b>: ", signif(comp_sum, 4),
        "<br>", comp_html,
        "<br><b>stock_years mean</b>: ", signif(mean_cum_stock_year, 5),
        "<br><b>80% CI</b>: [", signif(lwr_cum_stock_year_80, 5), ", ", signif(upr_cum_stock_year_80, 5), "]"
      )
    )

  plotly::plot_ly(
    data = x,
    x = ~production_target,
    y = ~mean_cum_stock_year,
    type = "scatter",
    mode = "markers",
    text = ~hover,
    hoverinfo = "text",
    marker = list(
      color = ~scen_col,
      symbol = ~if_else(has_plantation == "Plantation", "x-open", "circle"),
      size = 7,
      line = list(width = 0.6)
    ),
    error_y = list(
      type = "data",
      symmetric = FALSE,
      array = ~(upr_cum_stock_year_80 - mean_cum_stock_year),
      arrayminus = ~(mean_cum_stock_year - lwr_cum_stock_year_80),
      thickness = 1
    )
  ) %>%
    plotly::layout(
      dragmode = "lasso",
      xaxis = list(title = "Production target", range = c(0, 1)),
      yaxis = list(title = "Carbon stock-years (Mg C)")
    )
}

# Save plots per scenario set, naming files by scenarioName (unique to each set).
for (k in seq_along(posterior_summary_all)) {
  df_k <- as.data.frame(posterior_summary_all[[k]])
  scen_names <- unique(df_k$scenarioName)
  scen_label <- if (length(scen_names) == 1) scen_names else paste0("multiple_", k)
  scen_label <- sub("\\.csv$", "", scen_label, ignore.case = TRUE)

  p1 <- plot_scc_master(df_k)
  ggsave(
    file.path(fig_dir, paste0("scc_master__", safe_filename(scen_label), ".png")),
    p1, width = 10, height = 6, units = "in", dpi = 220
  )

  p2 <- plot_stock_years_master(df_k)
  ggsave(
    file.path(fig_dir, paste0("stock_years_master__", safe_filename(scen_label), ".png")),
    p2, width = 12, height = 8, units = "in", dpi = 220
  )

  # Interactive HTML (click/hover/select)
  p1_html <- plotly_scc_master(df_k)
  save_widget_in_dir(
    widget = p1_html,
    file.path(fig_dir, paste0("scc_i__", safe_filename(scen_label), ".html")),
    selfcontained = FALSE
  )

  p2_html <- plotly_stock_years_master(df_k)
  save_widget_in_dir(
    widget = p2_html,
    file.path(fig_dir, paste0("sy_i__", safe_filename(scen_label), ".html")),
    selfcontained = FALSE
  )
}

# Also save combined master plots (useful if you run multiple scenario sets together).
p1_all <- plot_scc_master(as.data.frame(posterior_summary_combined))
ggsave(file.path(fig_dir, "scc_master__ALL.png"), p1_all, width = 12, height = 8, units = "in", dpi = 220)

p2_all <- plot_stock_years_master(as.data.frame(posterior_summary_combined))
ggsave(file.path(fig_dir, "stock_years_master__ALL.png"), p2_all, width = 12, height = 8, units = "in", dpi = 220)

# Interactive combined views
save_widget_in_dir(
  widget = plotly_scc_master(as.data.frame(posterior_summary_combined)),
  file.path(fig_dir, "scc_i__ALL.html"),
  selfcontained = FALSE
)
save_widget_in_dir(
  widget = plotly_stock_years_master(as.data.frame(posterior_summary_combined)),
  file.path(fig_dir, "sy_i__ALL.html"),
  selfcontained = FALSE
)

log_line("Saved: ", file.path(rds_dir, "carbon_outcomes__all_trajectories.rds"))
log_line("Done.")

