# =============================================================================
# Nature Revision 2 — Step 04 (TEMP refactor): propagate carbon through scenarios
# =============================================================================
#
# Purpose
# - Clean, modular rewrite of `04_propage_carbon_thru_scenarios.R` for iterative development.
# - Produces the same *type* of outputs, but is easier to test and reason about.
#
# Key features
# - Deterministic outputs under: Outputs/Nature_Revision_Outputs/NR2/current/04_temp/
# - TEST_MODE: run a small subset of scenario-sets (and optionally fewer posterior draws)
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

# =============================================================================
# 0) OUTPUTS + RUN SETTINGS
# =============================================================================

source(file.path("Scripts", "Nature_Revision_2", "_config.R"))
paths <- nr2_step_paths("04_temp")
nr2_ensure_dirs(paths)

# Base folder for this temp step (shared across runs)
out_root <- paths$root

# ---------------------------------
# Per-run labelled output subfolder
# ---------------------------------
RUN_LABEL <- "run_01"  # change this each time you want separate figures/outputs
if (identical(RUN_LABEL, "") || is.null(RUN_LABEL)) {
  RUN_LABEL <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
}

run_root <- file.path(out_root, "runs", RUN_LABEL)
fig_dir <- file.path(run_root, "figures")
tab_dir <- file.path(run_root, "tables")
rds_dir <- file.path(run_root, "rds")

dir.create(run_root, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(run_root, "run_log.txt")
log_line <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(txt, "\n", file = log_path, append = TRUE)
  message(txt)
}

# -----------------------------
# Main knobs for fast iteration
# -----------------------------
TEST_MODE <- TRUE

# Cache the slowest step (scenario schedules prep) under the temp output folder.
# This makes repeated iterations much faster.
USE_CACHE <- TRUE

# Scenario-set selection
# - A "scenario set" = one element of `Inputs/MasterAllScenarios_withHarvestDelays.rds` (a list).
# - In full runs, you typically process ALL sets.
scenario_set_idxs <- c(1,5,8)       # e.g. c(1, 5, 9). If NULL: auto (TEST_MODE) or all (full).

# Within-set filtering (optional): keep only a few scenario `index` values per set
limit_indices_per_set <- NULL    # e.g. 5 (keep first 5 unique index values)

# Posterior draw subsetting (optional): use fewer draws for smoke tests
max_draws <- 50                # e.g. 50 (keep first 50 draws)

# Twice-logged recovery trajectory choice (matches current script)
twice_logged_slope_trajectory <- 1

set.seed(123)

if (isTRUE(TEST_MODE)) {
  # Very small default smoke-test (fast iteration).
  # Increase these gradually as you gain confidence.
  if (is.null(scenario_set_idxs)) scenario_set_idxs <- 1
 # if (is.null(limit_indices_per_set)) limit_indices_per_set <- 2
  if (is.null(max_draws)) max_draws <- 10
}

log_line("TEST_MODE = ", TEST_MODE)
log_line("RUN_LABEL = ", RUN_LABEL)
log_line("scenario_set_idxs = ", paste(scenario_set_idxs %||% "ALL", collapse = ", "))
log_line("limit_indices_per_set = ", limit_indices_per_set %||% "NULL")
log_line("max_draws = ", max_draws %||% "NULL")
log_line("twice_logged_slope_trajectory = ", twice_logged_slope_trajectory)
log_line("nr2_out_root = ", normalizePath(nr2_out_root, winslash = "/", mustWork = FALSE))
log_line("USE_CACHE = ", USE_CACHE)

cache_dir <- file.path(out_root, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
cache_scenarios_prepped <- file.path(cache_dir, "scenarios_prepped.rds")
cache_version <- "04_temp_scenarios_v1"

time_it <- function(label, expr) {
  t0 <- Sys.time()
  log_line(label, " ...")
  out <- force(expr)
  dt <- difftime(Sys.time(), t0, units = "secs")
  log_line(label, " done (", sprintf("%.1f", as.numeric(dt)), "s)")
  out
}

# =============================================================================
# 1) INPUTS
# =============================================================================

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

load_social_discount_rates <- function(path = file.path(nr2_out_root, "03_scc", "tables", "scc_dr_2_4_6.csv")) {
  stopifnot(file.exists(path))
  read.csv(path) %>%
    group_by(discount_rate) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(scc_discounted_ratio = scc_discounted / first(scc_discounted)) %>%
    ungroup()
}

load_scenario_composition <- function(path = file.path("Inputs", "MasterAllScenarios.rds")) {
  stopifnot(file.exists(path))
  scenarios2 <- readRDS(path)
  out <- rbindlist(scenarios2, use.names = TRUE) %>%
    group_by(scenarioName, index, production_target) %>%
    mutate(
      propOG = sum(num_parcels[habitat == "primary"]) / 1000,
      propPlant = sum(num_parcels[habitat %in% c(
        "eucalyptus_current", "albizia_current",
        "albizia_future", "eucalyptus_future"
      )]) / 1000
    )
  rm(scenarios2)
  out
}

load_scenario_schedules <- function(path = file.path("Inputs", "MasterAllScenarios_withHarvestDelays.rds")) {
  stopifnot(file.exists(path))
  readRDS(path)
}

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
# b. 1L -> 1L if parcel starts as once-logged and stays 1L, assume logging happened 15 yrs befor scnanario start
# c. 1L -> 2L if parcel starts as once logged before re-harvest, we assume 1L recovery
#    (already explicit in scenarios) that resets to yr 0 once-logged values.
#
# for b. we use the function update_1L_functional_habitat to ensure that scenarios beginning
# the scenario as once-logged were harvest 15yrs before
#
# if parcel starts as twice-logged, use the 15yr after 1L curves - with different slope factors,
# tracking the different slope assumption rates for twice-logged forest
# if parcel starts as primary and goes to twice-logged, use the 30yrafter 1L, with different slope factors
# ---------------------------------------------------------------------------

filter_original_eq_habitat <- function(df_list) {
  df_list %>% map(~ .x %>% filter(!(original_habitat == habitat & harvest_delay != 0)))
}

update_1L_functional_habitat <- function(df_list) {
  # scenarios that start as once-logged are treated as already 15y post-logging at scenario start
  df_list %>% map(~ .x %>% filter(!(original_habitat == "once-logged" & harvest_delay < 15)))
}

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

remove_no_timber_scenarios <- function(df_list) {
  df_list %>% map(~ .x %>% filter(production_target > 0))
}

prepare_scenarios <- function(scenarios) {
  scenarios <- scenarios %>%
    map(~ .x %>% mutate(harvest_delay = as.numeric(stringr::str_extract(harvest_delay, "\\d+"))))

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

prepare_draw_list <- function(acd_draws, twice_logged_slope_trajectory, max_draws = NULL) {
  # Match legacy naming and select the desired twice-logged slope scenario
  acd_draws <- acd_draws %>%
    mutate(habitat = if_else(habitat == "once_logged", "once-logged", habitat)) %>%
    filter(slope_factor == twice_logged_slope_trajectory) %>%
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

calculate_delays_per_transition <- function(x) {
  setDT(x)
  transitions <- x[, .(num_transition_delays = uniqueN(harvest_delay)),
                   by = .(index, production_target, original_habitat, habitat)]

  x[transitions, on = .(index, production_target, original_habitat, habitat),
    num_transition_delays := i.num_transition_delays]

  x[, parcels_per_delay := num_parcels / pmax(1, num_transition_delays)]
  x
}

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

summarise_totcarbon_draws <- function(df) {
  #----------------------------------------------
  # Function: summarize_totcarbon_draws  (original notes; adapted)
  # Purpose: Summarize TOTcarbon_ACD_impact across posterior draws
  #----------------------------------------------
  df %>%
    group_by(index, production_target, scenarioStart, scenarioName) %>%
    summarise(
      TOTcarbon_ACD_mean = mean(TOTcarbon_ACD_impact, na.rm = TRUE),
      # NOTE: these quantiles match the existing script even though the names suggest 95%
      TOTcarbon_ACD_lwr95 = quantile(TOTcarbon_ACD_impact, 0.25, na.rm = TRUE),
      TOTcarbon_ACD_upr95 = quantile(TOTcarbon_ACD_impact, 0.75, na.rm = TRUE),
      TOTcarbon_ACD_lwr80 = quantile(TOTcarbon_ACD_impact, 0.10, na.rm = TRUE),
      TOTcarbon_ACD_upr80 = quantile(TOTcarbon_ACD_impact, 0.90, na.rm = TRUE),
      .groups = "drop"
    )
}

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
      scenario_composition %>% select(scenarioName, index, production_target, propOG, propPlant) %>% distinct(),
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

acd_draws <- time_it("Load ACD draws (02_draws output)", load_acd_draws())
draw_list <- time_it(
  "Prepare per-draw lookup tables",
  prepare_draw_list(acd_draws, twice_logged_slope_trajectory, max_draws = max_draws)
)
rm(acd_draws)

if (!is.null(scenario_set_idxs)) {
  scenario_set_idxs <- scenario_set_idxs[scenario_set_idxs >= 1 & scenario_set_idxs <= length(scenarios)]
  scenarios <- scenarios[scenario_set_idxs]
}

posterior_summary_all <- vector("list", length(scenarios))
names(posterior_summary_all) <- paste0("scenario_set_", seq_along(scenarios))

log_line("Running ", length(scenarios), " scenario set(s) ...")

for (i in seq_along(scenarios)) {
  log_line("Scenario set ", i, " / ", length(scenarios))
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
posterior_summary_combined <- rbindlist(posterior_summary_all, use.names = TRUE, fill = TRUE)

# =============================================================================
# 6) SAVE OUTPUTS
# =============================================================================

saveRDS(posterior_summary_all, file.path(rds_dir, "carbon_outcomes_temp.rds"))
writeLines(capture.output(sessionInfo()), con = file.path(out_root, "sessionInfo.txt"))

# Save a single combined CSV for easy inspection / plotting elsewhere
write.csv(
  as.data.frame(posterior_summary_combined),
  file.path(tab_dir, "carbon_outcomes_temp_combined.csv"),
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

plot_scc_master <- function(df) {
  x <- df %>%
    mutate(
      has_plantation = if_else(!is.na(propPlant) & propPlant > 0, "Plantation", "No plantation")
    )

  # Use the SAME jitter for points and error bars so they remain connected
  pos <- position_jitter(width = 0.05, height = 0)

  x %>%
    ggplot(aes(
      x = production_target,
      y = TOTcarbon_ACD_mean / 1e9,
      colour = propOG,
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
    scale_colour_gradient(
      low = "#FFE8C2",
      high = "#B30000",
      name = "Proportion OG"
    ) +
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

plot_stock_years_master <- function(df) {
  stock_df <- df %>%
    distinct(
      scenarioName, index, production_target, scenarioStart,
      propOG, propPlant,
      mean_cum_stock_year, lwr_cum_stock_year_80, upr_cum_stock_year_80
    ) %>%
    mutate(
      has_plantation = if_else(!is.na(propPlant) & propPlant > 0, "Plantation", "No plantation")
    )

  # Use the SAME jitter for points and error bars so they remain connected
  pos <- position_jitter(width = 0.015, height = 0)

  stock_df %>%
    ggplot(aes(
      x = production_target,
      y = mean_cum_stock_year,
      colour = propOG,
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
    scale_colour_gradient(
      name = "Proportion old growth",
      low = "#FFE8C2",
      high = "#B30000",
    ) +
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

safe_filename <- function(x, max_len = 120) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (nchar(x) > max_len) x <- substr(x, 1, max_len)
  if (!nzchar(x)) x <- "unnamed"
  x
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
}

# Also save combined master plots (useful if you run multiple scenario sets together).
p1_all <- plot_scc_master(as.data.frame(posterior_summary_combined))
ggsave(file.path(fig_dir, "scc_master__ALL.png"), p1_all, width = 12, height = 8, units = "in", dpi = 220)

p2_all <- plot_stock_years_master(as.data.frame(posterior_summary_combined))
ggsave(file.path(fig_dir, "stock_years_master__ALL.png"), p2_all, width = 12, height = 8, units = "in", dpi = 220)

log_line("Saved: ", file.path(rds_dir, "carbon_outcomes_temp.rds"))
log_line("Done.")

