# ----------------------------------------------------------------------------
# Nature Revision 2 — single-scenario debugger for carbon propagation
#
# When something looks wrong in step 04, I run this to walk one index through each intermediate join without the full wrapper functions.
# Inputs: same intermediate RDS as the main 04 script (paths near the top).
# Outputs: printed diagnostics and optional small tables; not used for production figures directly.
# ----------------------------------------------------------------------------

## =============================================================================
## Nature Revision 2 — Step-through ONE scenario (piece-by-piece; no wrappers)
## =============================================================================
##
## Goal
## - For ONE chosen scenario, run each step of the carbon pipeline explicitly,
##   storing the input/output tables at every stage so you can inspect/plot them
##   in the R environment.
##
## Usage (interactive R, from repo root):
##   source("Scripts/Nature_Revision_2/04_stepthrough_single_scenario_piece_by_piece.R")
##   View(sched_raw)
##   View(sched_prepped)
##   View(sched_carbon)
##   View(habitat_year)
##   View(scen_year)
##   View(sl_year)
##   View(joined)
##   View(changes)
##   View(flux)
##   View(annual_scc)
##

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

source(file.path("Scripts", "Nature_Revision2", "_config.R"))

## -----------------------------------------------------------------------------
## 0) SETTINGS (edit these)
## -----------------------------------------------------------------------------

SCENARIO_NAME <- "mostly_1L_CY_D.csv"
SCENARIO_INDEX <- "mostly_1L_CY_D.csv 312"
SCENARIO_START_EXPECTED <- "mostly_1L"

DRAW_ID <- 1
TWICE_LOGGED_SLOPE_TRAJECTORY <- 1
DISCOUNT_RATE <- "2%"  # one of: "2%", "4%", "6%"

YEARS <- 0:60

## -----------------------------------------------------------------------------
## 1) LOAD ONE SCENARIO SCHEDULE (raw)
## -----------------------------------------------------------------------------

scenarios_path <- file.path("Inputs", "MasterAllScenarios_withHarvestDelays.rds")
stopifnot(file.exists(scenarios_path))
scenarios_raw <- readRDS(scenarios_path)
stopifnot(is.list(scenarios_raw), length(scenarios_raw) >= 1)

scenario_set_names <- vapply(
  scenarios_raw,
  function(x) as.character(unique(as.data.frame(x)$scenarioName)[1]),
  character(1)
)

scenario_set_idx <- which(scenario_set_names == SCENARIO_NAME)[1]
stopifnot(!is.na(scenario_set_idx))

sched_all <- as.data.table(copy(scenarios_raw[[scenario_set_idx]]))

## Filter to the exact scenario (index + scenarioName). If production_target is not unique,
## we pick the first and print what happened.
sched_subset <- sched_all[scenarioName == SCENARIO_NAME & index == SCENARIO_INDEX]
stopifnot(nrow(sched_subset) > 0)

pts <- sort(unique(sched_subset$production_target))
if (length(pts) > 1) {
  message("Multiple production_target values detected for this index. Using the first: ", pts[[1]])
  sched_subset <- sched_subset[production_target == pts[[1]]]
}

sched_raw <- copy(sched_subset)
View(sched_raw)
## Sanity: scenarioStart should match what you think
scen_start_vals <- unique(sched_raw$scenarioStart)
if (length(scen_start_vals) != 1) {
  message("Multiple scenarioStart values present: ", paste(scen_start_vals, collapse = ", "))
}
if (length(scen_start_vals) >= 1 && !identical(as.character(scen_start_vals[[1]]), SCENARIO_START_EXPECTED)) {
  message("WARNING: scenarioStart is '", as.character(scen_start_vals[[1]]),
          "' but you expected '", SCENARIO_START_EXPECTED, "'.")
}

## -----------------------------------------------------------------------------
## 2) PREP schedule like 04_temp.R (explicit; no wrappers)
## -----------------------------------------------------------------------------

sched_prepped <- copy(sched_raw)

## Convert harvest_delay (if it's encoded like "delay_15")
if ("harvest_delay" %in% names(sched_prepped) && !is.numeric(sched_prepped$harvest_delay)) {
  sched_prepped[, harvest_delay := as.numeric(stringr::str_extract(as.character(harvest_delay), "\\d+"))]
}

## RDS schedules use once_logged / twice_logged; routing expects once-logged / twice-logged
if (all(c("original_habitat", "habitat") %in% names(sched_prepped))) {
  sched_prepped[, original_habitat := nr2_normalize_habitat_labels(original_habitat)]
  sched_prepped[, habitat := nr2_normalize_habitat_labels(habitat)]
}
if ("functional_habitat" %in% names(sched_prepped)) {
  sched_prepped[, functional_habitat := nr2_normalize_habitat_labels(functional_habitat)]
}

## Order matches 04_temp.R prepare_scenarios(): update_1L, update_2L, filter_original_eq_habitat, remove_no_timber

if (!("functional_habitat" %in% names(sched_prepped)) && "habitat" %in% names(sched_prepped)) {
  sched_prepped[, functional_habitat := habitat]
}

## update_1L_functional_habitat
## ---------------------------------------------------------------------------
## Plain-English intent for this block:
## 1) We only target rows where a parcel *starts* as once-logged.
## 2) For once-logged -> restored rows:
##    - Re-label `functional_habitat` to "restored_start" so carbon lookup uses
##      the custom restored-start trajectory (not the generic restored/1L curve).
##    - Set `functionalhabAge = true_year` so lookup age tracks scenario year.
##    - Keep ONLY rows with harvest_delay == 0.
## 3) For once-logged -> once-logged rows:
##    - Keep ONLY rows with harvest_delay == 0.
## 4) All other transition types are untouched by this block.
##
## Why this is split into explicit booleans (`bad_rest`, `bad_1l`):
## - It makes it easy to inspect exactly which rows are being dropped.
## - You can `View(sched_prepped[bad_rest | bad_1l])` during debugging.
## ---------------------------------------------------------------------------
if (all(c("original_habitat", "habitat", "functional_habitat") %in% names(sched_prepped))) {
  # Ensure text comparisons/assignments are stable (factor columns can be awkward).
  sched_prepped[, functional_habitat := as.character(functional_habitat)]

  # Rule A: route once-logged -> restored onto restored_start curves.
  sched_prepped[original_habitat == "once-logged" & habitat == "restored", functional_habitat := "restored_start"]

  # Rule A (continued): when routed to restored_start, tie curve age to true_year.
  # This keeps the schedule and draw lookup on the same time axis.
  if ("true_year" %in% names(sched_prepped) && "functionalhabAge" %in% names(sched_prepped)) {
    sched_prepped[original_habitat == "once-logged" & habitat == "restored", functionalhabAge := as.numeric(true_year)]
  }

  # Convert once to numeric so delay checks behave consistently for character/factor input.
  hd <- as.numeric(sched_prepped$harvest_delay)

  # "Bad" once-logged -> restored rows are any with missing delay or delay != 0.
  bad_rest <- sched_prepped$original_habitat == "once-logged" & sched_prepped$habitat == "restored" &
    (is.na(hd) | hd != 0)

  # "Bad" once-logged -> once-logged rows are any with missing delay or delay != 0.
  bad_1l <- sched_prepped$original_habitat == "once-logged" & sched_prepped$habitat == "once-logged" &
    (is.na(hd) | hd != 0)

  # Drop exactly those rows that violate the two once-logged delay rules above.
  sched_prepped <- sched_prepped[!(bad_rest | bad_1l)]
}

## update_2L_functional_habitat
if (all(c("original_habitat", "functional_habitat", "harvest_delay") %in% names(sched_prepped))) {
  sched_prepped[original_habitat == "twice-logged" & functional_habitat == "twice-logged", functional_habitat := "twice_logged_start"]
  sched_prepped <- sched_prepped[!(original_habitat == "twice-logged" & !is.na(harvest_delay) & harvest_delay < 15)]
}

## filter_original_eq_habitat: if original_habitat == habitat, keep only harvest_delay == 0
if (all(c("original_habitat", "habitat", "harvest_delay") %in% names(sched_prepped))) {
  sched_prepped <- sched_prepped[!(original_habitat == habitat & (is.na(harvest_delay) | harvest_delay != 0))]
}

## remove_no_timber_scenarios
if ("production_target" %in% names(sched_prepped)) {
  sched_prepped <- sched_prepped[production_target > 0]
}

stopifnot(nrow(sched_prepped) > 0)

## -----------------------------------------------------------------------------
## 3) LOAD ONE DRAW of carbon curves (ACD lookup)
## -----------------------------------------------------------------------------

## nr2_out_root from _config.R (default: Desktop/carbon_data_plantation_models)
draws_path <- file.path(nr2_out_root, "02_draws", "rds", "acdraws_aboveground.rds")
stopifnot(file.exists(draws_path))
acd_draws <- readRDS(draws_path)
draw_dt <- as.data.table(copy(acd_draws))

## Match naming used elsewhere
if ("habitat" %in% names(draw_dt)) {
  draw_dt[habitat == "once_logged", habitat := "once-logged"]
}
if ("slope_factor" %in% names(draw_dt)) {
  draw_dt <- draw_dt[as.character(slope_factor) == as.character(TWICE_LOGGED_SLOPE_TRAJECTORY)]
}
if ("start_age" %in% names(draw_dt)) {
  draw_dt[, start_age := NULL]
}
setnames(draw_dt, old = "habitat", new = "functional_habitat", skip_absent = TRUE)
draw_dt <- draw_dt[draw == DRAW_ID]
stopifnot(nrow(draw_dt) > 0)

## -----------------------------------------------------------------------------
## 4) JOIN carbon curves onto schedule (explicit replication of 04_temp logic)
## -----------------------------------------------------------------------------

sched_dt <- as.data.table(copy(sched_prepped))

## Ensure functional_habitat exists (some schedules already have it; if not, default to habitat)
if (!("functional_habitat" %in% names(sched_dt))) {
  sched_dt[, functional_habitat := habitat]
}

## Split primary vs non-primary joins
x_primary <- sched_dt[functional_habitat == "primary"]
x_other <- sched_dt[functional_habitat != "primary"]

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

sched_carbon <- rbindlist(list(joined_primary, joined_other), use.names = TRUE, fill = TRUE)
stopifnot(nrow(sched_carbon) > 0)

## Optional helper (may not exist in all sessions)
if (exists("scenario_row_order_visualisation_fun", mode = "function")) {
  sched_carbon <- scenario_row_order_visualisation_fun(sched_carbon)
}

## -----------------------------------------------------------------------------
## 5) SPLIT parcels across transition delays (parcels_per_delay)
## -----------------------------------------------------------------------------
sched_carbon <- data.table(sched_carbon)
if (!("harvest_delay" %in% names(sched_carbon))) {
  ## If harvest_delay doesn't exist, create it as NA (meaning: no explicit delays)
  sched_carbon[, harvest_delay := NA_real_]
}

transitions <- sched_carbon[, .(num_transition_delays = uniqueN(harvest_delay)),
  by = .(index, production_target, original_habitat, habitat)
]
sched_carbon[transitions, on = .(index, production_target, original_habitat, habitat),
  num_transition_delays := i.num_transition_delays
]
sched_carbon[, parcels_per_delay := num_parcels / pmax(1, num_transition_delays)]

## Parcel accounting invariant check (should equal parcel_total = 1000 each year)
parcel_total_expected <- 1000
parcels_by_year <- sched_carbon[, .(parcels_total = sum(parcels_per_delay, na.rm = TRUE)), by = .(true_year)]

## -----------------------------------------------------------------------------
## 6) SCENARIO carbon stock by year (habitat → scenario)
## -----------------------------------------------------------------------------

## ACD (Mg/ha) → Mg per 10 km^2 parcel (1000 ha)
sched_carbon[, ACD_per_parcel := ACD * 1000]
sched_carbon[, ACD_10km2_stag := ACD_per_parcel * parcels_per_delay]

habitat_year <- sched_carbon[, .(
  hab_ACD_year = sum(ACD_10km2_stag, na.rm = TRUE)
), by = .(index, production_target, original_habitat, habitat, functional_habitat, harvest_delay, true_year)][
  order(true_year, habitat, harvest_delay, functional_habitat)
]

scen_year <- habitat_year[, .(
  scen_ACD_year = sum(hab_ACD_year, na.rm = TRUE)
), by = .(index, production_target, true_year)][order(true_year)]

## Stock-years (sum across years)
stock_years <- scen_year[, .(cumaltive_stock_year = sum(scen_ACD_year, na.rm = TRUE)),
  by = .(index, production_target)
]

## -----------------------------------------------------------------------------
## 7) STARTING LANDSCAPE carbon stock by year (for this scenarioStart)
## -----------------------------------------------------------------------------

fixed_path <- file.path("Inputs", "FixedScenarioParmams.R")
stopifnot(file.exists(fixed_path))
env_fixed <- new.env(parent = globalenv())
sys.source(fixed_path, envir = env_fixed)
stopifnot(exists("all_start_landscape", envir = env_fixed, inherits = FALSE))
all_start_landscape <- as.data.table(get("all_start_landscape", envir = env_fixed, inherits = FALSE))

scen_start <- unique(sched_carbon$scenarioStart)
scen_start <- as.character(scen_start[[1]])

## Expand starting landscape through time (explicit; no uncount)
years_vec <- YEARS
sl0 <- copy(all_start_landscape)[scenarioStart == scen_start]
stopifnot(nrow(sl0) > 0)

sl_exp <- sl0[rep(seq_len(.N), each = length(years_vec))]
sl_exp[, functional_habitat := habitat]
sl_exp[habitat == "once-logged", functional_habitat := "once-logged_start"]
sl_exp[habitat == "twice-logged", functional_habitat := "twice-logged_start"]
sl_exp[, functionalhabAge := rep(years_vec, times = nrow(sl0))]

## Join starting landscape parcels to draw curve and compute starting-landscape carbon
primary_dt <- draw_dt[functional_habitat == "primary"]
non_primary_dt <- draw_dt[functional_habitat != "primary"]

primary_joined <- sl_exp[primary_dt, on = .(functional_habitat), nomatch = 0L]
non_primary_joined <- sl_exp[non_primary_dt, on = .(functional_habitat, functionalhabAge), nomatch = 0L]

sl_merged <- rbindlist(list(primary_joined, non_primary_joined), use.names = TRUE, fill = TRUE)
if ("i.functionalhabAge" %in% names(sl_merged)) sl_merged[, `i.functionalhabAge` := NULL]

sl_merged[, SL_ACD_year := ACD * 1000 * num_parcels]
sl_year <- sl_merged[, .(SL_ACD_year = sum(SL_ACD_year, na.rm = TRUE)),
  by = .(scenarioStart, true_year = functionalhabAge)
][order(true_year)]

## -----------------------------------------------------------------------------
## 8) JOIN scenario vs starting landscape, then changes + flux
## -----------------------------------------------------------------------------

joined <- merge(
  scen_year[, scenarioStart := scen_start][],
  sl_year,
  by = c("scenarioStart", "true_year"),
  all.x = TRUE,
  sort = FALSE
)

setorder(joined, index, production_target, true_year)
changes <- copy(joined)
changes[, `:=`(
  scen_ACD_change = scen_ACD_year - shift(scen_ACD_year),
  SL_ACD_change = SL_ACD_year - shift(SL_ACD_year)
), by = .(index, production_target)]

mw <- 44.01 / 12.01
flux <- copy(changes)
flux[, `:=`(
  scen_flux_ACD = scen_ACD_change * mw,
  SL_flux_ACD = SL_ACD_change * mw
)]

## -----------------------------------------------------------------------------
## 9) SCC monetisation (annual + total) for ONE discount rate
## -----------------------------------------------------------------------------

scc_path <- file.path(nr2_out_root, "03_scc", "tables", "scc_dr_2_4_6.csv")
stopifnot(file.exists(scc_path))
socialDR <- read.csv(scc_path, stringsAsFactors = FALSE)

scc_dt <- socialDR %>%
  filter(discount_rate == DISCOUNT_RATE) %>%
  rename(true_year = year)

annual_scc <- as.data.frame(flux) %>%
  left_join(scc_dt, by = "true_year") %>%
  mutate(annual_carbon_impact_ACD = (scen_flux_ACD - SL_flux_ACD) * scc_discounted) %>%
  mutate(scenarioName = SCENARIO_NAME, scenarioStart = scen_start) %>%
  as.data.table()

total_scc <- annual_scc[, .(
  TOTcarbon_ACD_impact = sum(annual_carbon_impact_ACD, na.rm = TRUE)
), by = .(index, production_target)]

## -----------------------------------------------------------------------------
## 10) VISUAL CHECKS (optional plots)
## -----------------------------------------------------------------------------

## Area accounting (should be ~1000; if not, you’ve found a core accounting issue)
p_parcels <- ggplot(as.data.frame(parcels_by_year), aes(x = true_year, y = parcels_total)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = parcel_total_expected, linetype = 2) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Parcel accounting by year (expected ~", parcel_total_expected, ")"),
    x = "true_year",
    y = "sum(parcels_per_delay)"
  )

## Scenario vs starting landscape carbon (stock) over time
p_stock <- ggplot(as.data.frame(joined), aes(x = true_year)) +
  geom_line(aes(y = scen_ACD_year / 1e9, colour = "Scenario"), linewidth = 0.7, na.rm = TRUE) +
  geom_line(aes(y = SL_ACD_year / 1e9, colour = "Starting landscape"), linewidth = 0.7, na.rm = TRUE) +
  theme_bw(base_size = 12) +
  labs(
    title = "Carbon stock through time (one draw)",
    x = "true_year",
    y = "Carbon stock (billion Mg C; ACD*area)",
    colour = NULL
  )

## Restored-only (once-logged → restored) ACD, coloured by functional curve
p_restored <- sched_carbon %>%
  as.data.frame() %>%
  filter(original_habitat == "once-logged", habitat == "restored") %>%
  ggplot(aes(
    x = true_year,
    y = ACD,
    colour = functional_habitat,
    group = interaction(harvest_delay, functional_habitat)
  )) +
  geom_line(alpha = 0.7, linewidth = 0.5, na.rm = TRUE) +
  theme_bw(base_size = 12) +
  labs(
    title = "Restored ACD (once-logged → restored), by functional curve + delay",
    x = "true_year",
    y = "ACD (Mg/ha)",
    colour = "functional_habitat"
  )

## Annual SCC impacts
p_annual_scc <- ggplot(as.data.frame(annual_scc), aes(x = true_year, y = annual_carbon_impact_ACD / 1e9)) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Annual SCC impact (", DISCOUNT_RATE, ")"),
    x = "true_year",
    y = "Annual impact (billion $; ACD only)"
  )

message("Loaded step-through objects into your environment:")
message("  sched_raw, sched_prepped, draw_dt, sched_carbon, parcels_by_year, habitat_year, scen_year, stock_years, sl_year, joined, changes, flux, annual_scc, total_scc")
message("Plot objects:")
message("  p_parcels, p_stock, p_restored, p_annual_scc")

