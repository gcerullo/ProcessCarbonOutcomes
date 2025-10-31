# =============================================================================
# SCRIPT OVERVIEW: Scenario-by-Scenario Carbon and ACD Analysis
# =============================================================================
#NEED TO ADD BELOWGROUND?NECROMASS IMPACTS 31.10.25

# This loop iterates over each scenario individually to manage memory usage efficiently.
# For each scenario `i`:
#   1. Extract the scenario schedule (`scenarios_test`) and combine it with 
#      posterior draws of carbon estimates (`all_draws`) using `add_carbon_fun`.
#   2. Calculate annual ACD at multiple scales:
#        - Parcel-level (10 km² units)
#        - Habitat-transition by year
#        - Full scenario-level by year using `scenario_ACD_fun`
#   3. Compute total carbon stock-years for the scenario (`carbon_stock_years_fun`).
#   4. Incorporate starting landscape ACD trajectories (`all_SL_carbon`) and join 
#      them to scenario trajectories for comparative analysis.
#   5. Calculate annual changes in ACD and convert to carbon fluxes (`ACD_change_function_dt` 
#      and `flux_conversion_function`).
#   6. Monetize carbon fluxes using the Social Cost of Carbon for multiple discount 
#      rates (2%, 4%, 6%) via `socialDR_fun`.
#   7. Summarize posterior draws across all 500 draws for median and credible 
#      intervals (`summarize_totcarbon_draws`).
#   8. Store the final scenario-level summary in `posterior_summary_all[[i]]` for 
#      later consolidation and plotting.
#
# At the end of the loop, each element of `posterior_summary_all` corresponds to a 
# single scenario with all carbon, ACD, and discounted social cost outputs.
# =============================================================================


#propagate the uncertainty and estimates of scenario ACD 
#NOTE - need to go back and correct belowground process.
library(dplyr)
library(purrr)
library(tidyverse)
library(data.table)

#___________________________
#hardcoded decision
#___________________________
#SELECT TWICE-LOGGED RECOVERY TRAJECTORY- HARD-CODED
#1 = same slope as 1L, 0.8 is 20% slower recovery, 1.2 = 20% faster recovery

twice_logged_slope_trajectory = 1

#___________________________
#Inputs 
#___________________________
#Inputs params including the makeup of starting landscapes
source('Inputs/FixedScenarioParmams.R')
all_start_landscape <- all_start_landscape %>%  
  mutate(functional_habitat = habitat) %>%  
  mutate(
    functional_habitat = case_when(
      #deals with our rule that scenarios starting once-logged were logged at t-15 and twice-logged were at t-1 (see rationale below)
      habitat == "once-logged"  ~ "once-logged_start",
      habitat == "twice-logged" ~ "twice-logged_start",
      TRUE                       ~ habitat  # keep all other values unchanged
    )
  )

#social discount rate caluclate for 2,4,6% 
#built in the CalculateSocialDiscountRates.R script 
socialDR <- read.csv("Outputs/SocialDiscountRates_2_4_6pc_DR.csv")

 
#----------------read in scenarios -------------------------------

#yield matched scenarios where 1/30th of plantation conversion happens annually - with no time delay
#temporarily read in to get composition
scenarios_rm <- readRDS("Inputs/MasterAllScenarios.rds")
scenario_composition <- rbindlist(scenarios_rm, use.names=TRUE)

#yield matched scenarios where 1/30th of plantation conversion happens annually - WITH TIME DELAY 
scenarios <- readRDS("Inputs/MasterAllScenarios_withHarvestDelays.rds")

#function for making sure we can visualise the order of scenarios correctly
scenario_row_order_visualisation_fun <- function(x){
  x %>%  
    group_by(index, production_target, habitat, original_habitat, harvest_delay) %>%  
    arrange(true_year, .by_group = TRUE) %>%
    ungroup()
}

#___________________________
#read in posterior draws
#___________________________

#Build final posterior draws across habs ####
p_l_r <- readRDS("Outputs/primary_restored_once_logged_ACD_draws.rds") 
#currently has different twice-logged trajectories - from either clearing 15 yr or 30 once-logged.
twice <- readRDS("Outputs/twice_logged_draws_diff_assumptions.rds")
plant <- readRDS("Outputs/plantation_carbon_draws.rds") %>%  rename(habitat = species, functionalhabAge = plantationAge)

#_______________________________________________________________________________
#transition rules 
#_______________________________________________________________________________
#primary - no rules 

#once logged
#a. P -> 1L  If parcel goes to from primary to once-logged, we start at 1L yr 0
#b. 1L -> 1L if parcel starts as once-logged and stays 1L, assume logging happened 15 yrs befor scnanario start
#c. 1L -> 2L if parcel starts as once logged before re-harvest, we assume 1L recovery (already explicit in scenarios) that resets to yr 0 once-logged values.

#for b. we use the function update_1L_functional_habitat to ensure that scenarios beginning the scenario as once-logged were harvest 15yrs before

#if parcel starts as twice-logged, use the 15yr after 1L curves - with different slope factors, tracking the different slope assumption rates for twice-logged forest 
#if parcel starts as primary and goes to twice-logged, use the 30yrafter 1L, with different slope factors 
#_________________________________________________________________
# modify scenarios to accommodate different harvest trajectories
#_____________________________________________________________

#if original habitat remains as habitat (ie no transition), select one delay year 
filter_original_eq_habitat <- function(df_list) {
  df_list %>%
    map(~ .x %>%
          filter(!(original_habitat == habitat & harvest_delay != "delay 0"))
    )
}


#___________________________________________________________________
#A) 1L- 1L modification. Forest that starts a scenario as once-logged was first harvested at yr t-15

# changes `functional_habitat` to "once_logged_start" when both
# `original_habitat` and `functional_habitat` are "once_logged".
# Otherwise, leaves `functional_habitat` unchanged.

# Function: update_functional_habitat()
#_______________________________________
# Takes a list of data frames and applies two operations to each:
# (1) If both `original_habitat` and `functional_habitat` are "once_logged",
#     change `functional_habitat` to "once_logged_start".
# (2) If `original_habitat` is "once_logged", remove rows where `harvest_delay > 15`.
update_1L_functional_habitat <- function(df_list) {
  df_list %>%
    map(~ .x %>%
          mutate(
            functional_habitat = if_else(
              original_habitat == "once-logged" & functional_habitat == "once-logged",
              "once-logged_start",
              functional_habitat
            )
          ) %>%
          filter(!(original_habitat == "once-logged" & harvest_delay < 15))
    )
}

#___________________________________________________________________
# B)2L - 2L modification. Forest that starts with twice-logged was re-harvest just before scenario start
# Takes a list of data frames and applies two operations to each:
# (1) If both `original_habitat` and `functional_habitat` are "twice-logged",
#     change `functional_habitat` to "twice_logged_start".
# (2) If `original_habitat` is "once_logged", remove rows where `harvest_delay < 15`.- because you have to wait 30 yrs before can reharvest once-logged forest 
update_2L_functional_habitat <- function(df_list) {
  df_list %>%
    map(~ .x %>%
          mutate(
            functional_habitat = if_else(
              original_habitat == "twice-logged" & functional_habitat == "twice-logged",
              "twice-logged_start",
              functional_habitat
            )
          ) %>%
          filter(!(original_habitat == "twice-logged" & harvest_delay < 15))
    )
}


#remove scenarios that produce zero timber 

remove_no_timber_scenarios <- function(df_list) {
  df_list %>%
    map(~ .x %>%
          filter(production_target > 0)  # keep only scenarios producing timber
    )
}

scenarios <- scenarios %>%  update_1L_functional_habitat() %>% 
  update_2L_functional_habitat() %>%  
  filter_original_eq_habitat() %>% 
  remove_no_timber_scenarios()

#convert to datatable format
scenarios <- lapply(scenarios, as.data.table)


# #check that we have the correct number of harvest delays per transition!!!
# m <- scenarios[[9]]
# p <- m %>% group_by(index, production_target, functionalhabAge, habitat, original_habitat) %>% count()
# unique(p$n)
# hist(p$n)
#q <- m %>% filter(index == "mostly_2L_deforested_CY_D.csv 1") %>% unique() 
#q2 <- m %>% filter(index == "mostly_2L_deforested_CY_D.csv 291") 
  
#________________________________________________________________
#make correct modifications to each draw
#__________________________________________________________________
primary <- p_l_r %>% filter((habitat == "primary"))
restored <- p_l_r %>% filter((habitat == "restored"))

once_logged <-  p_l_r %>% filter((habitat == "once_logged")) %>%  mutate(habitat = if_else(habitat == "once_logged", "once-logged", habitat))
#for scenarios starting once-logged assume the forest was logged 15yrs before scneario data
once_logged_start <- p_l_r %>% filter((habitat == "once_logged")) %>% 
  mutate(functionalhabAge = functionalhabAge -15) %>%  
  filter(functionalhabAge >-1) %>%  
  mutate(habitat = "once-logged_start")

#2L ->2L forest that starts as twice-logged was first logged fifteen years previously then reharvested just before scenario start
twice_logged_start <-twice %>% filter(start_age == "15yrAfter1L - e.g if parcel starts scenario 2L") %>%  
  mutate(functionalhabAge = functionalhabAge -15) %>%  
  filter(functionalhabAge >-1) %>% 
  mutate(habitat = "twice-logged_start")

#1L-2L forest that goes from once-logged to twice-logged during scenarios
#not allowed to be harvest for 15yrs)

#second harvest is not allowed for first 15 yrs, then it goes to twice-logged(30 yr after once-logging version)
twice_logged <-twice %>% filter(start_age == "30yrAfter1L - e.g if primary goes to 2L during scenario") %>% 
  mutate(habitat = "twice-logged")

#____________________________________________
#combine all draws
#____________________________________________
#Note primary forest has NA for functionalhabAge
all_draws <- primary %>% rbind(restored) %>% rbind(once_logged) %>% rbind(once_logged_start) %>%  
   rbind(plant)

#SELECT TWICE-LOGGED RECOVERY TRAJECTORY- HARD-CODE
#!!!!!select which twice-logged recovery trajectory to consider using 'slope factor' - 
#1 = same slope as 1L, 0.8 is 20% slower recovery, 1.2 = 20% faster recovery
diff_2Ls <- twice_logged %>% rbind(twice_logged_start) %>% 
  filter(slope_factor == twice_logged_slope_trajectory) %>% 
  select( draw, functionalhabAge, ACD, habitat) %>% unique()

all_draws<- all_draws %>% rbind(diff_2Ls)

#turn draws into a list where each draw is it's own df
all_draws <- all_draws %>%
  rename(functional_habitat = habitat) %>%  # rename first
  group_split(draw)     

#convert to datatable
all_draws <- lapply(all_draws, as.data.table)


#____________________________________________________________________________________________________________________________________
#We are now ready to run the full workflow below
# For each scenario set i, we extract the scenario schedule and sequentially apply the full workflow:
#   1) Add carbon estimates from all posterior draws to the scenario landscape.
#   2) Calculate ACD (aboveground carbon density) at the parcel, habitat-transition, and scenario-year scales.
#   3) Compute carbon stock-years for the scenario landscape.
#   4) Add starting landscape carbon trajectories and calculate annual ACD changes for both scenario and starting landscapes.
#   5) Convert annual ACD changes into CO2-equivalent fluxes and apply Social Cost of Carbon (SCC) discounting for multiple discount rates (2%, 4%, 6%).
#   6) Summarize posterior draws to calculate mean and 95% / 80% credible intervals of carbon fluxes.
# At the end of each iteration, the processed results for scenario i are stored in the list `posterior_summary_all[[i]]`,
# containing SCC-discounted carbon impacts with uncertainty and cumulative carbon stock-years.
#___________________________________________________________________________________________________________________________________

#make empty list to store results
posterior_summary_all <- vector("list", length(scenarios))


#to test running for 1 scenario uncomment the following: 
# i = 1

for(i in seq_along(scenarios)) {
  
  cat("Running scenario", i, "of", length(scenarios), "\n")
  
  # pick the i-th scenario
  scenario_i <- scenarios[[i]]

   #pick just one scenario-set for  now 
    #scenario_i <- scenarios[[1]]
    #but use all draws of the carbon estimates
     draw <- all_draws

#__________________________________________
#combine scenarios with carbon values
#__________________________________________

#function to add carbon estimates to scenario schedule
add_carbon_fun <- function(x, draw_dt) {
  # Ensure both are data.tables
  setDT(x)
  setDT(draw_dt)
  
  # Separate primary and non-primary to handle joins cleanly
  x_primary <- x[functional_habitat == "primary"]
  x_other   <- x[functional_habitat != "primary"]
  
  # Join scenario primary carbon only - only by functional habitat as functionalhabAge is NA for primary 
  joined_primary <- x_primary[
    draw_dt[functional_habitat == "primary"],
    on = .(functional_habitat),
    nomatch = 0
  ]
  
  # Remove additional column created by join
  joined_primary[, i.functionalhabAge := NULL]
  
  # Join non-primary habitats to carbon values 
  joined_other <- x_other[
    draw_dt[functional_habitat != "primary"],
    on = .(functional_habitat, functionalhabAge),
    nomatch = 0
  ]
  
  # Combine primary and non-primary to get carbon results of each row
  result <- rbindlist(list(joined_primary, joined_other), use.names = TRUE)
  
  return(result)
}

# Now apply the add carbon function for each posterior draw of the datatable list:
# Apply add_carbon_fun() to each draw
scenario_i <- lapply(draw, \(d) add_carbon_fun(scenario_i, d))
#check <- scenario_i[[1]] %>% scenario_row_order_visualisation_fun()

#________________________________________________________________
# Calculate ACD per scenario year 
#________________________________________________________________
#next - apply the next sttep in the carbon calculation for the scenario

# --- Function to calculate number of unique delays per transition ---
calculate_delays_per_transition_fun <- function(x) {
  setDT(x)  # ensure it's a data.table
  
  # Count unique harvest_delay values per transition
  transitions <- x[, .(num_transition_delays = uniqueN(harvest_delay)), 
                   by = .(index, production_target,original_habitat, habitat)]
  
  # Join back to original table
  x[transitions, on = .(index, production_target,original_habitat, habitat), num_transition_delays := i.num_transition_delays]
  
  return(x)
}

# --- Apply function to each list element ---
scenario_i <- lapply(scenario_i, calculate_delays_per_transition_fun)
#check <- scenario_i[[1]] %>% scenario_row_order_visualisation_fun()

#_______________________________________________
# Function: scenario_ACD_fun
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

#scenario_test_lowP <- check %>% filter(index == "all_primary_CY_D.csv 10")
#scenario_test_highP <- check %>% filter(index == "unq2_1055")
#x <- scenario_test_lowP

scenario_ACD_fun <- function(x) {
  
  # --- Step 0: Ensure input is a data.table ---
  setDT(x)
  
  # --- Step 1: Calculate ACD values at 10 km² parcel scale --------------------
  # ACD_10km:  Aboveground carbon density (Mg/ha) scaled to Mg per 10 km² area
  # ACD_10km2_stag:  ACD_10km adjusted by number of parcels and divided by
  #                  the number of unique transition delays for this transition.
  #                  This distributes carbon proportionally across delays for areas undergoing harvest.
  #                  Note that parcels that remain unharvested don't have a delsy- so this area remains the same throughout ( ie. is only /1)
  x[, ACD_10km := ACD * 1000]  # Convert ACD (Mg/ha) to Mg per 10 km²
  x[, ACD_10km2_stag := ACD_10km * (num_parcels / num_transition_delays)]
  
  
  # --- Step 2: Aggregate to habitat-transition level per year -----------------
  # For each combination of scenario index, production target, habitat type,
  # and year, sum up all the parcel-level ACD to get total habitat-level ACD.
  # Keep metadata for reference (scenario name and start year).
  habitat_year_summary <- x[, .(
    
    # Sum all parcel-level ACD within the transition-year group
    hab_ACD_year = sum(ACD_10km2_stag, na.rm = TRUE),
    
    # Retain scenario metadata (assuming constant within group)
    scenarioName = first(scenarioName),
    scenarioStart = first(scenarioStart)
    
  ), by = .(index, production_target, original_habitat, habitat, true_year)]
  
  
  # --- Step 3: Aggregate to full-scenario level per year ----------------------
  # Sum across all habitat transitions to get the total ACD in each scenario
  # and production target combination for a given year.
  scenario_year_summary <- habitat_year_summary[, .(
    
    # Sum all habitat-level ACD values for this year
    scen_ACD_year = sum(hab_ACD_year, na.rm = TRUE),
    
    # Retain metadata
    scenarioName = first(scenarioName),
    scenarioStart = first(scenarioStart)
    
  ), by = .(index, production_target, true_year)]
  
  
  # --- Step 4: Return scenario-level summary ---------------------------------
  # This data.table contains total ACD per year for each scenario and production target.
  return(scenario_year_summary)
}

# Apply scenario_ACD_fun to all scenarios
scenario_i <- lapply(scenario_i, scenario_ACD_fun)
#x <- scenario_i[[1]]

#CALCULATE CARBON STOCK YEARS ####
#____________________________________________________________
# Function: carbon_stock_years_fun
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
#   production target:
#        - production_target
#        - index
#        - stock_year : sum of annual ACD values across all years in scenario
#___________________________________________________________________________

carbon_stock_years_fun <- function(x) {
  
  # --- Step 0: Ensure data.table format --------------------------------------
  setDT(x)
  
  # --- Step 1: Summarize across years for each scenario ----------------------
  # Group by scenario ID and production target to sum ACD across all years.
  # Retain scenario metadata (assumed constant within each group).
  carbon_summary <- x[, .(
    
    # Sum of annual aboveground carbon across all years (Mg)
    cumaltive_stock_year = sum(scen_ACD_year, na.rm = TRUE),
    
    # Retain descriptive metadata for easier tracking
    scenarioName = first(scenarioName),
    scenarioStart = first(scenarioStart)
    
  ), by = .(production_target, index)]
  
  # --- Step 2: Return summarized data ---------------------------------------
  return(carbon_summary)
}

carbon_yrs <- lapply(scenario_i, carbon_stock_years_fun)
#rbindlist
carbon_yrs <- rbindlist(carbon_yrs, use.names = TRUE, fill = TRUE, idcol = "draw")

#get mean scenario stock years with uncertainty:

# --- Summarize across draws for each scenario (index × production_target) ---
carbon_stock_years_summary <- carbon_yrs[, .(
  
  # Mean of stock-years across draws
  mean_cum_stock_year = mean(cumaltive_stock_year, na.rm = TRUE),
  
  # 95% credible interval
  lwr_cum_stock_year_95 = quantile(cumaltive_stock_year, probs = 0.025, na.rm = TRUE),
  upr_cum_stock_year_95 = quantile(cumaltive_stock_year, probs = 0.975, na.rm = TRUE),
  
  # 80% credible interval
  lwr_cum_stock_year_80 = quantile(cumaltive_stock_year, probs = 0.10, na.rm = TRUE),
  upr_cum_stock_year_80 = quantile(cumaltive_stock_year, probs = 0.90, na.rm = TRUE),
  
  # Optional metadata (if constant across draws)
  scenarioName = first(scenarioName),
  scenarioStart = first(scenarioStart)
  
), by = .(index, production_target)]

#__________________________________________
# Get starting landscape carbon trajectories 
#__________________________________________
#print(all_start_landscape)
#print(all_draws)

all_start_landscape<- all_start_landscape %>%
  # Repeat each row 61 times (0:60)
  uncount(weights = 61) %>%
  # Add functionalhabAge column from 0 to 60
  group_by(across(everything())) %>%   # preserve other columns
  mutate(functionalhabAge = 0:60) %>%
  ungroup()


#add carbon trjaectories for each starting landscape 
# Ensure all_start_landscape_expanded is a data.table
setDT(all_start_landscape)

all_SL_carbon <- lapply(all_draws, function(draw_dt) {
  
  setDT(draw_dt)
  setDT(all_start_landscape)
  
  # Split primary vs non-primary
  primary_dt     <- draw_dt[functional_habitat == "primary"]
  non_primary_dt <- draw_dt[functional_habitat != "primary"]
  
  # Join primary: only by functional_habitat (because functionalAge is NA but we still need year-level estimates)
  primary_joined <- all_start_landscape[primary_dt, 
                                        on = .(functional_habitat), 
                                        nomatch = 0L]
  
  # Join non-primary: by functional_habitat AND functionalhabAge
  non_primary_joined <- all_start_landscape[non_primary_dt, 
                                            on = .(functional_habitat, functionalhabAge), 
                                            nomatch = 0L]
  
  # Combine both
  merged_dt <- rbindlist(list(primary_joined, non_primary_joined), use.names = TRUE, fill = TRUE)
  
  # Remove unwanted i.* column
  if ("i.functionalhabAge" %in% names(merged_dt)) merged_dt[, `i.functionalhabAge` := NULL]
  
  return(merged_dt)
})

#x <- all_SL_carbon[[1]]
SL_ACD_by_age <- map(all_SL_carbon, ~ .x %>%
                       #  Calculate scaled annual carbon per parcel
                       mutate(
                         SL_ACD_year = ACD * 1000 * num_parcels) %>%
                       # Group by functionalhabAge, draw, and scenarioStart
                       group_by(functionalhabAge, draw, scenarioStart) %>%
                       #  Sum ACD across parcels within group
                       summarise(
                         SL_ACD_year = sum(SL_ACD_year, na.rm = TRUE),
                         .groups = "drop"
                       ) %>% rename(true_year = functionalhabAge)
)
#x <- SL_ACD_by_age[[1]]

# calculate SL cumulative stock  years 
SL_stockyears_list <- map(all_SL_carbon, ~ .x %>%
                            # Step 1: Calculate annual carbon per landscape
                            mutate(
                              yearACD = ACD * 1000 * num_parcels#
                            ) %>%
                            #  Group by scenarioStart and sum stock-years
                            group_by(scenarioStart) %>%
                            summarise(
                              cumulative_stock_yearSL = sum(yearACD, na.rm = TRUE),
                              .groups = "drop"
                            )
)
SL_stockyears <- bind_rows(SL_stockyears_list, .id = "draw") %>%
  mutate(draw = as.integer(draw))
setDT(SL_stockyears)

# Summarise across draws
SL_carbon_stock_years_summary <- SL_stockyears[, .(
  
  # Mean of stock-years across draws
  mean_cum_stock_year = mean(cumulative_stock_yearSL, na.rm = TRUE),
  
  # 95% credible interval (2.5th and 97.5th percentiles)
  lwr_cum_stock_year  = quantile(cumulative_stock_yearSL, probs = 0.025, na.rm = TRUE),
  upr_cum_stock_year  = quantile(cumulative_stock_yearSL, probs = 0.975, na.rm = TRUE),
  
  # Optional metadata (if constant across draws)
  scenarioStart = first(scenarioStart)
  
), by = .(scenarioStart )]

#__________________________________________________
#Calculate fluxes within starting and scenario landscape ####
#__________________________________________________
#print(SL_ACD_by_age[[1]]) #starting landscape annual ACD
#print(scenario_i[[1]]) #scenario landscape annual ACD 
#unique(SL_ACD_by_age[[1]]$scenarioStart)
#x <- SL_ACD_by_age[[1]]
#y <- scenario_i[[1]]

#____________________________________________________
# Function: join_scenario_with_start_ACD
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

join_scenario_with_start_ACD <- function(scenario_list, start_ACD_list,
                                         join_cols = c("scenarioStart", "true_year")) {
  
  # Ensure input lists have the same length
  if (length(scenario_list) != length(start_ACD_list)) {
    stop("Both lists must have the same length")
  }
  
  # Convert all elements to data.table if not already
  scenario_list <- lapply(scenario_list, setDT)
  start_ACD_list <- lapply(start_ACD_list, setDT)
  
  # Join each draw with corresponding starting landscape
  joined_list <- Map(function(scenario_dt, start_dt) {
    
    # Remove draw column if present to avoid confusion
    scenario_dt[, draw := NULL]
    start_dt[, draw := NULL]
    
    # Perform left join by specified columns
    merged_dt <- merge(
      scenario_dt,
      start_dt,
      by = join_cols,
      all.x = TRUE,   # keep all scenario rows
      sort = FALSE
    )
    
    # Remove duplicate rows if any
    unique(merged_dt)
    
  }, scenario_list, start_ACD_list)
  
  return(joined_list)
}


# Apply the function to  500-draw lists
scenario_i <- join_scenario_with_start_ACD(
  scenario_list = scenario_i,
  start_ACD_list = SL_ACD_by_age
)

#x <- scenario_i

#________________________________________________
#CHANGES IN ACD #####
#________________________________________________

# --------------------------------------------------------------------------------
# Function: ACD_change_function_dt
# Purpose:  Calculate annual changes in Aboveground Carbon Density (ACD) and total 
#           carbon for scenario landscapes and the starting landscape.
# Reason:   Clean, efficient version for posterior draws where error columns are 
#           no longer needed. Works on a data.table and can be applied to lists of
#           scenarios/draws.
# --------------------------------------------------------------------------------
ACD_change_function_dt <- function(dt) {
  
  # Ensure the input is a data.table (fast in-place operations)
  setDT(dt)
  
  # Step 1: Sort data by scenario index, production target, and year
  # This ensures that the lag calculations are correct within each group
  setorder(dt, index, production_target, true_year)
  
  # Step 2: Calculate annual changes by group
  # Using `shift()` to get the previous year's value
  # Grouped by `index` and `production_target` to calculate changes within each scenario
  dt[, `:=`(
    
    # Aboveground carbon change in the scenario landscape
    scen_ACD_change = scen_ACD_year - shift(scen_ACD_year),
    
    # Aboveground carbon change in the starting landscape
    SL_ACD_change = SL_ACD_year - shift(SL_ACD_year)
    
  ), by = .(index, production_target)]
  
  # Step 3: Keep only the relevant columns for downstream analysis
  dt[, .(
    index,              # unique scenario identifier
    production_target,  # target production intensity for this scenario
    true_year,          # year of the scenario
    scen_ACD_change,    # annual change in aboveground carbon (scenario)
    SL_ACD_change,      # annual change in aboveground carbon (starting landscape)
    scenarioName,       # descriptive scenario name
    scenarioStart       # starting landscape name
  )]
}

#Apply the function to a list of scenario draws
scenario_i<- lapply(scenario_i, ACD_change_function_dt)


#_________________________________________________
#visualise ACD changes for a subset of scenarios 
#_________________________________________________
# scenarios_df <- rbindlist(scenario_i, idcol = "draw")
# 
# 
# # --- Summarise posterior draws across draws for each scenario, year, and index ---
# scenarios_summary <- scenarios_df %>%
#   filter(production_target == 0.5) %>%
#   filter(scenarioName %in% c("all_primary_CY_D.csv", 
#                              "mostly_1L_CY_D.csv", 
#                              "mostly_2L_CY_D.csv")) %>%
#   group_by(scenarioName, true_year, index) %>%
#   summarise(
#     # posterior summaries (median and 95% credible interval), with NA handling
#     scen_med = median(scen_ACD_change, na.rm = TRUE) / 1e6,
#     scen_lwr = quantile(scen_ACD_change, 0.025, na.rm = TRUE) / 1e6,
#     scen_upr = quantile(scen_ACD_change, 0.975, na.rm = TRUE) / 1e6,
#     SL_med   = median(SL_ACD_change, na.rm = TRUE) / 1e6,
#     SL_lwr   = quantile(SL_ACD_change, 0.025, na.rm = TRUE) / 1e6,
#     SL_upr   = quantile(SL_ACD_change, 0.975, na.rm = TRUE) / 1e6,
#     .groups = "drop"
#   )
# 
# # --- Plot posterior median with 95% credible interval ribbons ---
# ggplot(scenarios_summary, aes(x = true_year, y = scen_med, colour = as.factor(index))) +
#   
#   # 95% credible interval ribbon for scenario trajectories
#   geom_ribbon(
#     aes(ymin = scen_lwr, ymax = scen_upr, fill = as.factor(index)),
#     alpha = 0.25, colour = NA
#   ) +
#   
#   # posterior median line for scenarios
#   geom_line(linewidth = 1) +
#   
#   # dashed black posterior median for SL reference trajectory
#   geom_line(aes(y = SL_med), linetype = "longdash", colour = "black", linewidth = 1) +
#   
#   # ribbon for SL reference 95% credible interval — explicitly map x
#   geom_ribbon(
#     aes(x = true_year, ymin = SL_lwr, ymax = SL_upr),
#     fill = "black", alpha = 0.15, inherit.aes = FALSE
#   ) +
#   
#   facet_wrap(~scenarioName, ncol = 3, scales = "free_y") +
#   theme_bw(base_size = 16) +
#   theme(
#     strip.background = element_blank(),
#     strip.text = element_text(hjust = 0, face = "bold"),
#     legend.position = "none",
#     panel.grid = element_blank(),
#     axis.text = element_text(colour = "black"),
#     axis.text.x = element_text(angle = 45, hjust = 1.05, colour = "black")
#   ) +
#   labs(
#     y = "Posterior ACD change (million Mg C)",
#     x = "",
#     title = "Posterior median and 95% credible intervals of annual ACD change"
#   )

#________________________________________________________
#convert changes in ACD into carbon fluxes
#________________________________________________________
#________________________________________________________
# Function: flux_conversion_function
# Purpose: Convert annual ACD (Aboveground Carbon Density) changes 
#          into corresponding carbon fluxes (CO2-equivalent)
#________________________________________________________

# Molecular weight ratio to convert C → CO2 (gCO2 per gC)
mw <- 44.01 / 12.01

flux_conversion_function <- function(x) {
  #' Convert ACD change to annual carbon fluxes
  #'
  #' @param x data.frame or data.table containing at least:
  #'   - `scen_ACD_change`: annual change in ACD for the scenario (Mg C)
  #'   - `SL_ACD_change`: annual change in ACD for the starting landscape (Mg C)
  #'
  #' @return A modified data.frame/data.table with two new columns:
  #'   - `scen_flux_ACD`: annual CO2-equivalent flux (Mg CO2)
  #'   - `SL_flux_ACD`: annual CO2-equivalent flux for starting landscape (Mg CO2)
  #'
  #' @details
  #' The conversion is based on the molecular weight ratio of CO2 to C:
  #' 44.01 (CO2) / 12.01 (C) = 3.664.
  #' This assumes complete oxidation of carbon to CO2.
  #'
  #' @examples
  #' df <- data.frame(
  #'   scen_ACD_change = c(10, -5),
  #'   SL_ACD_change = c(8, -3)
  #' )
  #' flux_conversion_function(df)
  #' 
  #' # Returns CO2-equivalent fluxes in Mg CO2
  #' # scen_flux_ACD = scen_ACD_change * 3.664
  #' # SL_flux_ACD   = SL_ACD_change * 3.664
  
  x %>%
    mutate(
      # Annual CO2-equivalent fluxes (Mg CO2 yr⁻¹)
      scen_flux_ACD = scen_ACD_change * mw,
      SL_flux_ACD   = SL_ACD_change   * mw
    )
}

scenario_i<- lapply(scenario_i, flux_conversion_function)
#x <- scenario_i[[1]]

#______________________________________________________
#Determine monetary social cost of carbon of  fluxes
#______________________________________________________


#  Prepare Social Cost of Carbon (SCC) tables
# SCC datasets contain discounted values of the social cost of carbon (scc_discounted)
# for different discount rates (2%, 4%, 6%) and years.
# We'll rename 'year' -> 'true_year' for clean joins later.

socialDR_2 <- socialDR %>%
  filter(discount_rate == "2%") %>%
  rename(true_year = year)

socialDR_4 <- socialDR %>%
  filter(discount_rate == "4%") %>%
  rename(true_year = year)

socialDR_6 <- socialDR %>%
  filter(discount_rate == "6%") %>%
  rename(true_year = year)

#_______________________________________
# Define core function: socialDR_fun()
#_______________________________________
# This function merges scenario fluxes with the chosen SCC table (SDR),
#  then applies the SCC discounting to compute annual and total monetized
# carbon impacts (for both total and aboveground carbon).

socialDR_fun <- function(x, SDR) {
  x %>%
    left_join(SDR, by = "true_year") %>%
    
    mutate(annual_carbon_impact_ACD = (scen_flux_ACD - SL_flux_ACD) * scc_discounted) %>%
    
    # Summarise across all years per scenario
    group_by(index, production_target, scenarioStart, scenarioName, true_year) %>%
    summarise(
      TOTcarbon_ACD_impact = sum(annual_carbon_impact_ACD, na.rm = TRUE),
      .groups = "drop"
    )
}



# Wrap socialDR_fun() for each discount rate
final_carbon_fun <- function(x, SDR) {
  socialDR_fun(x, SDR) 
}

# Apply to all scenarios and discount rates
final_carbon2 <- lapply(scenario_i, final_carbon_fun, SDR = socialDR_2) %>% rbindlist(idcol = "draw")
final_carbon4 <- lapply(scenario_i, final_carbon_fun, SDR = socialDR_4) %>% rbindlist(idcol = "draw")
final_carbon6 <- lapply(scenario_i, final_carbon_fun, SDR = socialDR_6) %>% rbindlist(idcol = "draw")

#now calculate summary (mean and 95% CIs across the 500 posterior draws)
#names(final_carbon2)

#----------------------------------------------
# Function: summarize_totcarbon_draws
# Purpose: Summarize TOTcarbon_ACD_impact across posterior draws
#----------------------------------------------
summarize_totcarbon_draws <- function(df) {
  df %>%
    group_by(index, production_target, scenarioStart, scenarioName) %>%
    summarise(
      # Posterior mean of total carbon impact
      TOTcarbon_ACD_mean = mean(TOTcarbon_ACD_impact, na.rm = TRUE),
      
      # 95% credible interval
      TOTcarbon_ACD_lwr95 = quantile(TOTcarbon_ACD_impact, 0.025, na.rm = TRUE),
      TOTcarbon_ACD_upr95 = quantile(TOTcarbon_ACD_impact, 0.975, na.rm = TRUE),
      
      # 80% credible interval
      TOTcarbon_ACD_lwr80 = quantile(TOTcarbon_ACD_impact, 0.10, na.rm = TRUE),
      TOTcarbon_ACD_upr80 = quantile(TOTcarbon_ACD_impact, 0.90, na.rm = TRUE),
      
      .groups = "drop"
    )
}


# Summarise across posterior draws
posterior_summary2 <- summarize_totcarbon_draws(final_carbon2)
posterior_summary4 <- summarize_totcarbon_draws(final_carbon4)
posterior_summary6 <- summarize_totcarbon_draws(final_carbon6)


# combine all discount rates into a single df
posterior_summary_comb <- bind_rows(
  posterior_summary2 %>% mutate(discount_rate = "2%"),
  posterior_summary4 %>% mutate(discount_rate = "4%"),
  posterior_summary6 %>% mutate(discount_rate = "6%")
)

#add back in cumulative carbon stock years alongside SCC-discounted carbon impacts 
posterior_summary_comb <- posterior_summary_comb %>%  left_join(carbon_stock_years_summary)

# at the very end, store the output
posterior_summary_all[[i]] <- posterior_summary_comb

# optionally remove large objects to free memory before next iteration
#rm(scenario_i, final_carbon2, final_carbon4, final_carbon6)
#gc()

}

posterior_summary_all
#___________________________________________
#Save outputs for consolidated final figure
#___________________________________________


#-----EXPORT OUTCOME PERFORMANCE for consolidated figure of all outcomes -----
getwd()
names(posterior_summary_comb)

#save final outputs....




