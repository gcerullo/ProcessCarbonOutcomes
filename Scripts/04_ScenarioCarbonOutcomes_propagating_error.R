
# 08.04.25; 
#NB; this code:
#this code calculates carbon consequences of different scenarios
#commented out code-blocks enable to assessment of a single scenario, to see how the code works.
#This code version propagates error properly 

library(tidyr)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)
library(stringr) 
library(cowplot)

#Reading in carbon by hab outcomes from  CalculateAllHabCarbonVals.R script 

#Read in Inputs ####
#read in the scenario parametres containing conversion factors for converting from point to parcel/entire landscape  
source('Inputs/FixedScenarioParmams.R')

# social discount rate caluclate for 2,4,6% 
#built in the CalculateSocialDiscountRates.R script 
socialDR <- read.csv("Outputs/SocialDiscountRates_2_4_6pc_DR.csv")

#define all starting landscapes: 
all_start_landscape 

#Import inputs

#Read in Data ####
#----------------read in scenarios -------------------------------

#yield matched scenarios where 1/30th of plantation conversion happens annually - with no time delay
#temporarily read in to get composition
scenarios_rm <- readRDS("Inputs/MasterAllScenarios.rds")
scenario_composition <- rbindlist(scenarios_rm, use.names=TRUE)

#yield matched scenarios where 1/30th of plantation conversion happens annually - WITH TIME DELAY 
scenarios <- readRDS("Inputs/MasterAllScenarios_withHarvestDelays.rds")

#------read in carbon by year per hab -----------------------
#read in carbon, with delays already calculated. ACD refers to just above-ground carbon 
#all_carbon incorporates belowground/necromass processes. 
#This is the output of the CalculateAllHabCarbonVals.R script
hab_carbon <- read.csv("Outputs/allHabCarbon_60yrACD.csv") %>% select(-X)

# Define function for twice-logged correction
# This fixes the fact that there is always a missing true year on the year 1L becomes 2L as Primary -> 2L transitions
adjust_twice_logged <- function(df) {
  df %>%
    filter(habitat == "twice-logged") %>%
    # Identify rows matching the condition
    filter(original_habitat == "primary",
           habitat == "twice-logged",
           functional_habitat == "twice-logged",
           functionalhabAge == 1) %>%
    # Duplicate these rows and modify the necessary columns
    mutate(
      functionalhabAge = functionalhabAge - 1,  # Adjust the age
      true_year = true_year - 1  # Adjust the year
    ) %>%
    # Bind the duplicated and modified rows back to the original dataset
    bind_rows(df) %>%  
    group_by(index, production_target, original_habitat, habitat) %>% 
    arrange(true_year, harvest_delay) 
}

# Convert hab_carbon to data.table for better performance
hab_carbon <- as.data.table(hab_carbon)

# Compute standard errors from upper/lower bounds assuming symmetric errors
hab_carbon <- hab_carbon %>%
  mutate(
    full_carbon_err = (full_carbon_upr - full_carbon_lwr) / 2,
    ACD_err = (upr_ACD - lwr_ACD) / 2
  )

# Correct scenarios ####
adjust_scenarios <- function(data) {
  data %>%
    mutate(
      functionalhabAge = if_else(functional_habitat %in% c("primary", "deforested"), 0, functionalhabAge),
      harvest_delay_numeric = as.numeric(stringr::str_extract(harvest_delay, "\\d+"))
    ) %>%
    filter(!(original_habitat == "once-logged" & habitat == "twice-logged" & harvest_delay_numeric < 15)) %>%
    select(-harvest_delay_numeric)
}

scenarios <- lapply(scenarios, adjust_scenarios)

# Join carbon information to scenarios ####
carbon_fun <- function(x) {
  x <- as.data.table(x)
  # Perform the left join operation in data table 
  result <- x[hab_carbon, on = .(habitat = habitat,
                                 original_habitat = original_habitat, 
                                 functional_habitat = functional_habitat, 
                                 functionalhabAge = functionalhabAge), nomatch = 0]
  
  return(result)
}

scenarios <- lapply(scenarios, carbon_fun)

# Apply our correction to ensure primary - 2L transitions have the correct number of years
scenarios <- lapply(scenarios, adjust_twice_logged)

# Check we have the same number of years (not allowed to harvest in first 15 yrs)
scenario_yrs_fun <- function(x) {
  x %>% filter(harvest_delay == "delay 15") %>% 
    group_by(index, production_target, original_habitat, habitat) %>%  
    count()
}

# Check our correction worked
count <- lapply(scenarios, scenario_yrs_fun)

# Get number of staggered harvests to define harvest window
J <- scenarios[[4]] 
harvest_window <- J$harvest_delay %>% unique %>% length()

# Get number of staggered harvests for once-logged to twice-logged transitions
harvest_window_short <- scenarios[[4]] %>% 
  filter(original_habitat == "once-logged" & habitat == "twice-logged") 
harvest_window_short <- harvest_window_short$harvest_delay %>% unique %>% length()

scenario_ACD_fun <- function(x) {
  # Convert to data.table for efficiency if not already
  if (!is.data.table(x)) x <- as.data.table(x)
  
  # Step 1: Calculate at 10km² parcel scale with delta method for error propagation
  x <- x %>%
    mutate(
      # For linear transformations (multiplication by constant),
      # delta method gives: Var[c*X] = c²*Var[X]
      # So SE[c*X] = c*SE[X]
      all_10km = full_carbon * 1000,
      all_10km_err = full_carbon_err * 1000,  
      
      ACD_10km = ACD * 1000,
      ACD_10km_err = ACD_err * 1000,
      
      # Determine divisor based on habitat transition
      harvest_divisor = if_else(original_habitat == "once-logged" & habitat == "twice-logged", 
                                harvest_window_short, harvest_window),
      
      # For division by constant, delta method gives: Var[X/c] = Var[X]/c²
      # For multiplication by constant: Var[c*X] = c²*Var[X]
      # Combined for X*n/d: Var[X*n/d] = (n/d)²*Var[X]
      # So SE[X*n/d] = (n/d)*SE[X]
      all_10km2_stag = all_10km * num_parcels / harvest_divisor,
      all_10km2_stag_err = all_10km_err * num_parcels / harvest_divisor,
      
      ACD_10km2_stag = ACD_10km * num_parcels / harvest_divisor,
      ACD_10km2_stag_err = ACD_10km_err * num_parcels / harvest_divisor
    )
  
  # Step 2: Group by year and habitat transition with delta method for sums
  habitat_year_summary <- x %>%
    group_by(index, production_target, original_habitat, habitat, true_year) %>%
    summarize(
      # Sum values across harvest schedules
      hab_all_year = sum(all_10km2_stag, na.rm = TRUE),
      
      # Delta method for sum of independent random variables:
      # Var[X₁ + X₂ + ... + Xₙ] = Var[X₁] + Var[X₂] + ... + Var[Xₙ]
      # So SE[X₁ + X₂ + ... + Xₙ] = sqrt(SE[X₁]² + SE[X₂]² + ... + SE[Xₙ]²)
      hab_all_year_err = sqrt(sum(all_10km2_stag_err^2, na.rm = TRUE)),
      
      hab_ACD_year = sum(ACD_10km2_stag, na.rm = TRUE),
      hab_ACD_year_err = sqrt(sum(ACD_10km2_stag_err^2, na.rm = TRUE)),
      
      # Keep necessary metadata
      scenarioName = first(scenarioName),
      scenarioStart = first(scenarioStart),
      .groups = "drop"
    )
  
  # Step 3: Group by year to calculate total carbon across all habitat transitions
  scenario_year_summary <- habitat_year_summary %>%
    group_by(index, production_target, true_year) %>%
    summarize(
      # Sum across all habitat transitions
      scen_all_year = sum(hab_all_year, na.rm = TRUE),
      
      # Apply delta method for sum of independent variables again
      scen_all_year_err = sqrt(sum(hab_all_year_err^2, na.rm = TRUE)),
      
      scen_ACD_year = sum(hab_ACD_year, na.rm = TRUE),
      scen_ACD_year_err = sqrt(sum(hab_ACD_year_err^2, na.rm = TRUE)),
      
      # Keep metadata
      scenarioName = first(scenarioName),
      scenarioStart = first(scenarioStart),
      .groups = "drop"
    )
  
  # Calculate confidence intervals based on standard errors
  # For normally distributed variables, 95% CI is approximately ±1.96*SE
  final_summary <- scenario_year_summary %>%
    mutate(
      # 68% confidence intervals (±1 SE)
      scen_all_year_lwr = scen_all_year - scen_all_year_err,
      scen_all_year_upr = scen_all_year + scen_all_year_err,
      
      scen_ACD_year_lwr = scen_ACD_year - scen_ACD_year_err,
      scen_ACD_year_upr = scen_ACD_year + scen_ACD_year_err,
      
      # # Optionally add 95% confidence intervals (±1.96 SE)
      # scen_all_year_lwr_95 = scen_all_year - 1.96 * scen_all_year_err,
      # scen_all_year_upr_95 = scen_all_year + 1.96 * scen_all_year_err,
      # 
      # scen_ACD_year_lwr_95 = scen_ACD_year - 1.96 * scen_ACD_year_err,
      # scen_ACD_year_upr_95 = scen_ACD_year + 1.96 * scen_ACD_year_err
    )
  
  return(final_summary)
}

# Apply scenario_ACD_fun to all scenarios
scenarios <- lapply(scenarios, scenario_ACD_fun)

###Determine carbon stock years -----------------------------------------

# Combine all scenarios for carbon stock calculation
carbon_stock_years <- rbindlist(scenarios) %>% 
  group_by(production_target, index, scenarioName, scenarioStart) %>%  
  reframe(
    all_carbon_stock = sum(scen_all_year, na.rm = TRUE), 
    # Proper error propagation for sum
    all_carbon_stock_err = sqrt(sum(scen_all_year_err^2, na.rm = TRUE)),
    all_carbon_stock_lwr = all_carbon_stock - all_carbon_stock_err,
    all_carbon_stock_upr = all_carbon_stock + all_carbon_stock_err,
    
    aboveground_carbon_stock = sum(scen_ACD_year, na.rm = TRUE),
    aboveground_carbon_stock_err = sqrt(sum(scen_ACD_year_err^2, na.rm = TRUE)),
    aboveground_carbon_stock_lwr = aboveground_carbon_stock - aboveground_carbon_stock_err,
    aboveground_carbon_stock_upr = aboveground_carbon_stock + aboveground_carbon_stock_err
  )


# Calculate starting landscape carbon ####
# Get carbon of starting landscape
starting_carbon <- hab_carbon %>% 
  filter(original_habitat == habitat) %>%  
  filter(!(original_habitat == "primary" & habitat == "primary")) %>%
  filter(!(original_habitat == "deforested" & habitat == "deforested")) 

# Calculate carbon per year per primary and deforested
starting_carbon_def <- hab_carbon %>%
  filter(original_habitat == "deforested" & habitat == "deforested") %>%
  uncount(61) %>% 
  # Add age 0-60 for primary
  mutate(functionalhabAge = row_number() - 1)

starting_carbon_prim <- hab_carbon %>%
  filter(original_habitat == "primary" & habitat == "primary") %>%
  uncount(61) %>% 
  # Add age 0-60 for primary
  mutate(functionalhabAge = row_number() - 1)

starting_carbon <- starting_carbon %>% 
  rbind(starting_carbon_def) %>%  
  rbind(starting_carbon_prim)  

# Calculate all-start landscape ACD for a given year  
all_start_landscapeCarbon <- all_start_landscape %>%  
  left_join(starting_carbon, relationship = "many-to-many") %>% 
  
  # Calculate ACD per 10km2 and for each habitat type with  error propagation
  mutate(
    # Include belowground processes with  error propagation
    all_SL = full_carbon * 1000 * num_parcels, 
    all_SL_err = full_carbon_err * 1000 * num_parcels,  # Error propagation for multiplication
    lwr_all_SL = all_SL - all_SL_err,
    upr_all_SL = all_SL + all_SL_err,
    
    # Aboveground processes only with proper error propagation
    ACD_SL = ACD * 1000 * num_parcels, 
    ACD_SL_err = ACD_err * 1000 * num_parcels,
    lwr_ACD_SL = ACD_SL - ACD_SL_err,
    upr_ACD_SL = ACD_SL + ACD_SL_err
  ) %>% 
  
  # Calculate ACD for a given year across habitat types with proper error propagation
  group_by(scenarioStart, functionalhabAge) %>% 
  summarize(
    # Sum with proper error propagation
    all_SL_tot = sum(all_SL, na.rm = TRUE),
    all_SL_tot_err = sqrt(sum(all_SL_err^2, na.rm = TRUE)),
    lwr_all_SL_tot = all_SL_tot - all_SL_tot_err,
    upr_all_SL_tot = all_SL_tot + all_SL_tot_err,
    
    ACD_SL_tot = sum(ACD_SL, na.rm = TRUE), 
    ACD_SL_tot_err = sqrt(sum(ACD_SL_err^2, na.rm = TRUE)),
    lwr_ACD_SL_tot = ACD_SL_tot - ACD_SL_tot_err,
    upr_ACD_SL_tot = ACD_SL_tot + ACD_SL_tot_err,
    .groups = "drop"
  ) %>%
  rename(true_year = functionalhabAge)

# Select relevant columns for joining
all_start_landscape_totalACD <- all_start_landscapeCarbon %>% 
  select(scenarioStart, all_SL_tot, all_SL_tot_err, lwr_all_SL_tot, upr_all_SL_tot,
         ACD_SL_tot, ACD_SL_tot_err, lwr_ACD_SL_tot, upr_ACD_SL_tot, true_year)

# Add starting landscape ACD to each scenario
add_SL_ACD_fun <- function(x) {
  x %>% 
    left_join(all_start_landscape_totalACD, by = c("scenarioStart", "true_year")) %>% 
    distinct()
}

scenarios <- lapply(scenarios, add_SL_ACD_fun)

# Create combined dataframe of all scenarios for plotting
plot_data_df <- rbindlist(scenarios)
plot_data_df_05 <- plot_data_df %>% filter(production_target == 0.5)

# Plot ACD through time for each scenario
ggplot(plot_data_df_05, aes(true_year, scen_ACD_year/1000000, group = index, colour = index)) +
  geom_point() +
  geom_line(aes(y = ACD_SL_tot/1000000), linetype = "longdash", colour = "black") +
  # Add error bands
  geom_ribbon(aes(ymin = scen_ACD_year_lwr/1000000, ymax = scen_ACD_year_upr/1000000, 
                  fill = index), alpha = 0.2, colour = NA) +
  facet_wrap(~scenarioName) +
  theme_bw(base_size = 16) + 
  theme(strip.background = element_blank(),
        legend.position = "none",
        strip.text = element_text(hjust=0, face="bold"), 
        panel.grid = element_blank(), 
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle=45, hjust=1.05, colour="black")) +
  labs(y = "ACD in millions", x="")

# Function to calculate annual changes in ACD with proper error propagation
ACD_change_function <- function(x){
  x %>%  
    #calculate change in ACD annually for the scenario landscape
    group_by(index, production_target) %>%  
    arrange(true_year, by_group = TRUE) %>% 
    
    #all carbon 
    mutate(scen_all_change = scen_all_year - lag(scen_all_year), 
           scen_all_change_lwr = scen_all_year_lwr - lag(scen_all_year_lwr),
           scen_all_change_upr = scen_all_year_upr- lag(scen_all_year_upr), 
           
           #Aboveground       
           scen_ACD_change = scen_ACD_year - lag(scen_ACD_year), 
           scen_ACD_change_lwr = scen_ACD_year_lwr - lag(scen_ACD_year_lwr),
           scen_ACD_change_upr = scen_ACD_year_upr- lag(scen_ACD_year_upr)
    ) %>% 
    
    #calculate ACD change annually for the starting landscape 
    mutate(
      SL_all_change = all_SL_tot- lag(all_SL_tot), 
      SL_all_change_lwr =lwr_all_SL_tot- lag(lwr_all_SL_tot),
      SL_all_change_upr = upr_all_SL_tot- lag(upr_all_SL_tot),
      
      SL_ACD_change = ACD_SL_tot- lag(ACD_SL_tot), 
      SL_ACD_change_lwr =lwr_ACD_SL_tot- lag(lwr_ACD_SL_tot),
      SL_ACD_change_upr = upr_ACD_SL_tot- lag(upr_ACD_SL_tot)) %>% 
    ungroup() %>%  
    select(index, production_target,true_year, 
           # scen_ACD_year,scen_ACD_year_lwr,scen_ACD_year_upr,
           scen_all_change,scen_all_change_lwr,scen_all_change_upr,
           scen_ACD_change,scen_ACD_change_lwr,scen_ACD_change_upr,
           # ACD_SL_tot, lwr_ACD_SL_tot,upr_ACD_SL_tot,
           SL_all_change,SL_all_change_lwr,SL_all_change_upr,
           SL_ACD_change,SL_ACD_change_lwr,SL_ACD_change_upr,
           scenarioName,scenarioStart)
  
  
}

scenarios <-lapply(scenarios, ACD_change_function)

scenarios_df <- rbindlist(scenarios)

# Plot scenario ACD annual change
scenarios_df %>%
  filter(production_target == 0.5) %>% 
  filter(!(true_year > 43 & true_year < 47)) %>% 
  filter(scenarioName %in% c("all_primary_CY_D.csv", "mostly_1L_CY_D.csv", "mostly_2L_CY_D.csv")) %>% 
  ggplot(aes(true_year, scen_ACD_change/1000000, colour = as.factor(index))) +
  geom_line() +
  # Add error bands
  geom_ribbon(aes(ymin = scen_ACD_change_lwr/1000000, ymax = scen_ACD_change_upr/1000000, 
                  fill = as.factor(index)), alpha = 0.2, colour = NA) +
  geom_line(aes(y = SL_ACD_change/1000000), linetype = "longdash", colour = "black") +
  facet_wrap(~scenarioName, ncol = 4, scales = "free_y") +
  theme_bw(base_size = 16) + 
  theme(strip.background = element_blank(), 
        legend.position = "none",
        strip.text = element_text(hjust=0, face="bold"), 
        panel.grid = element_blank(), 
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle=45, hjust=1.05, colour="black")) +
  labs(y = "ACD change in millions", x="")

# # Plot changes in ACD in starting landscape
# scenarios_df %>% 
#   ungroup() %>% 
#   filter(production_target == 0.5) %>% 
#   select(true_year, SL_ACD_change, SL_ACD_change_lwr, SL_ACD_change_upr, scenarioName) %>% 
#   group_by(true_year, scenarioName) %>%
#   slice(1) %>% 
#   ungroup() %>% 
#   ggplot(aes(true_year, SL_ACD_change)) +
#   geom_line() +
#   # Add error bands
#   geom_ribbon(aes(ymin = SL_ACD_change_lwr, ymax = SL_ACD_change_upr), 
#               alpha = 0.2, fill = "grey") +
#   facet_wrap(~scenarioName) +
#   theme_bw() + 
#   theme(strip.background = element_blank(), 
#         legend.position = "none",
#         strip.text = element_text(hjust=0, face="bold"), 
#         panel.grid = element_blank(), 
#         axis.text = element_text(colour = "black"),
#         axis.text.x = element_text(angle=45, hjust=1.05, colour="black")) +
#   labs(y = "drawdown", x="")

# Define molecular weight (mw) conversion for CO2
mw <- 44.01/12.01

flux_converted_function <- function(x){
  x %>% 
    mutate(
      
      
      #anual fluxes in scenario
      scen_flux_all = scen_all_change *mw, 
      scen_flux_all_lwr = scen_all_change_lwr *mw,
      scen_flux_all_upr = scen_all_change_upr *mw,
      
      scen_flux_ACD = scen_ACD_change *mw, 
      scen_flux_ACD_lwr = scen_ACD_change_lwr *mw,
      scen_flux_ACD_upr = scen_ACD_change_upr *mw,
      
      #anual fluxes in starting landscapes 
      
      SL_flux_all = SL_all_change *mw, 
      SL_flux_all_lwr =SL_all_change_lwr *mw,
      SL_flux_all_upr = SL_all_change_upr * mw, 
      
      SL_flux_ACD = SL_ACD_change *mw, 
      SL_flux_ACD_lwr =SL_ACD_change_lwr *mw,
      SL_flux_ACD_upr = SL_ACD_change_upr * mw
      
    )
  
  
  
}

#!!!! This allows running a short test for a single scenario 
#Single scenario test 
# scenario_flux <- flux_converted_function(scenarios)
# flux_per_scenario_and_SL <- scenario_flux %>%
#   ggplot( aes(x = true_year, y = scen_flux/1000000)) +
#   geom_line(colour = "red")+
#   geom_line(aes(x = true_year , y = SL_flux/1000000),linetype = "dashed") +
#   ylab("Scenario Fluxes in millions
#        (involes converting carbon to CO2 equivalents)")+
#   theme_bw(base_size = 14)
#!!!!

scenario_flux <- lapply(scenarios,flux_converted_function)
scenario_flux_df <- rbindlist(scenario_flux)

scenario_flux_df <- rbindlist(scenario_flux)

#------ incorporate social cost of carbon adjust discount rate -------

#CONSIDER FOR A 4 % DISCOUNT RATE (CAN FILTER TO GET 2 OR 6)
socialDR_4 <- socialDR %>% filter(discount_rate == "4%") %>%
  rename(true_year = year) 

socialDR_2 <- socialDR %>% filter(discount_rate == "2%") %>%
  rename(true_year = year) 

socialDR_6 <- socialDR %>% filter(discount_rate == "6%") %>%
  rename(true_year = year) 


socialDR_fun <- function(x, SDR){
  
  x %>% left_join(SDR, by = "true_year") %>% 
    
    #apply final equation to calculte socially discounted carbon impacts 
    mutate(carbon_impact_all = (scen_flux_all - SL_flux_all) * scc_discounted, 
           carbon_impact_all_lwr = (scen_flux_all_lwr - SL_flux_all_lwr) * scc_discounted, 
           carbon_impact_all_upr = (scen_flux_all_upr - SL_flux_all_upr) * scc_discounted,
           
           carbon_impact_ACD = (scen_flux_ACD - SL_flux_ACD) * scc_discounted, 
           carbon_impact_ACD_lwr = (scen_flux_ACD_lwr - SL_flux_ACD_lwr) * scc_discounted, 
           carbon_impact_ACD_upr = (scen_flux_ACD_upr - SL_flux_ACD_upr) * scc_discounted) %>%   
    
    #take the sum to give one final value for the whole scenario 
    group_by(index, production_target) %>%  
    
    mutate(TOTcarbon_all_impact = sum(carbon_impact_all,na.rm = TRUE), 
           TOTcarbon_impact_all_lwr = sum(carbon_impact_all_lwr,na.rm = TRUE), 
           TOTcarbon_impact_all_upr = sum(carbon_impact_all_upr,na.rm = TRUE),
           
           TOTcarbon_ACD_impact = sum(carbon_impact_ACD,na.rm = TRUE), 
           TOTcarbon_impact_ACD_lwr = sum(carbon_impact_ACD_lwr,na.rm = TRUE), 
           TOTcarbon_impact_ACD_upr = sum(carbon_impact_ACD_upr,na.rm = TRUE)
    ) %>%  
    
    ungroup()
  
}


# #---------Sense check ------- 
# As a sense check, run two equation formats (Tom Swinfield vs Andrew Balmford) to check they are the same
#good their identical 
# print(single) 
# 
# single <- single %>% left_join(socialDR_2, by = "true_year")
# 
# #select relevant columns only 
# single <- single %>% select(index, true_year, scen_flux, SL_flux, scc_discounted, discount_rate)
# 
# 
# #Andrew's way #-22135273
# andrew <- single %>%  
#   mutate(
#     #get discounted annual impact for scenario and starting landscape
#     SC_carbon_impact = (scen_flux * scc_discounted),
#     SL_carbon_impact = (SL_flux * scc_discounted), 
#     
#     #sum over these values
#     SC_summedC = sum(SC_carbon_impact,na.rm = TRUE), 
#     SL_summedC = sum(SL_carbon_impact,na.rm = TRUE), 
#     
#     
#     #minus SL from Scenario to get full impact 
#     totalSCC = SC_summedC - SL_summedC
#   )
# 
# 
# #Swinfield's way #-22135273
# 
# Swinfield <- single %>%  
#   mutate(
#     #get discounted annual impact for scenario and starting landscape
#     carbon_impact = (scen_flux - SL_flux) * scc_discounted, 
#     totalSCC = sum(carbon_impact,na.rm = TRUE))


##!!!
# #carbon impact = 5,129,924,965   -5129924965
# final_carbon2 <- socialDR_fun(scenario_flux,socialDR_2) %>%
#   mutate(netFlux = scen_flux - SL_flux) %>%
#   select(index, production_target, true_year,scen_flux, SL_flux,netFlux, scc_discounted, discount_rate,carbon_impact, TOTcarbon_impact)
# #carbon impact = = - 4,916,152,909
# final_carbon4 <- socialDR_fun(scenario_flux,socialDR_4) %>%
#   mutate(netFlux = scen_flux -  SL_flux) %>%
#   select(index, production_target, true_year,scen_flux, SL_flux,netFlux, scc_discounted, discount_rate,carbon_impact, TOTcarbon_impact)
# #carbon impact = 4,307,306,658
# final_carbon6 <- socialDR_fun(scenario_flux,socialDR_6) %>%
#   mutate(netFlux = scen_flux - SL_flux) %>%
#   select(index, production_target, true_year,scen_flux, SL_flux,netFlux, scc_discounted, discount_rate,carbon_impact, TOTcarbon_impact)
# final_carbon_comb <- rbind(final_carbon2,final_carbon4,final_carbon6)
# ##!!!!!

#Net flux #### 

final_carbon2_fun <- function(x){
  socialDR_fun(x,socialDR_2) %>% 
    mutate(netFlux_all = SL_flux_all -scen_flux_all, 
           netFlux_ACD = SL_flux_ACD -scen_flux_ACD
    ) #%>% 
  #    select(index, scenarioName,production_target, true_year,scen_flux, SL_flux,netFlux, scc_discounted, discount_rate,carbon_impact, TOTcarbon_impact)
}
final_carbon4_fun <- function(x){
  socialDR_fun(x,socialDR_4) %>% 
    mutate(netFlux_all = SL_flux_all -scen_flux_all, 
           netFlux_ACD = SL_flux_ACD -scen_flux_ACD 
    )# %>% 
  #select(index,scenarioName, production_target, true_year,scen_flux, SL_flux,netFlux, scc_discounted, discount_rate,carbon_impact, TOTcarbon_impact)
}
final_carbon6_fun <- function(x){
  socialDR_fun(x,socialDR_6) %>% 
    mutate(netFlux_all = SL_flux_all -scen_flux_all, 
           netFlux_ACD = SL_flux_ACD -scen_flux_ACD 
    )# %>% 
  #select(index,scenarioName, production_target, true_year,scen_flux, SL_flux,netFlux, scc_discounted, discount_rate,carbon_impact, TOTcarbon_impact)
}

final_carbon2 <- lapply(scenario_flux,final_carbon2_fun)
final_carbon4 <- lapply(scenario_flux,final_carbon4_fun)
final_carbon6 <- lapply(scenario_flux,final_carbon6_fun)

final_carbon2 <- rbindlist(final_carbon2) 
final_carbon4 <- rbindlist(final_carbon4) 
final_carbon6 <- rbindlist(final_carbon6)
final_carbon_comb <- final_carbon2 %>%  rbind(final_carbon4) %>% rbind(final_carbon6)

# Plot netFlux vs true_year for different discount rates


# Plot carbon_impact vs true_year for different discount rates
carbonImpact <- final_carbon_comb %>% drop_na %>%  
  
  group_by(production_target,index) %>% 
  filter(production_target == 0.3) %>% 
  filter(scenarioName %in% c("all_primary_CY_D.csv", "mostly_1L_CY_D.csv", "mostly_2L_CY_D.csv")) %>% 
  
  ggplot( aes(x = true_year, y = carbon_impact_all/10000000, color = index)) +
  geom_line() +
  # Add shaded regions
  annotate("rect", xmin = min(final_carbon_comb$true_year), xmax = max(final_carbon_comb$true_year), ymin = 0, ymax = Inf, fill = "honeydew") +
  annotate("rect", xmin = min(final_carbon_comb$true_year), xmax = max(final_carbon_comb$true_year), ymin = -Inf, ymax = 0, fill = "mistyrose" ) +
  # Add horizontal line at y = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line() +
  facet_grid(scenarioName ~discount_rate)+
  labs(title = "Carbon impact 
  (ScenarioFlux - StartingLandscapeFlux) * scc_discounted)
       for different discount rates",
       x = "True Year",
       y = "Carbon Impact ^6",
       color = "Discount Rate") +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, hjust = 0.5)) 

#sFlux <- scenario_flux%>%
# sFlux <- scenario_flux_df%>%
#   filter(production_target == 0.3) %>% 
#   
#   drop_na %>% 
#   
#   group_by(production_target,index) %>% 
#   filter(scenarioName %in% c("all_primary_CY_D.csv", "mostly_1L_CY_D.csv", "mostly_2L_CY_D.csv")) %>% 
#   
#   mutate(net_flux =  scen_flux/1000000 - SL_flux/1000000 ) %>% 
#   ggplot( aes(x = true_year, y = net_flux, colour = index)) +
#   geom_line()+
#   geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
#   
#   facet_wrap(~scenarioName)+
#   
#   ylab("Fluxes (ScenarioFlux^6 - StartingLandscapeFlux^6)")+
#   theme_minimal(base_size = 14)+
#   theme(legend.position = "none", 
#         axis.text = element_text(size = 12), 
#         axis.title = element_text(size = 14, face = "bold"),
#         plot.title = element_text(size = 16, hjust = 0.5)) 


scc <- final_carbon_comb %>% select(true_year,discount_rate, scc_discounted) %>% 
  drop_na %>%  
  ggplot( aes(x = true_year, y = scc_discounted, color = discount_rate)) +
  geom_line() +
  facet_wrap(~discount_rate)+
  ylab("Discounted social cost of carbon")+
  theme_bw()+
  theme(legend.position = "none", 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, hjust = 0.5))


##!!!
#DISCOUNT RATE 4%
final_carbon_4DR <- lapply(scenario_flux, socialDR_fun, socialDR_4)
final_carbon_df_4DR <- rbindlist(final_carbon_4DR) 

#DISCOUNT RATE 2%
final_carbon_2DR <- lapply(scenario_flux, socialDR_fun, socialDR_2)
final_carbon_df_2DR <- rbindlist(final_carbon_2DR)

#DISCOUNT RATE 6%
final_carbon_6DR <- lapply(scenario_flux, socialDR_fun, socialDR_6)
final_carbon_df_6DR <- rbindlist(final_carbon_6DR)

all_discount_rates <- final_carbon_df_2DR %>% rbind(final_carbon_df_4DR) %>% rbind(final_carbon_df_6DR)

#-----EXPORT OUTCOME PERFORMANCE for consolidated figure of all outcomes -----
getwd()
names(all_discount_rates)
output <- all_discount_rates %>%
  #remove year 0 as this will have NA values 
  filter(true_year> 0 ) %>% 
  select(index, production_target, scenarioName,scenarioStart,
         #all carbon (i.e. incorporates belowground carbon)
         TOTcarbon_all_impact,TOTcarbon_impact_all_lwr, TOTcarbon_impact_all_upr,
         #aboveground carbon only
         TOTcarbon_ACD_impact,TOTcarbon_impact_ACD_lwr, TOTcarbon_impact_ACD_upr,
         
         discount_rate) %>% 
  unique() %>% 
  cbind(outcome = "carbon")
output <- unique(output)

#save outputs ####
saveRDS(output, "FinalPerformanceOutput/MasterCarbonPerformance_withuncertainty.rds")

carbon_stock_years %>% cbind(outcome = "carbon")
saveRDS(carbon_stock_years, "FinalPerformanceOutput/carbonstock_years__withuncertainty.rds")
