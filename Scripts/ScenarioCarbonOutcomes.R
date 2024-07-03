# 09.06.24; 
#NB; this code:
#this code calculates carbon consequences of different scenarios
#commented out code-blocks enable to assessment of a single scenario, to see how the code works.


#!!27.06.24
#Need to check because there's a 40 yr spike in carbon that's happening for some reason 

library(tidyr)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)
library(stringr) 
library(cowplot)

#READ!!!: there are sign-posted hard-coded param decisions in this script!!

#Hard-coded param decisions: 
#Reading in carbon by hab outcomes from  CalculateAllHabCarbonVals.R script, 
#beware that ACD refers to just above-ground carbon, whereas
#all_carbon incorporates belowground/necromass processes. 


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
hab_carbon <-read.csv("Outputs/allHabCarbon_60yrACD.csv") %>% select(-X)

#--- run pipeline for a single scenario -------
# # 
 J <- scenarios[[1]]
 test_scen<-  J %>% filter(index == "all_primary_CY_D.csv 253") 


#convert hab_harbon to data-table 
hab_carbon <- as.data.table(hab_carbon)

#carry out scenario corrections ####

#1. Correct for the fact that habitat starting as 1L cannot 
#be harvested in the first 15 years
#2. For efficiency in joins, if functionalHab = primary or deforested, make functionalhabAge = 0.
#I changed the hab carbon for these fix values to only give 
#values for 0 functionalHabAge

adjust_scenarios <- function(data) {
 
   data %>%
    mutate(
      functionalhabAge = if_else(functional_habitat %in% c("primary", "deforested"), 0, functionalhabAge),
      harvest_delay_numeric = as.numeric(stringr::str_extract(harvest_delay, "\\d+"))
    ) %>%
    filter(!(original_habitat == "once-logged" & habitat == "twice-logged" & harvest_delay_numeric > 15)) %>%
    select(-harvest_delay_numeric)
}

scenarios <- lapply(scenarios, adjust_scenarios)

#!!!!
#single scenario 
test_scen <- adjust_scenarios(test_scen)
#!!!!

#join carbon information to scenarios #### 

carbon_fun <- function(x){
  x <- as.data.table(x)
  # Perform the left join operation in data table 
  result <- x[hab_carbon, on = .(habitat = habitat,
                                 original_habitat = original_habitat, 
                                 functional_habitat = functional_habitat, 
                                 functionalhabAge = functionalhabAge), nomatch = 0]
  
  #equivalent to: 
  # x %>%  left_join(hab_carbon, by = c("habitat" = "habitat",
  #                                     "original_habitat" = "original_habitat",
  #                                      "functional_habitat" = "functional_habitat", 
  #                                       "functionalhabAge" = "functionalhabAge") 
}
#!!!!
#single scenario 
test_scen <- carbon_fun(test_scen)
#!!!!

scenarios <- lapply(scenarios,carbon_fun)



# -------- ACD for a true year----------
#Make the hard coded decision below of whether to use all_carbon
#(includes belowground processes), or just above ground carbon (ACD)


##1. Take Carbon measures from ACD per hectare to 10km2 parcel levels
#2. For each staggered harvest, assume that 1/30th of each occurs (and so 1/30th annually)
#3 calculate hab ACD for a given year, across staggered harvests
#4 calculate scenario ACD for a given year 

#get number of staggered harvests to define harvest window (this must match harvests)
J <- scenarios[[4]] 
harvest_window <- J$harvest_delay %>% unique %>% length()

#get number of staggered harvests (fewer) for once-logged to twice-logged transitions
harvest_window_short <-   scenarios[[4]] %>% filter(original_habitat == "once-logged" & habitat == "twice-logged") 
harvest_window_short<- harvest_window_short$harvest_delay %>% unique %>% length()
#!!!!!!
#use for single scenario test
 harvest_window <- 30
#!!!!!!


scenario_ACD_fun <- function(x){
  x %>% 
    #1. go to 10km2 parcel scale for carbon
    
  #  HARD-CODED DECISION!!!!! #### 
  #include aboveground carbon only, or also belowground 

  #INCLUDE BELOWGROUND PROCESSES 
  mutate(ACD_10km = full_carbon*1000, 
         lwr_ACD_10km =full_carbon_lwr *1000, 
         upr_ACD_10km = full_carbon_upr*1000) %>% 
  
  # #ABOVEGROUND PROCESSES ONLY 
  #  mutate(ACD_10km = ACD*1000, 
  #          lwr_ACD_10km =lwr_ACD *1000, 
  #          upr_ACD_10km = upr_ACD*1000) %>% 
    
    #2. assuming 1/30th of of each habitat type is applied to each harvesting delay schedule
    #, calculate the total ACD for a given habitat type in a given year
    #NB- if there is no habitat transition, then don't need to divide by harvest window 
    mutate(
      ACD_10km2_stag =  ACD_10km * num_parcels / harvest_window,
      lwr_ACD_10km2_stag = lwr_ACD_10km * num_parcels / harvest_window,
      upr_ACD_10km2_stag = upr_ACD_10km * num_parcels / harvest_window)  %>% 
    

  #3. for each true year and habitat transition, calculate carbon combined across the staggered
  #harvesting schedule (i.e. the carbon in a given habitat transition for a given year) 
  group_by(index,production_target, original_habitat, habitat, true_year) %>%  
    mutate(hab_ACD_year = sum(ACD_10km2_stag), 
           hab_ACD_year_lwr= sum(lwr_ACD_10km2_stag), 
           hab_ACD_year_upr = sum(upr_ACD_10km2_stag)) %>%  ungroup %>%  
    
    #select a single harvest delay worth of data, as we now have calculated ACD across harvesting schedules
    filter(harvest_delay == "delay 15") %>% 
    select(-harvest_delay) %>% 
    
    
    
    #4. Across habitat type transitions (e.g for all hab_parcel transitions) in a scenario, calculate ACD for a given year
    group_by(index, production_target, true_year) %>%  
    mutate(scen_ACD_year = sum(hab_ACD_year), 
           scen_ACD_year_lwr= sum(hab_ACD_year_lwr), 
           scen_ACD_year_upr = sum(hab_ACD_year_upr)) %>%  ungroup() %>% 
    
    
    #now we make sure we only have one row for each scenario and year, showing scen_ACD_year
    select(index, production_target,scenarioName,scenarioStart, true_year, 
           scen_ACD_year,scen_ACD_year_lwr, scen_ACD_year_upr) %>%  
    group_by(true_year,index) %>%  slice(1) %>% 
    ungroup()
}

# #!!!!
#Single scenario 
test_scen <- scenario_ACD_fun(test_scen)
#

scenarios <- lapply(scenarios, scenario_ACD_fun)
#JJ <- scenarios[[11]]
#test <- JJ %>%  group_by(index, true_year) %>% count()

#CONTUNINEU FROM HERE ####
#---------- calculate carbon in starting landscape ------------

#get carbon of starting landscape
starting_carbon <- hab_carbon %>% 
  filter(original_habitat == habitat) %>%  
  filter(!(original_habitat == "primary" & habitat == "primary")) %>%
  filter(!(original_habitat == "deforested" & habitat == "deforested")) 

#calculate carbon per year per primary and deforested  
starting_carbon_def <- hab_carbon %>%
  filter(original_habitat == "deforested" & habitat == "deforested") %>%
  uncount(61) %>% 
  #add age 0-60 for primary
  mutate(functionalhabAge = row_number() - 1)

starting_carbon_prim <- hab_carbon %>%
  filter(original_habitat == "primary" & habitat == "primary") %>%
  uncount(61) %>% 
  #add age 0-60 for primary
  mutate(functionalhabAge = row_number() - 1)

starting_carbon <- starting_carbon %>% 
  rbind(starting_carbon_def) %>%  
  rbind(starting_carbon_prim)  


#----- calculate all-start landscape ACD for a given year ------  
all_start_landscapeCarbon <- all_start_landscape %>%  
  left_join(starting_carbon, relationship = "many-to-many") %>% 
  
  #calculate ACD per 10km2 and for each habitat type
  
  #  HARD-CODED DECISION!!!!! #### 
#include aboveground carbon only, or also belowground 
  
  #above and belowground carbon 
   mutate(ACD_SL =  full_carbon* 1000 *num_parcels, 
       lwr_ACD_SL = full_carbon_lwr * 1000 *num_parcels,
       upr_ACD_SL = full_carbon_upr * 1000 *num_parcels) %>% 
  

  # #aboveground only - [uncomment]    
  # mutate(ACD_SL = ACD * 1000 *num_parcels, 
  #        lwr_ACD_SL = lwr_ACD * 1000 *num_parcels,
  #        upr_ACD_SL = upr_ACD * 1000 *num_parcels) %>% 
  
  
  ungroup %>% 
  #calcualte ACD for a given year across habitat types (e.g. across ALL parcels) 
  group_by(scenarioStart, functionalhabAge) %>% 
  mutate(ACD_SL_tot =sum(ACD_SL), 
         lwr_ACD_SL_tot = sum(lwr_ACD_SL),
         upr_ACD_SL_tot = sum(upr_ACD_SL)) %>% 
  #we have now summarise across multiple hab transitions - make sure we only have 1 
  group_by(scenarioStart, functionalhabAge) %>%  
  slice(1) %>% 
  ungroup() %>% 
  rename(true_year = functionalhabAge)



all_start_landscape_totalACD <- all_start_landscapeCarbon %>% select(
  scenarioStart, ACD_SL_tot,lwr_ACD_SL_tot,upr_ACD_SL_tot,true_year)

# --------- #add starting landscape ACD to each scenario ---------

add_SL_ACD_fun <- function(x){
  x %>% left_join(all_start_landscape_totalACD, by = c("scenarioStart","true_year")) %>% 
    distinct()
  
}
#!!!!!!
#single scenario
#scenarios <- add_SL_ACD_fun(scenarios)
#!!!!!!

scenarios <- lapply(scenarios, add_SL_ACD_fun)


# -------- plot ACD thru time for each scenario  ------------


# plot_data <- lapply(scenarios, get_plot_data_fun)
#plot_data_df <- rbindlist(plot_data)

plot_data_df <- rbindlist(scenarios)

plot_data_df_05 <- plot_data_df %>% filter(production_target == 0.5)


plot_data_df_05 %>% 
  # scenarios %>% 
  # filter(scenarioName == "all_primary_CY_D.csv") %>% 
  #!!!
  #scenarios %>% 
  #!!!
 # filter(scenarioName %in% c("all_primary_CY_D.csv", "mostly_1L_CY_D.csv", "mostly_2L_CY_D.csv")) %>% 
  
  ggplot(aes(true_year, scen_ACD_year/1000000, group = index,colour = index)) +
  geom_line() +
  geom_line(aes(y = ACD_SL_tot/1000000),linetype = "longdash", colour = "black") +
  #geom_point() +
  #geom_linerange() + 
  facet_wrap(~scenarioName) +
  theme_bw(base_size = 16) + 
  theme(strip.background = element_blank(),
        legend.position = "none",
        strip.text = element_text(hjust=0, face="bold"), 
        panel.grid = element_blank(), 
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle=45, hjust=1.05, colour="black")) +
  labs(y = "ACD in millions", x="")

#fig_scale <- 1.5
# ggsave("figures/occupancy_estimates.png", units="mm", 
#        height=120*fig_scale, width=115*fig_scale)

##WHAT'S WITH THE 40 YR DIP?
#something is happening at yr 45 for primary - twice-logged transitions, 
#for a single year, there is a big dip. 

spike_scenarios <- plot_data_df %>% group_by(scenarioStart) %>% 
filter(production_target == 0.5) %>%  
  filter(true_year>40) %>%  
   filter(scen_ACD_year == min(scen_ACD_year, na.rm = TRUE)) %>%  
    unique() %>% 
  ungroup() %>%
  select(index)

problem_compositions <- scenarios %>% rbindlist() %>%  
  filter(production_target == 0.5) %>%  
  right_join(scenario_composition, relationship = "many-to-many") %>% 
  right_join(spike_scenarios) %>%  
  mutate(scen_ACD_year = scen_ACD_year/1000000)

# --------  Convert ACD into changes in ACD per scenario ----------
names(scenarios)
ACD_change_function <- function(x){
  x %>%  
    # #filter only a single delay year(as scenario and SL ACD are now summarised across delay years, so these rows are identical)
    # select(index,production_target,true_year, original_habitat, habitat, scen_ACD_year,scen_ACD_year_lwr,scen_ACD_year_upr,
    #        ACD_SL_tot,lwr_ACD_SL_tot,upr_ACD_SL_tot,scenarioName) %>% 
    # group_by(index, production_target, true_year) %>% 
    # slice(1) %>%  
    # ungroup %>% 
    
    #calculate change in ACD annually for the scenario landscape
    group_by(index, production_target) %>%  
    arrange(true_year, by_group = TRUE) %>% 
    
    mutate(scen_ACD_change = scen_ACD_year - lag(scen_ACD_year), 
           scen_ACD_change_lwr = scen_ACD_year_lwr - lag(scen_ACD_year_lwr),
           scen_ACD_change_upr = scen_ACD_year_upr- lag(scen_ACD_year_upr)) %>% 
    
    #calculate ACD change annually for the starting landscape 
    mutate(SL_ACD_change = ACD_SL_tot- lag(ACD_SL_tot), 
           SL_ACD_change_lwr =lwr_ACD_SL_tot- lag(lwr_ACD_SL_tot),
           SL_ACD_change_upr = upr_ACD_SL_tot- lag(upr_ACD_SL_tot)) %>% 
    ungroup() %>%  
    select(index, production_target,true_year, 
           scen_ACD_year,scen_ACD_year_lwr,scen_ACD_year_upr,
           scen_ACD_change,scen_ACD_change_lwr,scen_ACD_change_upr,
           ACD_SL_tot, lwr_ACD_SL_tot,upr_ACD_SL_tot,
           SL_ACD_change,SL_ACD_change_lwr,SL_ACD_change_upr,
           scenarioName,scenarioStart)
  
  
}

#!!!!!!
#scenarios <- ACD_change_function(scenarios)
#!!!!!!

scenarios <-lapply(scenarios, ACD_change_function)
scenarios_df <- rbindlist(scenarios)


#plot scenario ACD annual change
scenarios_df %>%
  
  filter(production_target == 0.5) %>% 
  
  filter(scenarioName %in% c("all_primary_CY_D.csv", "mostly_1L_CY_D.csv", "mostly_2L_CY_D.csv")) %>% 
  #filter(scenarioName %in% c( "mostly_1L_CY_D.csv")) %>% 
  
  
  #!!!
  #scenarios %>% 
  #!!!
  ggplot(aes(true_year, scen_ACD_change/1000000,colour = as.factor(index))) +
  geom_line()+
  geom_line(aes(y = SL_ACD_change/1000000),linetype = "longdash", colour = "black") +
  #geom_point() +
  #geom_linerange() + 
  facet_wrap(~scenarioName, ncol =4, scales ="free_y" ) +
  theme_bw(base_size = 16) + 
  theme(strip.background = element_blank(), 
        legend.position = "none",
        strip.text = element_text(hjust=0, face="bold"), 
        panel.grid = element_blank(), 
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle=45, hjust=1.05, colour="black")) +
  labs(y = "ACD change in millions", x="")
#J <- scenarios[[1]] %>%  filter(index == "all_primary_CY_D.csv 251" | index =="all_primary_CY_D.csv 252"|index =="all_primary_CY_D.csv 253"|index =="all_primary_CY_D.csv 254"|index =="all_primary_CY_D.csv 255")
#Tx <- scenario_composition %>% filter(index == "all_primary_CY_D.csv 251" | index =="all_primary_CY_D.csv 252"|index =="all_primary_CY_D.csv 253"|index =="all_primary_CY_D.csv 254"|index =="all_primary_CY_D.csv 255")


# ------ change in ACD in starting landscape  #--------------
#good - starting landscapes with deforested land at beginning have lower ACD change, as have a smaller amount of recovering forest

scenarios_df %>% ungroup %>% filter(production_target ==0.5) %>% 
  select(true_year,SL_ACD_change,scenarioName) %>% group_by(true_year,scenarioName) %>%
  slice(1) %>% ungroup %>% 
  #!!!
  # scenarios %>% 
  #!!!
  ggplot(aes(true_year, SL_ACD_change)) +
  geom_line()+
  #geom_point() +
  #geom_linerange() + 
  facet_wrap(~scenarioName) +
  theme_bw() + 
  theme(strip.background = element_blank(), 
        legend.position = "none",
        strip.text = element_text(hjust=0, face="bold"), 
        panel.grid = element_blank(), 
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle=45, hjust=1.05, colour="black")) +
  labs(y = "drawdown", x="")


#NOTES ON VIEWING: filter a single iddex and then order inview finder by original_hab,habitat and true_year
#xx <- scenarios[[17]]
#xxx <-  xx %>% filter(index=="mostly_2L_deforested_IY_D.csv 5") %>%  distinct()
#View(xxx)

#---------- convert to fluxes across habitat parcels --------
#The atomic mass of carbon (C) is approximately 12.01 g/mol, and the molecular weight of carbon dioxide (CO2) is the sum of the atomic masses of one carbon atom and two oxygen (O) atoms, which is approximately 12.01 + 2 * 16.00 = 44.01 g/mol.
#Convert Carbon Change to CO2 Change:
#To convert the net change in above-ground carbon density to carbon dioxide flux, you need to multiply the carbon change by the ratio of the molecular weights:
#  CO2 Flux = Net Change in AGCD * (44.01 g/mol / 12.01 g/mol)

#define molecular weight (mw) conversion
mw <- 44.01/12.01

flux_converted_function <- function(x){
  x %>% 
    mutate(
      
      
      #anual fluxes in scenario
      scen_flux = scen_ACD_change *mw, 
      scen_flux_lwr = scen_ACD_change_lwr *mw,
      scen_flux_upr = scen_ACD_change_upr *mw,
      
      #anual fluxes in starting landscapes 
      
      SL_flux = SL_ACD_change *mw, 
      SL_flux_lwr =SL_ACD_change_lwr *mw,
      SL_flux_upr = SL_ACD_change_upr * mw)
  
  
}

#!!!!
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

#names(scenario_flux[[12]])
#x <- scenario_flux[[13]]
#single <- x %>% filter(index == "mostly_2L_CY_D.csv 10")


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
    mutate(carbon_impact = (scen_flux - SL_flux) * scc_discounted, 
           carbon_impact_lwr = (scen_flux_lwr - SL_flux_lwr) * scc_discounted, 
           carbon_impact_upr = (scen_flux_upr - SL_flux_upr) * scc_discounted) %>%  
    
    #take the sum to give one final value for the whole scenario 
    group_by(index, production_target) %>%  
    
    mutate(TOTcarbon_impact = sum(carbon_impact,na.rm = TRUE), 
           TOTcarbon_impact_lwr = sum(carbon_impact_lwr,na.rm = TRUE), 
           TOTcarbon_impact_upr = sum(carbon_impact_upr,na.rm = TRUE)) %>%  
    
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
    mutate(netFlux = SL_flux -scen_flux) %>% 
    select(index, scenarioName,production_target, true_year,scen_flux, SL_flux,netFlux, scc_discounted, discount_rate,carbon_impact, TOTcarbon_impact)
}
final_carbon4_fun <- function(x){
  socialDR_fun(x,socialDR_4) %>% 
    mutate(netFlux = SL_flux -scen_flux) %>% 
    select(index,scenarioName, production_target, true_year,scen_flux, SL_flux,netFlux, scc_discounted, discount_rate,carbon_impact, TOTcarbon_impact)
}
final_carbon6_fun <- function(x){
  socialDR_fun(x,socialDR_6) %>% 
    mutate(netFlux = SL_flux -scen_flux) %>% 
    select(index,scenarioName, production_target, true_year,scen_flux, SL_flux,netFlux, scc_discounted, discount_rate,carbon_impact, TOTcarbon_impact)
}

final_carbon2 <- lapply(scenario_flux,final_carbon2_fun)
final_carbon4 <- lapply(scenario_flux,final_carbon4_fun)
final_carbon6 <- lapply(scenario_flux,final_carbon6_fun)

final_carbon2 <- rbindlist(final_carbon2) 
final_carbon4 <- rbindlist(final_carbon4) 
final_carbon6 <- rbindlist(final_carbon6)
final_carbon_comb <- final_carbon2 %>%  rbind(final_carbon4) %>% rbind(final_carbon6)
# Plot netFlux vs true_year for different discount rates
# Define a color palette for the discount rates

# Plot netFlux vs true_year for different discount rates

# Plot carbon_impact vs true_year for different discount rates
carbonImpact <- final_carbon_comb %>% drop_na %>%  
  
  group_by(production_target,index) %>% 
  filter(production_target == 0.3) %>% 
  filter(scenarioName %in% c("all_primary_CY_D.csv", "mostly_1L_CY_D.csv", "mostly_2L_CY_D.csv")) %>% 
  
  ggplot( aes(x = true_year, y = carbon_impact/1000000, color = index)) +
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
                                        TOTcarbon_impact,TOTcarbon_impact_lwr, TOTcarbon_impact_upr, discount_rate) %>% 
  unique() %>% 
  cbind(outcome = "carbon")
output <- unique(output)

saveRDS(output, "Outputs/MasterCarbonPerformance.rds")

