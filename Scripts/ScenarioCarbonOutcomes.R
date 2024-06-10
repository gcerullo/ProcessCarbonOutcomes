# 09.06.24; 
#NB; this code:
#this code calculates carbon consequences of different scenarios 

library(tidyr)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)
library(stringr) 
library(cowplot)

#read in the scenario parametres containing conversion factors for converting from point to parcel/entire landscape  
source('R_code/BuildScenarios/BiolerCodeScenarioConstruction.R')

#read in social discount rate caluclate for 2,4,6% 
#built in the CalculateSocialDiscountRates.R script 

socialDR <- read.csv("Outputs/SocialDiscountRates_2_4_6pc_DR.csv")

#define all starting landscapes: 
all_start_landscape 


#----------------read in scenarios -------------------------------
#scenarios where all plantation conversion happens in year 0 
#scenarios <- readRDS("BuildScenarios/BuildingHarvestingScenarios/allScenarios.rds") 

#yield matched scenarios where 1/30th of plantation conversion happens annually
scenarios <- readRDS("R_code/BuildScenarios/BuildingHarvestingScenarios/allScenariosStaggered.rds")
scenario_composition <- rbindlist(scenarios, use.names=TRUE)

#------read in carbon by year per hab -----------------------
#read in carbon, with delays already calculated and with 
#below ground carbon and necromass incorporated
#WARNING; AS CURRENTLY IMPORTED, THE UNCERTAINTY (e.g. lwr and upr ACD) have not incoroprated belowground carbon dynamics)
# DO NOT USE OR TRUST THESE UNCERTAINTY VALUES. 
write.csv(allHabCarbon_60yrACD_withDelays, "R_code/CarbonAnalysis/Inputs/allHabCarbon_60yrACD_withDelays.csv")


hab_carbon <- read.csv("R_code/CarbonAnalysis/Inputs/allHabCarbon_60yrACD_withDelays.csv") %>% select(-X)

#quick checks 
x <- hab_carbon %>% group_by(original_habitat,habitat,harvest_delay) %>% count()
names(hab_carbon)
#--- run pipeline a single scenario -------
# # 
#  J <- rbindlist(scenarios,use.names=TRUE)
#  test_scen_comp <- scenario_composition %>% filter(index == "unq2_1793") # a scenario with a good mixture of logging
#  test_scen<-  J %>% filter(index == "all_primary_CY_D.csv 349") #a restored + primary scenario
# 
# ------Add temporal carbon data to scenairos -----------------
names(scenario_composition)
names(hab_carbon)

#convert hab_harbon to data-table 
hab_carbon <- as.data.table(hab_carbon)

carbon_fun <- function(x){
  x <- as.data.table(x)
  # Perform the left join operation in data table 
  result <- x[hab_carbon, on = .(habitat = habitat, original_habitat = original_habitat), allow.cartesian = TRUE, nomatch = 0]
  
  # x %>%  left_join(hab_carbon, by = c("habitat" = "habitat",
  #                                     "original_habitat" = "original_habitat"),
  #                                         relationship = "many-to-many")  
}
#!!!!
#scenarios <- carbon_fun(test_scen)
#!!!!
scenarios <- lapply(scenarios,carbon_fun)

# -------- CALCULATE SCENARIO ACD FOR A GIVEN TRUE YEAR ----------
##1. Take Carbon measures from ACD per hectare to 10km2 parcel levels
#2. For each staggered harvest, assume that 1/30th of each occurs (and so 1/30th annually)
#3 calculate hab ACD for a given year, across staggered harvests
#4 calculate scenario ACD for a given year 

#get number of staggered harvests to define harvest window (this must match harvests)
J <- scenarios[[19]] 
harvest_window <- J$harvest_delay %>% unique %>% length()

#use for single scenario test
# harvest_window <- 30

scenario_ACD_fun <- function(x){
  x %>% 
    #1. go to 10km2 parcel scale for carbon
    mutate(ACD_10km = ACD*1000, 
           lwr_ACD_10km =lwr_ACD *1000, 
           upr_ACD_10km = upr_ACD*1000) %>% 
    
    #2. assuming 1/30th of of each habitat type is applied to each harvesting delay schedule
    #, calculate the total ACD for a given habitat type in a given year
    #NB- if there is no habitat transition, then don't need to divide by harvest window 
    mutate(
      ACD_10km2_stag =  ACD_10km * num_parcels / harvest_window,
      lwr_ACD_10km2_stag = lwr_ACD_10km * num_parcels / harvest_window,
      upr_ACD_10km2_stag = upr_ACD_10km * num_parcels / harvest_window)  %>% 
    
    
    #!!!!####
  #Nb- to incorporate here - once-L to twice-L transitions are not allowed to occur in the first 15 yrs, as forests stating as 1L were harvested in yr t-1
  
  #3. for each true year and habitat transition, calculate ACD combined across the staggered
  #harvesting schedule (i.e. the ACD in a given habitat transition for a given year) 
  group_by(index,production_target, original_habitat, habitat, true_year) %>%  
    mutate(hab_ACD_year = sum(ACD_10km2_stag), 
           #----PROPOGATING ERROR INCORRECTLY?!!!!------
           hab_ACD_year_lwr= sum(lwr_ACD_10km2_stag), 
           hab_ACD_year_upr = sum(upr_ACD_10km2_stag)) %>%  ungroup %>%  
    
    #select a single harvest delay worth of data, as we now have calculated ACD across harvesting schedules
    filter(harvest_delay == 15) %>% select(-harvest_delay) %>% 
    
    
    
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
# scenTEST <- scenario_ACD_fun(scenarios)
# singleYear <- scenTEST %>% filter(true_year == 9)
# scenarios <- scenario_ACD_fun(scenarios)
#!!!!

scenarios <- lapply(scenarios, scenario_ACD_fun)
JJ <- scenarios[[11]]
test <- JJ %>%  group_by(index, true_year) %>% count()

#---------- calculate carbon in starting landscape ------------

#get ages of starting landscape thru time
time_forSL <- hab_carbon %>% filter(harvest_delay == 0) %>% 
  select(original_habitat, true_year) %>% unique() %>%  
  rename(habitat = original_habitat, 
         functionalhabAge = true_year)

#add hab ages to starting landscape 
all_start_landscape <- all_start_landscape %>%
  left_join(time_forSL, by = c("habitat"), relationship = "many-to-many")

#add carbon to starting landscape
starting_carbon <- hab_carbon %>% filter(harvest_delay == 0 & 
                                           original_habitat == habitat)

all_start_landscape <- all_start_landscape %>% 
  left_join(starting_carbon, by = c("habitat", "functionalhabAge")) 


#----- calculate all-start landscape ACD for a given year ------  
all_start_landscape <- all_start_landscape %>%  
  #calculate ACD per 10km2 and for each habitat type
  mutate(ACD_SL = ACD * 1000 *num_parcels, 
         lwr_ACD_SL = lwr_ACD * 1000 *num_parcels,
         upr_ACD_SL = upr_ACD * 1000 *num_parcels) %>%  ungroup %>% 
  #calcualte ACD for a given year across habitat types (e.g. across ALL parcels) 
  group_by(scenarioStart, functionalhabAge) %>% 
  mutate(ACD_SL_tot =sum(ACD_SL), 
         lwr_ACD_SL_tot = sum(lwr_ACD_SL),
         upr_ACD_SL_tot = sum(upr_ACD_SL)) %>% 
  #we have now summarise across multiple hab transitions - make sure we only have 1 
  group_by(scenarioStart, functionalhabAge) %>%  
  slice(1) %>% 
  ungroup() 



all_start_landscape_totalACD <- all_start_landscape %>% select(
  scenarioStart, ACD_SL_tot,lwr_ACD_SL_tot,upr_ACD_SL_tot,functionalhabAge) %>% 
  rename(true_year = functionalhabAge)


# --------- #add starting landscape ACD to each scenario ---------

add_SL_ACD_fun <- function(x){
  x %>% left_join(all_start_landscape_totalACD, by = c("scenarioStart","true_year")) %>% 
    distinct()
  
}
#!!!!!!
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
  filter(scenarioName %in% c("all_primary_CY_D.csv", "mostly_1L_CY_D.csv", "mostly_2L_CY_D.csv")) %>% 
  
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
names(scenario_flux[[13]])
x <- scenario_flux[[13]]
single <- x %>% filter(index == "mostly_2L_CY_D.csv 10")


#------ incorporate social cost of carbon adjust discount rate -------

#CONSIDER FOR A 4 % DISCOUNT RATE (CAN FILTER TO GET 2 OR 6)
socialDR_4 <- socialDR %>% filter(discount_rate == "4%") %>%
  rename(true_year = year) %>% 
  select(-X)

socialDR_2 <- socialDR %>% filter(discount_rate == "2%") %>%
  rename(true_year = year) %>% 
  select(-X) 


socialDR_6 <- socialDR %>% filter(discount_rate == "6%") %>%
  rename(true_year = year) %>% 
  select(-X) 


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

################
#############

#---------TEST ANDREW VS TOM SWNIFIELD WAY OF DOING THINGS: ------- 
#good their identical 
print(single) 

single <- single %>% left_join(socialDR_2, by = "true_year")

#select relevant columns only 
single <- single %>% select(index, true_year, scen_flux, SL_flux, scc_discounted, discount_rate)

#----- ANDREW's WAY -----------
#Andrew's way #-22135273
andrew <- single %>%  
  mutate(
    #get discounted annual impact for scenario and starting landscape
    SC_carbon_impact = (scen_flux * scc_discounted),
    SL_carbon_impact = (SL_flux * scc_discounted), 
    
    #sum over these values
    SC_summedC = sum(SC_carbon_impact,na.rm = TRUE), 
    SL_summedC = sum(SL_carbon_impact,na.rm = TRUE), 
    
    
    #minus SL from Scenario to get full impact 
    totalSCC = SC_summedC - SL_summedC
  )


#Swinfield's way #-22135273

Swinfield <- single %>%  
  mutate(
    #get discounted annual impact for scenario and starting landscape
    carbon_impact = (scen_flux - SL_flux) * scc_discounted, 
    totalSCC = sum(carbon_impact,na.rm = TRUE))

################
#############

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

##THIS BIT JUST TELLS US WHAT THE NET FLUX IS. 
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

x <- plot_grid(carbonImpact, sFlux, rel_heights = c(1,0.5))
plot_grid(x, scc, nrow = 2, rel_heights = c(1,0.5))


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
output <- all_discount_rates %>% select(index, production_target, scenarioName,scenarioStart,
                                        TOTcarbon_impact,TOTcarbon_impact_lwr, TOTcarbon_impact_upr, discount_rate) %>% 
  cbind(outcome = "carbon")
output <- unique(output)
saveRDS(output, "R_code/AllOutcomesFigure/Data/carbon.rds")



# final_carbon_df %>% 
#   ##!!
# #  final_carbon %>% 
#   #!!
#   
#   select(
#   index, scenarioName, 
#   production_target, 
#   TOTcarbon_impact,TOTcarbon_impact_lwr,TOTcarbon_impact_upr ) %>%  
#  # distinct() %>% 
#   
#  # filter(scenarioName == "primary_deforested_CY_D.csv") %>% 
#   ggplot(aes(x = production_target, y =TOTcarbon_impact))+
#   geom_point()+
#   #geom_line(aes(y = SL_ACD_change), colour ="black",lty="longdash", alpha=.2)+
#   facet_wrap(~scenarioName, ncol = 4)+#, scales = "free")+
#   theme_bw()+
#   theme(legend.position = "none")
#   

#------- Calculate scenario composition for plots---------

#get the amount of hab in each starting landscape 

habInStart <- all_start_landscape %>% select(scenarioStart) %>% unique() %>% 
  mutate(originalOG = c(1,0.2,0.2,0.2,0.2,0.8), 
         original1L = c(0,0.8,0.6,0,0,0), 
         original2L = c(0.2,0,0,0.8,0.6,0))


#build a function that calculate proportion of TOTAL and proportion of REMAINING habitat
#in each scenario 

prop_OG_fun <- function(x){
  
  #proportion of TOTAL landscape [1000 parcels] in different habitat type 
  x %>% group_by(index, production_target) %>% 
    #total OG
    mutate(propOG = sum(num_parcels[habitat == "primary"])/1000,
           propPlant = sum(num_parcels[habitat %in% c("eucalyptus_current", "albizia_current", "albizia_future","eucalyptus_future")])/1000,   
           #prop-1L in the scenario landscape
           prop1L = sum(num_parcels[habitat == "once-logged"])/1000,
           #proportion of 2-L in the scenario landscape
           prop2L = sum(num_parcels[habitat == "twice-logged"])/1000,
           propRestored = sum(num_parcels[habitat == "restored"])/1000) %>%  
    
    #get starting landscape
    mutate(scenarioStart = scenarioName) %>% 
    mutate(scenarioStart = str_remove(scenarioStart, "_IY_ND.csv")) %>%
    mutate(scenarioStart = str_remove(scenarioStart, "_CY_ND.csv")) %>%
    mutate(scenarioStart = str_remove(scenarioStart, "_IY_D.csv")) %>%
    mutate(scenarioStart = str_remove(scenarioStart, "_CY_D.csv")) %>% 
    ungroup %>% 
    
    #get total amount of each habitat in STARTING landscape for a scenario
    left_join(habInStart, by = "scenarioStart") %>% 
    
    #calculate PROPORTION of REMAINING original habitat type 
    #(nb there can actually be more once-logged or twice-logged forest in scenario than scenarioStart, if primary forest is logged)
    mutate(remainingOG = propOG/originalOG, 
           remaining1L = prop1L/original1L, 
           remaining2L = prop2L/original2L) %>%  
    #correct for INF values for if dividing by 0
    mutate_at(vars(remainingOG, remaining1L, remaining2L), ~ ifelse(is.infinite(.) | is.nan(.), 0, .)) %>%
    
    select(index, production_target, scenarioName,scenarioStart,
           propOG, propPlant,prop1L,prop2L,
           remainingOG,remaining1L,remaining2L,propRestored) %>% unique()
  
}

propOGcomp <- prop_OG_fun(scenario_composition) %>% ungroup


# -------combine carbon values and scenario composition ---------------------
addOGfun <- function(x){
  x %>% left_join(propOGcomp, by = c("index", "production_target", "scenarioName"))
}

final_carbon_df_4DR <- addOGfun(final_carbon_df_4DR)
final_carbon_df_2DR <- addOGfun(final_carbon_df_2DR)
final_carbon_df_6DR <- addOGfun(final_carbon_df_6DR)



#---- BIVARIATE PLOTTING PARAMETRES -------------
library(biscale)
COL <- "DkBlue2" # define colour pallete
COL <- "BlueOr"
#get colours for bivariate plotting
biv_pallete <- bi_pal(COL, dim =4 ) # for plotting
cols <- data.frame(bi_pal(COL, dim = 4, preview = FALSE))
colnames(cols) <- c("hex")
cols <- cols %>% mutate(bi_class = rownames(.))

textSize  <- 15

#make bivar legend
primary_legend <- bi_legend(pal = "BlueOr", dim = 4, 
                            xlab = "Proportion old-growth", 
                            ylab = "Proportion once-logged", size = textSize) 

onceL_legend <- bi_legend(pal = "BlueOr", dim = 4, 
                          xlab = "Proportion remaining old-growth", 
                          ylab = "Proportion remainng once-logged",size = textSize) 

twiceL_legend <- bi_legend(pal = "BlueOr", dim = 4, 
                           xlab = "Proportion remaining old-growth", 
                           ylab = "Proportion remainng twice-logged",size = textSize) 

all_legend <- plot_grid(primary_legend,onceL_legend,twiceL_legend, ncol =3)

#assign scenarios the colours from the bivariate plot for primary start
bivariate_colours_PRIM <- function(X){
  X %>%  bi_class(x = propOG, y = prop1L, dim = 4, style = "equal") %>%  
    left_join(cols, by = "bi_class") # add hex colours
}

final_carbon_df_4DR<- bivariate_colours_PRIM(final_carbon_df_4DR) %>% rename(hexP = hex)
final_carbon_df_2DR <- bivariate_colours_PRIM(final_carbon_df_2DR)  %>% rename(hexP = hex)
final_carbon_df_6DR <-bivariate_colours_PRIM(final_carbon_df_6DR)  %>% rename(hexP = hex)


#assign scenarios the colours from the bivariate plot for mostly 1L start
bivariate_colours_1L <- function(X){
  X %>%  bi_class(x = remainingOG, y = remaining1L, dim = 4, style = "equal") %>%  
    left_join(cols, by = "bi_class") # add hex colours
} 

final_carbon_df_4DR<- bivariate_colours_1L(final_carbon_df_4DR) %>% rename(hex1L = hex)
final_carbon_df_2DR <- bivariate_colours_1L(final_carbon_df_2DR) %>% rename(hex1L = hex)
final_carbon_df_6DR <-bivariate_colours_1L(final_carbon_df_6DR) %>% rename(hex1L = hex)


#assign scenarios the colours from the bivariate plot for mostly 2L start
bivariate_colours_2L <- function(X){
  X %>%  bi_class(x = remainingOG, y = remaining2L, dim = 4, style = "equal") %>%  
    left_join(cols, by = "bi_class") # add hex colours
}

final_carbon_df_4DR<- bivariate_colours_2L(final_carbon_df_4DR) %>% rename(hex2L = hex)
final_carbon_df_2DR <- bivariate_colours_2L(final_carbon_df_2DR) %>% rename(hex2L = hex)
final_carbon_df_6DR <-bivariate_colours_2L(final_carbon_df_6DR) %>% rename(hex2L = hex)

#hex shows colours for 1L vs primary.
#hex_2L shows colours for 2L vs primary

#------build plotting function -----

scenario_filters <- c("all_primary_CY_D.csv", "mostly_1L_CY_D.csv", "mostly_2L_CY_D.csv")

plot_fun <- function(x){
  
  x <- x %>%
    filter(scenarioName %in% scenario_filters) 
  
  
  #if scenario contains plantation add cross, if containts restored, add triangle  
  x <- x %>%
    mutate(is_cross = ifelse(propPlant > 0, "Cross", "Point"))
  
  #add conditional shapes
  # x <- x %>%
  #   mutate(is_cross = case_when(
  #     propPlant > 0 & propRestored > 0 ~ 3, #3 = cross
  #     propPlant > 0 & propRestored == 0 ~ 3, #cross
  #     propRestored > 0 & propPlant == 0 ~ 25,  # 6 = triangle 
  #     TRUE ~ 1)) #2 = point
  
  # #add conditional alpha 
  # x <- x %>%
  #   mutate(alpha = case_when(
  #     is_cross == Cross ~ 0.05,
  #     is_cross == 6 ~ 1,
  #     is_cross == 2 ~ 1 #cir
  #   ))
  
  # x <- x %>%
  #    mutate(alpha = case_when(
  #        is_cross == "Cross" ~ 0.05,TRUE ~ 1))
  # 
  x <- x %>% 
    mutate(
      size = case_when(
        is_cross == "Cross" ~ 4, TRUE ~ 1))  # Increase to size 4 for is_cross == 3 (crosses)
  #       TRUE ~ 1  # Set the default size for other points
  #     ))
  
  # #INSERT CASE_WHEN CONDITIONAL SIZE HERE
  # x <- x %>% 
  #   mutate(
  #     size = case_when(
  #       is_cross == 3 ~ 4,  # Increase to size 4 for is_cross == 3 (crosses)
  #       TRUE ~ 1  # Set the default size for other points
  #     ))
  
  #reorder facet order 
  
  x$scenarioName <- factor(x$scenarioName, levels = c(
    #select scenarios I want to plot 
    
    "all_primary_CY_D.csv", "mostly_1L_CY_D.csv",  "mostly_2L_CY_D.csv"))#,
  
  # "all_primary_CY_ND.csv","all_primary_IY_D.csv", "all_primary_IY_ND.csv",
  # "mostly_1L_CY_ND.csv", "mostly_1L_IY_D.csv", "mostly_1L_IY_ND.csv",
  # "mostly_1L_deforested_CY_D.csv", "mostly_1L_deforested_CY_ND.csv", "mostly_1L_deforested_IY_D.csv", "mostly_1L_deforested_IY_ND.csv", 
  #"mostly_2L_CY_ND.csv", "mostly_2L_IY_D.csv", "mostly_2L_IY_ND.csv",
  #"mostly_2L_deforested_CY_D.csv", "mostly_2L_deforested_CY_ND.csv", "mostly_2L_deforested_IY_D.csv","mostly_2L_deforested_IY_ND.csv",
  #"primary_deforested_CY_D.csv", "primary_deforested_CY_ND.csv", "primary_deforested_IY_D.csv","primary_deforested_IY_ND.csv"))
  
  
  x %>%  ggplot(aes(x = production_target, y = TOTcarbon_impact))+
    # conditionally colour so that if we plot bivariate between proportion of primary and proportion of least logged (either 1L or 2L depending on starting landscape) in the scenario    
    geom_point(aes(x = production_target,
                   y = TOTcarbon_impact/10000000000,
                   
                   colour = case_when(
                     scenarioStart %in% c("all_primary", "primary_deforested") ~ hexP,
                     scenarioStart %in% c("mostly_1L", "mostly_1L_deforested") ~ hex1L,
                     scenarioStart %in% c("mostly_2L", "mostly_2L_deforested") ~ hex2L),
                   
                   shape = as.factor(is_cross), # Add cross to plantation scenarios 
                   #  alpha = alpha,  # Calculate alpha values conditionally
                   alpha = ifelse(x$is_cross == "Cross", 0.1, 0.5),  # Calculate alpha values conditionally
                   
                   #  size = size
    ), 
    position = position_jitter(width = 0.05, height = 0.1)
    )+
    
    
    # position = position_jitter(width = 0.11, height = 0.01)
    
    
    
    scale_colour_identity()+
    scale_shape_manual(values = c("Point" = 19, "Cross" = 3)) + # Define shape mapping
    #remove tiny production targets 
    xlim(0.001, 1)+
    ylim(0, -2.5)+
    
    xlab("Production target")+
    ylab("Social cost of carbon (Billion USD$)")+
    
    #labs(colour = "Proportion of remaining old-growth forest spared")+
    # labs(colour = "Proportion of plantation in remaining landscape")+
    
    facet_wrap(~scenarioName, ncol = 4, scales = "free_y")+
    #   geom_hline(aes(yintercept = SL_geom_mean))+
    theme_bw(base_size = textSize)+
    theme(legend.position = "none")
}



#------ final plots ---------------
p4 <- plot_fun(final_carbon_df_4DR)
plot_grid(p4,all_legend, nrow = 2)
p2 <- plot_fun(final_carbon_df_2DR)
plot_grid(p2,all_legend, nrow = 2)
p6 <- plot_fun(final_carbon_df_6DR)
plot_grid(p6,all_legend, nrow = 2)

#----identify 2% discount rate high performing grey dot indexes -----

final_carbon_df_2DR %>% filter(hexP =="#d3d3d3"& scenarioStart == "all_primary") %>%  
  filter(TOTcarbon_impact/10000000000 >-0.5) %>% select(index) %>% unique

#ALL SEEM TO CONTAIN RESTORED FOREST
searchIndex <- "all_primary_CY_D.csv 349"
scenario_composition %>% filter(index ==searchIndex )
XX <- final_carbon_df_2DR %>% filter(index == searchIndex)

names(final_carbon_df_2DR)

#------remove any index with restored in composition (FOR NOW) -----

#this seems to get rid of the instances where plantation-containing and non-plantation 
#containing scenarioS overlap at low low discount rates 

filtered_df <- scenario_composition %>%
  filter(habitat == "restored")

# Extract the indices to be removed
indices_to_remove <- filtered_df$index 

# Remove rows with the selected indices
result_df <- final_carbon_df_2DR %>%
  filter(!index %in% indices_to_remove)
p2 <- plot_fun(result_df)
plot_grid(p2,all_legend, nrow = 2)



#---- unused code ----------

# 
# 

# #build plotting function

# plot_fun <- function(x){
#   
#   #reorder facet order 
#   
#   x$scenarioName <- factor(x$scenarioName, levels = c(
#     "all_primary_CY_D.csv","all_primary_CY_ND.csv","all_primary_IY_D.csv", "all_primary_IY_ND.csv",
#     "mostly_1L_CY_D.csv", "mostly_1L_CY_ND.csv", "mostly_1L_IY_D.csv", "mostly_1L_IY_ND.csv",
#     "mostly_1L_deforested_CY_D.csv", "mostly_1L_deforested_CY_ND.csv", "mostly_1L_deforested_IY_D.csv", "mostly_1L_deforested_IY_ND.csv", 
#     "mostly_2L_CY_D.csv", "mostly_2L_CY_ND.csv", "mostly_2L_IY_D.csv", "mostly_2L_IY_ND.csv",
#     "mostly_2L_deforested_CY_D.csv", "mostly_2L_deforested_CY_ND.csv", "mostly_2L_deforested_IY_D.csv","mostly_2L_deforested_IY_ND.csv",
#     "primary_deforested_CY_D.csv", "primary_deforested_CY_ND.csv", "primary_deforested_IY_D.csv","primary_deforested_IY_ND.csv"))
#   
#   #........................................................................................
#   #filter subset of scenrarios
#   filtered_scenarios <- c("all_primary_CY_D.csv", 
#                           "all_primary_CY_ND.csv", 
#                          "primary_deforested_CY_D.csv", "primary_deforested_CY_ND.csv",
#                          "mostly_1L_CY_D.csv", "mostly_1L_CY_ND.csv", 
#                            "mostly_1L_deforested_CY_D.csv", "mostly_1L_deforested_CY_ND.csv",
#                           "mostly_2L_CY_D.csv", "mostly_2L_CY_ND.csv", 
#                           "mostly_2L_deforested_CY_D.csv", "mostly_2L_deforested_CY_ND.csv")
#   
#   x <- x %>% filter(scenarioName %in%filtered_scenarios)
#   
#   x$scenarioName <- factor(x$scenarioName, levels = c(
#     "all_primary_CY_D.csv","all_primary_CY_ND.csv",
#     "primary_deforested_CY_D.csv", "primary_deforested_CY_ND.csv",
#     "mostly_1L_CY_D.csv", "mostly_1L_CY_ND.csv", 
#     "mostly_1L_deforested_CY_D.csv", "mostly_1L_deforested_CY_ND.csv", 
#     "mostly_2L_CY_D.csv", "mostly_2L_CY_ND.csv", 
#     "mostly_2L_deforested_CY_D.csv", "mostly_2L_deforested_CY_ND.csv"))
#   
#   
#   
#   
# #  ........................................................................................
#   x %>%  ggplot(aes(x = production_target, y = TOTcarbon_impact))+
#     #colour by proportion of remaining old-growth forest 
#    # geom_point(aes(x= production_target, y = TOTcarbon_impact, colour = propOriginalOG,alpha = 0.2), position = position_jitter(width = 0.05, height = 0.01)) + 
#     #colour by proportion of plantation in end scenarios
#    #geom_point(aes(x= production_target, y = TOTcarbon_impact, colour = propPlant,alpha = 0.2), position = position_jitter(width = 0.05, height = 0.01)) + 
#    # geom_point(aes(x= production_target, y = TOTcarbon_impact, colour = prop1L,alpha = 0.2), position = position_jitter(width = 0.05, height = 0.01)) + 
#     geom_point(aes(x= production_target, y = TOTcarbon_impact, colour = hex), position = position_jitter(width = 0.05, height = 0.01)) + 
#     scale_colour_identity()+
#     #scale_color_manual(values = color_scale) +
#     # xlab("Production target")+
#     # ylab("Proportion of landscape covered in megatrees (>50 m trees)")
#     #scale_color_viridis_d(option = "magma", alpha = 0.8, name = "Proportion of remaining old-growth forest spared") + # Use 'magma' color scale, change alpha if needed
# # 
# #      scale_colour_gradient(low = '#fee090', high = '#d73027',
# #                            breaks =c(0,1),
# #                            limits=c(0,1),#only show 0 and 4000 (parcels of legend
# #                            labels = c("0", "1"))+   #display as a percentage o concession
#     #remove tiny production targets 
#     xlim(0.001, 0.5)+
#     ylim(-10000000000, 40000000000)+
#     # labs(x = "Production target", 
#     #      y = "Carbon impact", 
#     #      colour = "Proportion of remaining old-growth forest spared")+
#     # labs(colour = "Proportion of plantation in remaining landscape")+
#     
#     facet_wrap(~scenarioName, ncol = 4)+#, scales = "free_y")+
#    # geom_hline(aes(yintercept = starting_l_prop_upr))+
#     theme_bw()+
#     theme(legend.position = "none")
#    }
# 
# p4 <- plot_fun(final_carbon_df_4DR)
# p2 <- plot_fun(final_carbon_df_2DR)
# p6 <- plot_fun(final_carbon_df_6DR)
# 

# #make legend

# #build plotting function
# legend_plot_fun <- function(x){
#  
#   
#   x %>%  ggplot(aes(x = production_target, y = TOTcarbon_impact))+
#     #colour by proportion of remaining old-growth forest 
#    # geom_point(aes(x= production_target, y = TOTcarbon_impact, colour = propOriginalOG), position = position_jitter(width = 0.05, height = 0.01)) + 
#   geom_point(aes(x= production_target, y = TOTcarbon_impact, colour = bi_class), position = position_jitter(width = 0.05, height = 0.01)) + 
#     # 
#     # scale_colour_gradient(low = '#fee090', high = '#d73027',
#     #                       breaks =c(0,1),
#     #                       limits=c(0,1),#only show 0 and 4000 (parcels of legend
#     #                       labels = c("0", "1"))+   #display as a percentage o concession
#     # 
#     labs(x = "Production target", 
#          y = "Carbon impact", 
#          colour = "Proportion of remaining old-growth forest spared")
# }
# legend_plot <-  final_carbon_df_4DR %>% filter(scenarioName == "all_primary_CY_D.csv") 
# legend_plot2 <- legend_plot_fun(legend_plot)
# #get legend
# legend <- get_legend(legend_plot2 + theme(legend.position = "bottom",         # c(0.5, 0.15),
#                                          legend.spacing.x = unit(1.0, 'cm'),
#                                          legend.title  = element_text(size  = 30, face = "bold"))) 
# #add legend function
# add_legend <-  function(x){
#   plot_grid(x, legend, 
#             nrow =2 , ncol = 1,
#             rel_heights = c(1, 0.1))
#  } 
# 
# #final figure 
# add_legend(p4)
# add_legend(p2)
# add_legend(p6)
# 
# #interactive plot 
# library(plotly)
# p <- ggplotly(x, tooltip = "index")
# 
# 
# x <- final_carbon_df_4DR %>% select(production_target,TOTcarbon_impact, hex, scenarioName)
# x %>% filter(scenarioName == "all_primary_CY_D.csv") %>% 
#   ggplot(aes(x = production_target, y = TOTcarbon_impact, colour = hex))+
#   geom_point()+
#   scale_colour_identity()
#   #colour by proportion of remaining old-growth forest 
#   # geom_point(aes(x= production_target, y = TOTcarbon_impact, colour = propOriginalOG), position = position_jitter(width = 0.05, height = 0.01)) + 
#   #geom_point(aes(x= production_target, y = TOTcarbon_impact, colour = bi_class), position = position_jitter(width = 0.05, height = 0.01))
# 
# #the lowest carbon impact is loads of albizia established solely on deforested land
# X <-final_carbon_df %>% filter(production_target > 0.01) %>% filter(TOTcarbon_impact == min(TOTcarbon_impact))# %>% select(index) %>% unique %>% slice(1)
# Z <- scenario_composition %>% filter(index =="primary_deforested_CY_D.csv 673" )
