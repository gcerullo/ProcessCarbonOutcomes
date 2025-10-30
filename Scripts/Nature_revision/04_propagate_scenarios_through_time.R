#propagate the uncertainty and estimates of scenario ACD 
#NOTE - need to go back and correct belowground process.
library(dplyr)
library(purrr)


#----------------read in scenarios -------------------------------

#yield matched scenarios where 1/30th of plantation conversion happens annually - with no time delay
#temporarily read in to get composition
scenarios_rm <- readRDS("Inputs/MasterAllScenarios.rds")
scenario_composition <- rbindlist(scenarios_rm, use.names=TRUE)

#yield matched scenarios where 1/30th of plantation conversion happens annually - WITH TIME DELAY 
scenarios <- readRDS("Inputs/MasterAllScenarios_withHarvestDelays.rds")

#check that we have the correct number of harvest delays per transition!!!
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

#___________________________________________________________________
#A) 1L- 1L modification. Forest that starts a scenario as once-logged was first harvest at yr t-15

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
# B)2L - 2L modification. Forest that starts with twice-logged was 
# Takes a list of data frames and applies two operations to each:
# (1) If both `original_habitat` and `functional_habitat` are "twice-logged",
#     change `functional_habitat` to "twice_logged_start".
# (2) If `original_habitat` is "once_logged", remove rows where `harvest_delay < 15`.
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

scenarios <- scenarios %>%  update_1L_functional_habitat %>% update_2L_functional_habitat()
#convret to datatable
scenarios <- lapply(scenarios, as.data.table)

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
all_draws <- primary %>% rbind(restored) %>% rbind(once_logged) %>% rbind(once_logged_start) %>%  
   rbind(plant)

#!!!!!select which twice-logged recovery trajectory to consider using 'slope factor' - 
#1 = same slope as 1L, 0.8 is 20% slower recovery, 1.2 = 20% faster recovery
diff_2Ls <- twice_logged %>% rbind(twice_lgged_start) %>% 
  filter(slope_factor == 1) %>% 
  select( draw, functionalhabAge, ACD, habitat) %>% unique()

all_draws<- all_draws %>% rbind(diff_2Ls)

#turn draws into a list where each draw is it's own df
all_draws <- all_draws %>%
  rename(functional_habitat = habitat) %>%  # rename first
  group_split(draw)     

#convret to datatable
all_draws <- lapply(all_draws, as.data.table)

#__________________________________________
#combine scenarios with carbon values
#__________________________________________

scenarios_test <- scenarios[[1]]

draw <- all_draws[[1]]

x <- scenarios_test

library(data.table)

add_carbon_fun <- function(x, draw) {
  # Ensure both are data.tables
  # setDT(x)
  # setDT(draw)
  
  # Separate primary and non-primary to handle joins cleanly
  x_primary <- x[functional_habitat == "primary"]
  x_other   <- x[functional_habitat != "primary"]
  
  # Join primary habitats
  joined_primary <- x_primary[
    draw[functional_habitat == "primary"],
    on = .(functional_habitat),
    nomatch = 0
  ]
  
  #remove addtional column created
  joined_primary[, i.functionalhabAge := NULL]
  
  # Join non-primary habitats
  joined_other <- x_other[
    draw[functional_habitat != "primary"],
    on = .(functional_habitat, functionalhabAge),
    nomatch = 0
  ]
  
  # Combine results
  result <- rbindlist(list(joined_primary, joined_other), use.names = TRUE)
  
  return(result)
}

scenarios <- lapply(scenarios_test, add_carbon_fun)


v <- result %>% filter(index == "all_primary_CY_D.csv 404")
#we have the incorrect number of harvest delays per transition!
x <-v %>% group_by(harvest_delay, functional_habitat) %>% count()


