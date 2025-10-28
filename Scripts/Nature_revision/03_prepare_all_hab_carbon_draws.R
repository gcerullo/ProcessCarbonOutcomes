#Combine posterior draws across habitat types 

library(brms)
library(tidyverse)
library(data.table)
library(broom)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(bayesplot)
library(tidybayes)
set.seed(123)

#Extract and organise ACD predictions for each of primary, once_logged, restored, twice-logged

#read in inputs ###

data_phil <- read.csv("RawData/Philipson20_PlotData2.csv")
Logged <- subset(data_phil, Forest == "Logged")
UnLogged <- subset(data_phil, Forest == "UnLogged")

# Ensure factors are as expected
Logged$FACE <- factor(Logged$FACE)
UnLogged$MeasureTime <- factor(UnLogged$MeasureTime)
unique(UnLogged$MeasureTime)


log_brm <- readRDS("Models/logged_restored_model.rds")
primary_brm <- readRDS("Models/primary_model.rds")
plantation <- readRDS("Outputs/plantation_carbon_draws.rds") %>%  
  rename(functionalhabAge = plantationAge, 
         habitat = species)

# ------------------------------------------------
# 1. Set up newdata for each model
# ------------------------------------------------
years <- 0:60
faces <- levels(Logged$FACE)

# Logged (once_logged + restored)
new_logged <- expand.grid(
  YearsSinceLogging = years,
  FACE = faces
)

# Primary forest (intercept-only model)
new_primary <- tibble(dummy = 1)  # just one row since it's intercept-only


# --------------------------------------------------------------------------
# 2. Extract 500 draws of expected ACD in primary, restored and once_logged
# --------------------------------------------------------------------------

# --- Logged habitats ---
logged_draws <- log_brm %>%
  add_epred_draws(
    newdata   = new_logged,
    re_formula = NA,     # exclude random effects
    ndraws     = 500     # number of posterior draws
  ) %>% ungroup %>% 
  mutate(
    habitat = recode(FACE,
                     "Baseline" = "once_logged",
                     "ProjectScenario" = "restored"),
    YearsSinceLogging = as.integer(YearsSinceLogging)
  ) %>%
  select(.draw, YearsSinceLogging, habitat, .epred)

#check correct num of draws
logged_draws %>%
  count(YearsSinceLogging, habitat) %>%
  summarise(
    min_draws = min(n),
    max_draws = max(n),
    total_combinations = n()
  )

# --- Primary forest ---
primary_draws <- primary_brm %>%
  add_epred_draws(
    newdata   = tibble(dummy = 1),  # one row since it's intercept-only
    re_formula = NA,
    ndraws     = 500
  ) %>%
  ungroup() %>%
  mutate(
    habitat = "primary",
    YearsSinceLogging = NA_integer_
  ) %>%
  select(.draw, YearsSinceLogging, habitat, .epred)

# ------------------------------------------------
#  Combine all draws
# ------------------------------------------------
ACD_draws <- bind_rows(logged_draws, primary_draws) %>% 
  rename(draw = .draw, 
                  functionalhabAge = YearsSinceLogging, 
                  ACD = .epred)
        
#also get slopes for once-logged forest recovery for each draw
delta <- 1  # 1-year slope
slopes_once_logged <- logged_draws %>%
  filter(habitat == "once_logged") %>%
  arrange(.draw, YearsSinceLogging) %>%
  group_by(.draw) %>%
  summarise(
    slope_once_logged = (.epred[YearsSinceLogging == delta] - .epred[YearsSinceLogging == 0]) / delta,
    .groups = "drop"
  ) %>%  
  rename(draw = .draw)

# ------------------------------------------------
# enforce plateau 
# ------------------------------------------------
# for each posterior draw, if the ACD estimate surpasses the primary estimate, replace with the primary estimate to enforce plateu
# Convert to data.table
dt <- as.data.table(ACD_draws)
primary_dt <- dt[habitat == "primary", .(draw, primary_ACD = ACD)]

# Merge primary draws on draw
dt <- merge(dt, primary_dt, by = "draw", all.x = TRUE)

# Apply plateau per draw
dt[habitat != "primary" & !is.na(functionalhabAge), 
   ACD := pmin(ACD, primary_ACD)]

# Drop helper column
dt[, primary_ACD := NULL]
ACD_draws <- as_tibble(dt)

# ------------------------------------------------
# Summary  plot
# ------------------------------------------------
ACD_summary <- ACD_draws %>%
  group_by(habitat, functionalhabAge) %>%
  summarise(
    mean = mean(ACD),
    lwr  = quantile(ACD, 0.025),
    upr  = quantile(ACD, 0.975),
    .groups = "drop"
  )

ggplot(ACD_summary, aes(x = functionalhabAge, y = mean,
                        color = habitat, fill = habitat)) +
  # ribbons for all habitats
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = habitat),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  
  # primary plateau line
  geom_hline(data = filter(ACD_summary, habitat == "primary"),
             aes(yintercept = mean),
             linetype = "dashed", color = "black") +
  
  # primary CI as rectangle spanning full x-axis
  geom_rect(data = filter(ACD_summary, habitat == "primary"),
            aes(xmin = -Inf, xmax = Inf,
                ymin = lwr, ymax = upr),
            fill = "grey80", alpha = 0.2, inherit.aes = FALSE) +

  theme_minimal(base_size = 14) +
  labs(x = "Years Since Logging",
       y = "Aboveground Carbon Density (ACD)",
       color = "Habitat type", fill = "Habitat type")


#INFER twice-logged dynamics ####


#DEFINE KEY PARAMS ####
#define amount of ACD loss from second rotation (ACD loss from getting 31.2m3 of timber in second harvest)

#RATIONALE - BEING USED for 2L -> 2L: assume forest that starts out as twice logged was logged fifteen years 
#after first logging (e.g. first harvest = t= -15, second harvest t = 0; tracks reality)

#1. For model of ACD recovery after once-logging,  it looks like the first harvesting rotation led to a decline of 148.8714 mg/ha   
ACD_l1_yr0 <-  ACD_summary %>% filter(functionalhabAge == 0 &habitat == "once_logged") %>% select(mean) %>% pull()
primary_ACD <- ACD_summary %>%  filter(habitat == "primary") %>% select(mean) %>% pull()
decline_from_first_logging <- primary_ACD - l1_yr0
#2 Thus 148.8714/112.96 =  1.317912 decline in ACD per m3 harvest
carbon_loss_pr_m3 <- decline_from_first_logging/112.96 
#3. Again, from models, once_logged forest has an ACD of 75.08815  after 15 years of recovery.  
ACD_l1_yr15 <- ACD_summary %>% filter(functionalhabAge == 15 &habitat == "once_logged") %>% select(mean) %>% pull()
ACD_l1_yr30 <- ACD_summary %>% filter(functionalhabAge == 30 &habitat == "once_logged") %>% select(mean) %>% pull()

#4. Twice-logging removes 31.24 m3 of timber (+/- 10.4) 
# If I assume the same ACD loss per m3 harvested as once_logged 
ACD_loss_from_2nd_harvest <-  31.24 * carbon_loss_pr_m3 #= 41.17158 = ACD loss from second harvest
#7. Thus given a 15yr once_logged ACD before logging, then:
ACD_2L_starting2L_15yrAfter1L <-  ACD_l1_yr15 - ACD_loss_from_2nd_harvest #= 40.62932 in twice-logged in yr 0   
ACD_2L_yr0_30yrAfter1L_30yrAfter1L <- ACD_l1_yr30 - ACD_loss_from_2nd_harvest #= 40.62932 in twice-logged in yr 0   


# --- constants ---
vol_first_rotation <- 112.96   # m³ removed in first harvest
vol_second_rotation <- 31.24   # m³ removed in second harvest

# --- extract per-draw ACD values ---
once_draws <- ACD_draws %>% filter(habitat == "once_logged")
primary_draws <- ACD_draws %>%
  filter(habitat == "primary") %>%
  group_by(draw) %>%
  summarise(primary_ACD = mean(ACD, na.rm = TRUE), .groups = "drop")

# --- get once-logged ACD at specific ages per draw ---
once_age0  <- once_draws %>% filter(functionalhabAge == 0)  %>%
  select(draw, once_ACD_age0  = ACD)
once_age15 <- once_draws %>% filter(functionalhabAge == 15) %>%
  select(draw, once_ACD_age15 = ACD)
once_age30 <- once_draws %>% filter(functionalhabAge == 30) %>%
  select(draw, once_ACD_age30 = ACD)

# --- per-draw ACD loss per m³ harvested during first rotation ---
loss_per_m3_df <- primary_draws %>%
  left_join(once_age0, by = "draw") %>%
  mutate(
    decline_first = primary_ACD - once_ACD_age0,
    loss_per_m3   = decline_first / vol_first_rotation
  )

# --- calculate per-draw ACD loss from second harvest ---
twice_loss <- loss_per_m3_df %>%
  mutate(ACD_loss_second = vol_second_rotation * loss_per_m3)

ACD_twice_draws <- twice_starting_ACD %>%
  left_join(slopes_once_logged, by = "draw") %>%
  tidyr::expand_grid(functionalhabAge = functional_habYears) %>%
  mutate(
    # shift slope from once-logged to new intercept
    ACD_twice_logged_15yrStart = ACD_2L_start_15yrAfter1L + slope_once_logged * functionalhabAge,
    ACD_twice_logged_30yrStart = ACD_2L_start_30yrAfter1L + slope_once_logged * functionalhabAge
  ) %>%
  select(draw, functionalhabAge, ACD_twice_logged_15yrStart, ACD_twice_logged_30yrStart)


# #PLOT 
# # Combine once-logged and twice-logged draws for plotting
# functional_habYears <- 0:60
# 
# # Expand twice-logged data if not already
# ACD_twice_draws_plot <- ACD_twice_draws %>%
#   pivot_longer(cols = starts_with("ACD_twice_logged"), 
#                names_to = "scenario", 
#                values_to = "ACD") %>%
#   mutate(
#     scenario = recode(scenario,
#                       "ACD_twice_logged_15yrStart" = "twice_logged_15yrStart",
#                       "ACD_twice_logged_30yrStart" = "twice_logged_30yrStart")
#   )
# 
# # Once-logged draws
# once_logged_plot <- ACD_draws %>%
#   filter(habitat == "once_logged") %>%
#   select(draw, functionalhabAge, ACD) %>%
#   mutate(scenario = "once_logged")
# 
# # Primary plateau draws (just mean line for clarity)
# primary_plot <- ACD_draws %>%
#   filter(habitat == "primary") %>%
#   group_by(functionalhabAge) %>%
#   summarise(ACD = mean(ACD), .groups = "drop") %>%
#   mutate(scenario = "primary")
# 
# # Combine all for plotting
# plot_df <- bind_rows(once_logged_plot, ACD_twice_draws_plot, primary_plot)
# 
# # Plot
# ggplot(plot_df, aes(x = functionalhabAge, y = ACD, group = interaction(draw, scenario), color = scenario)) +
#   geom_line(alpha = 0.15, size = 0.8) +  # posterior draws lightly
#   # Overlay mean lines for clarity
#   stat_summary(aes(group = scenario), fun = mean, geom = "line", size = 1.2, color = "black") +
#   scale_color_manual(values = c("once_logged" = "blue", 
#                                 "twice_logged_15yrStart" = "red", 
#                                 "twice_logged_30yrStart" = "orange",
#                                 "primary" = "darkgreen")) +
#   theme_minimal(base_size = 14) +
#   labs(
#     x = "Years Since Logging",
#     y = "Aboveground Carbon Density (ACD)",
#     color = "Scenario",
#     title = "ACD Recovery Curves: Once-logged, Twice-logged, and Primary Forest"
#   )
#Briefly: 
#1 = P -> 2L
#2 = 2L -> 2L 
#3 = 1L -> 2L 


#NOTE THERE ARE THREE TYPES OF TWICE-LOGGING CURVES THRU TIME BUILT HERE: 
#1. ACD_2L_yr0_30yrAfter1L_30yrAfter1L - second rotations happens 30 years after 1st logging. Always the case in all_primary and all_primary_deforested scenarios

#2, ACD_2L_starting2L_15yrAfter1L - 
#here the second rotation happens 15 years after 1st rotation 
#( the second rotation is assumed to have happened in t-1, before our scenarios start) 
# We incoroprate this for all scenarios that start 2L - as this more rapid re-harvest of forest more closely resembles the true 
#logging dynamics in the landscape where we surveyed twice logged forests.


#i.e in   t-16----t-15 --------------------------t-1--------t-0------------------------------------t-60
#         P        1L                             2L        2L                                      2L

#We also make this true for scenarios that go from 1L (at the beginning of the scenarios) -> 2L....
#....which is operationalised by not allowing 1L-2L transitions for 15 years....see below

#3. twice_logged_delays  - second rotation happens at varying intervals after the once-logging from 15-30 years, depending on the harvesting delay. 
# this is important for scenarios that begin with once_logged forest (eg. mostly_1L and mostly_1L_deforested scenarios) 
# the delay matters because if the second harvest happens after a big time delay (e.g. 25 yrs) then ACD was recovering all this time in once_logged, so the 2nd harvest leaves behind higher ACD 
#NB - no logging of 1L -> 2L is allowed in the first 15 yrs, as the forest was logged in t-1. 


