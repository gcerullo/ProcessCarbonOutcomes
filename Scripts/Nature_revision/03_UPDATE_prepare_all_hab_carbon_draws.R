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
# for each posterior draw, if the ACD estimate surpasses the primary estimate, replace with the primary estimate to enforce plateau
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

# --- define alternative slopes as factors of once-logged slope ---
slope_factors <- c(0.5, 1, 1.5)  # slower = 50%, same = 100%, faster = 150%

# --- expand twice-logged draws across slope scenarios ---
ACD_twice_draws <- twice_starting_ACD %>%
  left_join(slopes_once_logged, by = "draw") %>%
  tidyr::expand_grid(
    functionalhabAge = functional_habYears,
    slope_factor = slope_factors
  ) %>%
  # apply scaled slope to each draw
  mutate(
    slope_scaled = slope_once_logged * slope_factor,
    # raw ACD trajectories for different slopes
    ACD_twice_logged_15yrStart_raw = ACD_2L_start_15yrAfter1L + slope_scaled * functionalhabAge,
    ACD_twice_logged_30yrStart_raw = ACD_2L_start_30yrAfter1L + slope_scaled * functionalhabAge
  ) %>%
  select(draw, functionalhabAge, slope_factor,
         ACD_twice_logged_15yrStart_raw, ACD_twice_logged_30yrStart_raw)

dt <- as.data.table(ACD_twice_draws)
dt <- merge(dt, primary_draws, by = "draw", all.x = TRUE)
dt[, ACD_twice_logged_15yrStart := pmin(ACD_twice_logged_15yrStart_raw, primary_ACD)]
dt[, ACD_twice_logged_30yrStart := pmin(ACD_twice_logged_30yrStart_raw, primary_ACD)]
dt[, primary_ACD := NULL]
ACD_twice_draws <- as_tibble(dt)

#PLOT

library(ggplot2)
library(dplyr)

# --- reshape for plotting ---
ACD_twice_long <- ACD_twice_draws %>%
  pivot_longer(
    cols = c(ACD_twice_logged_15yrStart, ACD_twice_logged_30yrStart),
    names_to = "start_age",
    values_to = "ACD"
  ) %>%
  mutate(
    start_age = case_when(
      start_age == "ACD_twice_logged_15yrStart" ~ "15yrAfter1L",
      start_age == "ACD_twice_logged_30yrStart" ~ "30yrAfter1L"
    ),
    slope_factor = factor(slope_factor, levels = slope_factors)
  )

# --- summarize posterior draws ---
ACD_twice_summary <- ACD_twice_long %>%
  group_by(functionalhabAge, slope_factor, start_age) %>%
  summarise(
    mean_ACD = mean(ACD),
    lwr_ACD  = quantile(ACD, 0.025),
    upr_ACD  = quantile(ACD, 0.975),
    .groups = "drop"
  )

# --- plot ---
ggplot(ACD_twice_summary, aes(x = functionalhabAge, y = mean_ACD,
                              color = slope_factor, fill = slope_factor)) +
  geom_ribbon(aes(ymin = lwr_ACD, ymax = upr_ACD), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~start_age, ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Years Since 2nd Logging",
    y = "Aboveground Carbon Density (ACD)",
    color = "Slope factor",
    fill = "Slope factor",
    title = "Twice-logged recovery trajectories with different slopes"
  )


