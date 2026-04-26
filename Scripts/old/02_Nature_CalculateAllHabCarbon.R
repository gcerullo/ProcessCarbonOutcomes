#09.06.24
#This code calculates carbon (above and belowground) of each habitat type through time 

#Code notes:
#Plantation ACD is calculate is from SSB inventory data 
#Primary, 1L and R is from Philipson and uses edited code originally written by 
#Philipson 2020; Science: https://github.com/PhilipsonChristopher/CarbonRecovery/blob/master/Fig1/Figure1Code.R


#(I adjust Philipson's code to allow estimates for 60 years (they
#estimate out to 30/35 years; I assume that slope and intercept of their models stay the same 
#and that 1L and R values plateu once they reach primary forest values)

#I estimate ACD in twice-logged forest, and use three 2L typologies 

#This code also incorporates belowground carbon and necromass processes

#Code output
##NB this code outputs a master csv called "allHabCarbon_60yrACD.csv"
# where ACD refers only to the raw above ground carbon density values and 
#where full_carbon incorporate belowground processes.
# Bayesian version of your carbon pipeline using brms
# 2025-10-02 (adapted from user's original script)
# NOTE: requires brms, tidyverse, data.table, cowplot, ggplot2

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

# ---------- PARAMETERS (PROBS NEED TO CHANGE) ----------
#DEFINE KEY PARAMS ####
#define amount of ACD loss from second rotation (ACD loss from getting 31.2m3 of timber in second harvest)

#RATIONALE - BEING USED for 2L -> 2L: assume forest that starts out as twice logged was logged fifteen years 
#after first logging (e.g. first harvest = t= -15, second harvest t = 0; tracks reality)

#1. For Chris Philipson's model of ACD recovery after once-logging,  it looks like the first harvesting rotation led to a decline of 155.9869 ACD_ha-1  (1 yrs after once-logged = 44.01312 +/- 31.20968), down  from ~200 ACD_ha-1 for primary forest, and removing 112.96 m3 (+/- 22.42) of timber in first rotation  
#2 Thus 155.9869/112.96 =  1.380904 decline in ACD per m3 harvest
#3. Again, from Philipson's models, once-logged forest has an ACD of 84.56991 (+/- 30.47371) after 15 years of recovery.  
#4. Twice-logging removes 31.24 m3 of timber (+/- 10.4) 
#5. I assume the same ACD loss per m3 harvested as once-logged 
#6. Thus 31.24 * 1.380904 = 43.13944 = ACD loss from second harvest
#7. Thus given a 15yr once-logged ACD before logging, then:
#      84.56991 - 43.94059 = 40.62932 in twice-logged in yr 0   
harvest2ndACD_loss <- 43.13944
ACD_2L_yr0_30yrAfter1L_30yrAfter1L <- 84.88416
ACD_2L_starting2L_15yrAfter1L <- 40.62932

# ---------- READ INPUTS ----------
hab_by_year <- read.csv("Inputs/HabByYears.csv", strip.white = TRUE) %>%
  rename(true_year = year,
         functionalhabAge = functional_habAge,
         habitat = transition_habitat) %>%
  select(-X) %>%
  filter(!str_detect(habitat, "improved"))

plantation <- read.csv("Outputs/plantation_carbon.csv") %>%
  select(-X) %>%
  mutate(habitat = case_when(habitat == "albizia" ~ "albizia_current",
                             habitat == "eucalyptus" ~ "eucalyptus_current"))

data_phil <- read.csv("RawData/Philipson20_PlotData2.csv")

Logged <- subset(data_phil, Forest == "Logged")
UnLogged <- subset(data_phil, Forest == "UnLogged")

# Ensure factors are as expected
Logged$FACE <- factor(Logged$FACE)
UnLogged$MeasureTime <- factor(UnLogged$MeasureTime)
unique(UnLogged$MeasureTime)
# ---------- FIT BRMS MODELS ----------
# 1) Main model: ACD ~ YearsSinceLogging * FACE + (1|Plot) + (1|LoggingMethod:Coupe)
# Use weakly informative priors
priors_main <- c(
  prior(normal(0, 50), class = "Intercept"),
  prior(normal(0, 10), class = "b"),
  prior(student_t(3, 0, 10), class = "sd") # group-level sd
)

# Fit model (adjust iter/chains/cores per your compute)
log_brm <- brm(
  formula = bf(ACD ~ YearsSinceLogging * FACE + (1 | Plot) + (1 | LoggingMethod:Coupe)),
  data = Logged,
  family = gaussian(),
  prior = priors_main,
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95)
)

# 2) Primary forest (UnLogged): intercept model with group-level MeasureTime
# priors_primary <- c(
#   prior(normal(0, 50), class = "Intercept")
# )
# 
# primary_brm <- brm(
#   formula = bf(ACD ~ 1),
#   data = UnLogged,
#   family = gaussian(),
#   prior = priors_primary,
#   chains = 4,
#   iter = 4000,
#   warmup = 1000,
#   cores = 4,
#   seed = 123,
#   control = list(adapt_delta = 0.95)
# )

# Priors
priors_primary <- c(
  prior(normal(0, 50), class = "Intercept"),  # wide prior for intercept
  prior(normal(0, 10), class = "sd"),         # prior for group-level SD
  prior(student_t(3, 0, 10), class = "sigma") # prior for residual SD
)

# brms model
primary_brm <- brm(
  formula = bf(ACD ~ 1 + (1 | MeasureTime)),
  data = UnLogged,
  family = gaussian(),
  prior = priors_primary,
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95)
)

#save BRMS models 
saveRDS(log_brm, "Models/logged_restored_model.rds")
saveRDS(primary_brm, "Models/primary_model.rds")

# ---------- Posterior Predictive Checks ----------

#Can start here:
log_brm <- readRDS("Models/logged_restored_model.rds")
primary_brm <- readRDS("Models/primary_model.rds")

# For Logged model
pp_check(log_brm)  # default density overlay
pp_check(log_brm, type = "hist")   # histograms of observed vs simulated
pp_check(log_brm, type = "scatter_avg")  # scatter of observed vs avg predicted
pp_check(log_brm, type = "stat", stat = "mean")  # distribution of predicted means
pp_check(log_brm, type = "stat", stat = "sd")    # distribution of predicted SDs

# For UnLogged (primary forest) model
pp_check(primary_brm)
pp_check(primary_brm, type = "hist")
pp_check(primary_brm, type = "stat", stat = "mean")
pp_check(primary_brm, type = "stat", stat = "sd")

hist(UnLogged$ACD)
mean(UnLogged$ACD)
mean(Logged$ACD)


# ------------------------------------------------
# 1. Set up newdata for each model
# ------------------------------------------------
years <- 0:60
faces <- levels(Logged$FACE)

# Logged (once-logged + restored)
new_logged <- expand.grid(
  YearsSinceLogging = years,
  FACE = faces
)

# Primary forest (intercept-only model)
new_primary <- tibble(dummy = 1)  # just one row since it's intercept-only


# --------------------------------------------------------------------------
# 2. Extract 500 draws of expected ACD in primary, restored and once-logged
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
                     "Baseline" = "once-logged",
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
# 3. Combine all draws
# ------------------------------------------------
ACD_draws <- bind_rows(logged_draws, primary_draws)

# ------------------------------------------------
# 4. Optional: enforce plateau
# ------------------------------------------------
# Get primary mean and CI from draws
primary_stats <- primary_draws %>%
  summarise(
    mean = mean(.epred),
    lwr  = quantile(.epred, 0.025),
    upr  = quantile(.epred, 0.975)
  )

ACD_draws <- ACD_draws %>%
  mutate(.epred = pmin(.epred, primary_stats$upr))  # enforce plateau ceiling

# ------------------------------------------------
# 5. Summarise if needed for plotting
# ------------------------------------------------
ACD_summary <- ACD_draws %>%
  group_by(habitat, YearsSinceLogging) %>%
  summarise(
    mean = mean(.epred),
    lwr  = quantile(.epred, 0.025),
    upr  = quantile(.epred, 0.975),
    .groups = "drop"
  )

# ------------------------------------------------
# 6. Plot
# ------------------------------------------------
ggplot(ACD_summary, aes(x = YearsSinceLogging, y = mean,
                        color = habitat, fill = habitat)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = habitat),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = primary_stats$mean, linetype = "dashed") +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = primary_stats$lwr, ymax = primary_stats$upr),
            fill = "grey50", alpha = 0.1, inherit.aes = FALSE) +
  theme_minimal(base_size = 14) +
  labs(x = "Years Since Logging",
       y = "Aboveground Carbon Density (ACD)",
       color = "Habitat type", fill = "Habitat type")







# 
# # ---------- POSTERIOR PREDICTIONS FOR 0:60 YEARS ----------
# years <- 0:60
# # Build newdata for both FACE levels found in Logged dataset
# faces <- levels(Logged$FACE)
# NEWDAT <- expand.grid(YearsSinceLogging = years, FACE = faces)
# 
# # Population-level (fixed effect) predictions: re_formula = NA excludes group-level effects so we get fitted marginal curves
# # Use posterior_epred to get draws on response scale (including Gaussian link)
# pp <- posterior_epred(log_brm, newdata = NEWDAT, re_formula = NA) # matrix draws x rows(newdata)
# 
# # Compute posterior summaries per NEWDAT row
# pp_mean <- apply(pp, 2, mean)
# pp_lwr <- apply(pp, 2, quantile, probs = 0.025)
# pp_upr <- apply(pp, 2, quantile, probs = 0.975)
# 
# 
# # Map FACE labels to your habitat names (to match earlier code)
# # original code later relabels "A: Natural regeneration" -> once-logged etc.
# # Here I'm assuming Logged$FACE has "Baseline" and "ProjectScenario" as before.
# NEWDAT$FACE <- factor(NEWDAT$FACE)
# 
# 


# 
# # ---------- BUILD L1_R dataframe (once-logged & restored) ----------
# # ---------- PREDICTIONS FOR LOGGED HABITATS ----------
# years <- 0:60
# faces <- levels(Logged$FACE)
# 
# NEWDAT <- expand.grid(YearsSinceLogging = years, FACE = faces)
# 
# # Posterior expected values (on response scale, population-level only)
# pp <- posterior_epred(log_brm, newdata = NEWDAT, re_formula = NA)
# 
# # Summarise draws into mean + 95% CI
# pred_summary <- t(apply(pp, 2, function(x) {
#   c(mean = mean(x),
#     lwr  = quantile(x, 0.025),
#     upr  = quantile(x, 0.975))
# })) %>% as.data.frame()
# 
# NEWDAT_pred <- bind_cols(NEWDAT, pred_summary) %>% 
#   rename(
#     lwr = `lwr.2.5%`,
#     upr = `upr.97.5%`
#   )
# 
# # ---------- PRIMARY FOREST PLATEAU ----------
# 
# # Determine the max primary ACD to impose plateau 
# # Get posterior draws for prime intercept (population-level)
# primary_epred <- posterior_epred(primary_brm, re_formula = NA) # draws x 1
# primary_mean  <- mean(primary_epred)
# primary_median <- median(primary_epred)
# primary_ci    <- quantile(primary_epred, probs = c(0.025, 0.975))
# 
# # ---------- PLOT ----------
# ggplot(NEWDAT_pred, aes(x = YearsSinceLogging, y = mean,
#                         color = FACE, fill = FACE)) +
#   # ribbons for logged scenarios
#   geom_ribbon(aes(ymin = lwr, ymax = upr, fill = FACE),
#               alpha = 0.2, color = NA) +
#   geom_line(linewidth = 1.2) +
#   
#   # primary plateau: horizontal ribbon + line
#   geom_hline(yintercept = primary_mean, color = "black", linetype = "dashed") +
#   geom_rect(aes(xmin = -Inf, xmax = Inf,
#                 ymin = primary_ci[1], ymax = primary_ci[2]),
#             fill = "grey40", alpha = 0.15, inherit.aes = FALSE) +
#   
#   theme_minimal(base_size = 14) +
#   labs(x = "Years Since Logging",
#        y = "Aboveground Carbon Density (ACD)",
#        color = "Management Type",
#        fill  = "Management Type") +
#   scale_color_brewer(palette = "Dark2") +
#   scale_fill_brewer(palette = "Dark2")
# 
# 
# ##### ENFORCE PLATEU #####
# # Enforce plateau in 1L and R: replace predictions > max_val with max_val (and maintain CI similarly)
# primary_mean
# primary_ci
# 
# # Enforce primary plateau on predictions
# NEWDAT_pred <- NEWDAT_pred %>%
#   mutate(
#     mean = pmin(mean, primary_mean),
#     lwr  = pmin(lwr,  primary_ci[1]),
#     upr  = pmin(upr,  primary_ci[2])
#   ) %>% 
#   mutate(
#     habitat = case_when(
#       FACE == "ProjectScenario" ~ "restored",
#       FACE == "Baseline" ~ "once-logged"
#     )
#   ) %>% 
#   rename(functional_habitat = habitat) %>%
#   mutate(functionalhabAge = as.integer(YearsSinceLogging)) %>% 
#   select(-FACE, YearsSinceLogging)
# 
# ggplot(NEWDAT_pred, aes(x = functionalhabAge, y = mean,
#                         color = functional_habitat , fill = functional_habitat)) +
#   # primary plateau: horizontal ribbon + line
#   geom_hline(yintercept = primary_mean, color = "black", linetype = "dashed") +
#   geom_rect(aes(xmin = -Inf, xmax = Inf,
#                 ymin = primary_ci[1], ymax = primary_ci[2]),
#             fill = "grey79", alpha = 0.1, inherit.aes = FALSE) +
#   
#   # ribbons for logged scenarios
#   geom_ribbon(aes(ymin = lwr, ymax = upr, fill = functional_habitat),
#               alpha = 0.2, color = NA) +
#   geom_line(linewidth = 1.2) +
#   
#  
#   theme_minimal(base_size = 14) +
#   labs(x = "Years Since Logging",
#        y = "Aboveground Carbon Density (ACD)",
#        color = "Management Type",
#        fill  = "Management Type") +
#   scale_color_brewer(palette = "Dark2") +
#   scale_fill_brewer(palette = "Dark2")
# 
# #------------------------------------
# #organise P - 1L and 1L - 1L transitions 
# #------------------------------------
# 
# # extract once-logged first 30 years 
# oncelogged <- NEWDAT_pred %>% filter(functional_habitat == "once-logged", functionalhabAge < 31)
# 
# P_1L_df <- oncelogged %>% mutate(original_habitat = "primary", habitat = "once-logged")
# L1_1L_df <- oncelogged %>% mutate(original_habitat = "once-logged", habitat = "once-logged")
# L1_1L_R_df <- L1_1L_df %>% mutate(habitat = "restored")
# 
# # extract RESTORED first 30 years 
# 
# #stays once-logged during delay years of restoration
# Restored <- NEWDAT_pred %>% filter(functional_habitat == "restored")
# 
# #organise P - R and 1L-R transitions 
# 
# P_R_df <- Restored %>% mutate(original_habitat = "primary", habitat = "restored")
# L1_R_df <- Restored %>% mutate(original_habitat = "once-logged", habitat = "restored")
# 
# # ---------- PLANTATION (60 yr) ----------
# plantation_df <- plantation %>%
#   rename(functionalhabAge = plantationAge) %>%
#   rename(functional_habitat = habitat)
# 
# # add zero-age rows for plantations:
# zero_eucPlant <- data.frame(ACD = 0, lwr_ACD = 0, upr_ACD = 0, habitat = "eucalyptus_current", functionalhabAge = 0)
# zero_albPlant <- data.frame(ACD = 0, lwr_ACD = 0, upr_ACD = 0, habitat = "albizia_current", functionalhabAge = 0)
# plantation_df2 <- bind_rows(zero_eucPlant, zero_albPlant, plantation) %>% rename(functional_habitat = habitat)
# 
# # expand by possible transitions
# habcrossAlb <- hab_by_year %>% filter(functional_habitat == "albizia_current") %>% select(original_habitat, habitat) %>% unique()
# habcrossPlantEuc <- hab_by_year %>% filter(functional_habitat == "eucalyptus_current") %>% select(original_habitat, habitat) %>% unique()
# alb_df <- plantation_df2 %>% filter(functional_habitat == "albizia_current") %>% expand_grid(habcrossAlb)
# euc_df <- plantation_df2 %>% filter(functional_habitat == "eucalyptus_current") %>% expand_grid(habcrossPlantEuc)
# plant_df <- bind_rows(alb_df, euc_df)
# 
# # ---------- PRIMARY (60 yr) ----------
# primary_vals <- data.frame(ACD = round(prime_mean, 0),
#                            lwr_ACD = PrimeCIs[1],
#                            upr_ACD = PrimeCIs[2],
#                            habitat = "primary",
#                            functionalhabAge = 0)
# 
# primary_cross <- hab_by_year %>% filter(functional_habitat == "primary") %>% select(original_habitat, habitat) %>% unique()
# primary_df <- expand_grid(primary_vals, primary_cross)
# 
# # ---------- DEFORESTED ----------
# deforested <- data.frame(functionalhabAge = 0, ACD = 0, lwr_ACD = 0, upr_ACD = 0, functional_habitat = "deforested")
# deforested_cross <- hab_by_year %>% filter(functional_habitat == "deforested") %>% select(original_habitat, habitat) %>% unique()
# deforested_df <- expand_grid(deforested, deforested_cross)
# 
# # ---------- COMBINE (except twice-logged yet) ----------
# carbhabs <- bind_rows(deforested_df, primary_df, plant_df, P_R_df, L1_R_df, L1_1L_R_df, P_1L_df, L1_1L_df)
# 
# # ---------- TWICE-LOGGED: Bayesian slope extraction and predictions ----------
# # In original code you used slope from a simple lm on the once-logged subset.
# # Here we extract posterior slope for 'YearsSinceLogging' and the interaction from the brms fit.
# 
# # Identify parameter names for slopes:
# # population-level slope for YearsSinceLogging is "b_YearsSinceLogging"
# # interaction slope (YearsSinceLogging:FACEProjectScenario) will be "b_YearsSinceLogging:FACEProjectScenario"
# post <- as_draws_df(m1_brm) # requires posterior package; brms returns draws as a draws_df
# 
# # Extract relevant parameters: if FACE has specific naming, check draws columns
# # safe names:
# b_years <- post %>% select(starts_with("b_YearsSinceLogging")) %>% names()
# # We expect at least "b_YearsSinceLogging" and possibly "b_YearsSinceLogging:FACEProjectScenario"
# # Grab the main and interaction if present
# main_name <- "b_YearsSinceLogging"
# int_names <- names(post)[str_detect(names(post), "b_YearsSinceLogging:")]
# 
# # Create slope draws for baseline (Baseline) and for ProjectScenario (if present)
# if (length(int_names) == 0) {
#   # no interaction present (unlikely) -> slope is main only
#   slope_baseline_draws <- post[[main_name]]
#   slope_project_draws <- post[[main_name]]
# } else {
#   # when an interaction exists, slope for ProjectScenario = main + interaction
#   slope_baseline_draws <- post[[main_name]]
#   slope_project_draws <- post[[main_name]] + post[[int_names[1]]]
# }
# 
# # Assign which FACE corresponds to once-logged in your mapping:
# # earlier mapping: "A: Natural regeneration" -> once-logged
# # We'll treat baseline slope as once-logged slope if Baseline corresponds to "A: Natural regeneration".
# # To be explicit: check what factor reference is in the model:
# print(m1_brm)
# 
# # We'll compute posterior mean slope for once-logged (assuming "Baseline" corresponds to once-logged).
# slope_once_mean <- mean(slope_baseline_draws)
# slope_once_lwr <- quantile(slope_baseline_draws, probs = 0.025)
# slope_once_upr <- quantile(slope_baseline_draws, probs = 0.975)
# 
# # Build twice-logged predicted vectors using the same approach you used:
# predicted2L_values <- ACD_2L_yr0_30yrAfter1L_30yrAfter1L + slope_once_mean * years
# predicted2L_values_starting <- ACD_2L_starting2L_15yrAfter1L + slope_once_mean * years
# 
# # Build a data.frame maker that also propagates slope uncertainty into CI (approx):
# # For each year we can compute distribution: ACD_intercept + slope_draws * year
# twice_pred_draws_30 <- sapply(years, function(y) ACD_2L_yr0_30yrAfter1L_30yrAfter1L + slope_baseline_draws * y)
# twice_pred_draws_starting <- sapply(years, function(y) ACD_2L_starting2L_15yrAfter1L + slope_baseline_draws * y)
# 
# # Summaries:
# twice_df_30 <- data.frame(
#   ACD = apply(twice_pred_draws_30, 2, mean),
#   lwr_ACD = apply(twice_pred_draws_30, 2, quantile, .025),
#   upr_ACD = apply(twice_pred_draws_30, 2, quantile, .975),
#   true_age = years
# ) %>%
#   mutate(habitat = "twice-logged", original_habitat = "primary") %>%
#   rename(functionalhabAge = true_age) %>%
#   mutate(functionalhabAge = functionalhabAge + 1) # to match your +1 shift
# 
# twice_df_starting <- data.frame(
#   ACD = apply(twice_pred_draws_starting, 2, mean),
#   lwr_ACD = apply(twice_pred_draws_starting, 2, quantile, .025),
#   upr_ACD = apply(twice_pred_draws_starting, 2, quantile, .975),
#   true_age = years
# ) %>%
#   mutate(habitat = "twice-logged", original_habitat = "twice-logged",
#          functional_habitat = "twice-logged") %>%
#   rename(functionalhabAge = true_age)
# 
# # Now adapt the P -> 2L construction (one-step once-logged then 2L at t+30)
# yr1_30_primary_2L <- oncelogged %>% rename(true_year = time_since_intervention) %>%
#   filter(true_year >= 1, true_year <= 29) %>%
#   mutate(habitat = "twice-logged",
#          original_habitat = "primary",
#          functional_habitat = "once-logged") %>%
#   rename(functionalhabAge = true_year) %>%
#   select(-time_since_intervention)
# 
# twice_df_30_shifted <- twice_df_30 %>%
#   mutate(true_year = (0:60) + 30) %>%
#   rename(functionalhabAge = functionalhabAge) %>%
#   filter(true_year < 61) %>%
#   mutate(functionalhabAge = functionalhabAge + 1) %>%
#   select(-true_year)
# 
# P_2L_df <- bind_rows(yr1_30_primary_2L, twice_df_30_shifted)
# 
# # 2L -> 2L starting typology
# L2_2L_df <- twice_df_starting %>%
#   select(-functional_habitat) %>% mutate(functional_habitat = "twice-logged") %>%
#   rename(functionalhabAge = functionalhabAge)
# 
# # 1L -> 2L (not allowed in first 15 yrs)
# # Build oncelogged subset > 14 & <45 (like your original)
# oncelogged_for_L1_L2 <- L1_R %>%
#   filter(functional_habitat == "once-logged") %>%
#   filter(functionalhabAge > 14, functionalhabAge < 45) %>%
#   mutate(functionalhabAge = functionalhabAge - 15,
#          original_habitat = "once-logged",
#          habitat = "twice-logged") %>%
#   select(-time_since_intervention)
# 
# twicelogged_from_P2L <- P_2L_df %>% filter(functional_habitat == "twice-logged") %>% mutate(original_habitat = "once-logged")
# L1_L2_df <- bind_rows(oncelogged_for_L1_L2, twicelogged_from_P2L)
# 
# # ---------- FINAL COMBINE ----------
# carbhabs <- carbhabs %>% bind_rows(P_2L_df) %>% bind_rows(L2_2L_df) %>% bind_rows(L1_L2_df)
# 
# # ---------- BELOWGROUND CARBON FUNCTION ----------
# belowground_fun <- function(x) {
#   plantations <- x %>%
#     filter(str_detect(habitat, "albizia|eucalyptus")) %>%
#     mutate(full_carbon = ACD,
#            full_carbon_lwr = lwr_ACD,
#            full_carbon_upr = upr_ACD)
#   
#   x2 <- x %>% filter(!str_detect(habitat, "albizia|eucalyptus"))
#   
#   yr1_filt <- x2 %>% filter(functionalhabAge == 1) %>%
#     select(functional_habitat, original_habitat, habitat, ACD, upr_ACD, lwr_ACD) %>%
#     rename(ACDt1 = ACD, uprACDt1 = upr_ACD, lwrACDt1 = lwr_ACD) %>% unique()
#   
#   y <- x2 %>% left_join(yr1_filt, by = c("functional_habitat", "original_habitat", "habitat")) %>%
#     mutate(
#       full_carbon = case_when(
#         original_habitat != functional_habitat & functionalhabAge <= 10 ~ ACDt1,
#         original_habitat == functional_habitat | functionalhabAge > 10 ~ ACD + (ACD * 0.31),
#         TRUE ~ NA_real_
#       ),
#       full_carbon_lwr = case_when(
#         original_habitat != functional_habitat & functionalhabAge <= 10 ~ lwr_ACD,
#         original_habitat == functional_habitat | functionalhabAge > 10 ~ lwr_ACD + (lwr_ACD * 0.31),
#         TRUE ~ NA_real_
#       ),
#       full_carbon_upr = case_when(
#         original_habitat != functional_habitat & functionalhabAge <= 10 ~ upr_ACD,
#         original_habitat == functional_habitat | functionalhabAge > 10 ~ upr_ACD + (upr_ACD * 0.31),
#         TRUE ~ NA_real_
#       )
#     ) %>% select(-c(ACDt1, uprACDt1, lwrACDt1))
#   
#   bind_rows(plantations, y)
# }
# 
# allHabCarbon_60yrACD <- belowground_fun(carbhabs)
# 
# allHabCarbon_60yrACD <- allHabCarbon_60yrACD %>%
#   select(original_habitat, habitat, functional_habitat,
#          functionalhabAge,
#          full_carbon, full_carbon_lwr, full_carbon_upr,
#          ACD, lwr_ACD, upr_ACD)
# 
# # ---------- PLOTS (same as before) ----------
# p_ACD <- allHabCarbon_60yrACD %>%
#   ggplot(aes(x = functionalhabAge, y = ACD, colour = functional_habitat)) +
#   geom_line() +
#   facet_wrap(~original_habitat + habitat, scales = "free_y") +
#   labs(x = "True Year", y = "Aboveground Carbon Density (ACD)") +
#   theme_minimal()
# 
# p_full <- allHabCarbon_60yrACD %>%
#   ggplot(aes(x = functionalhabAge, y = full_carbon, colour = functional_habitat)) +
#   geom_line() +
#   facet_wrap(~original_habitat + habitat, scales = "free_y") +
#   labs(x = "True Year", y = "Full Carbon") +
#   theme_minimal()
# 
# plot_grid(p_ACD, p_full, ncol = 1)
# 
# # ---------- OUTPUT ----------
# write.csv(allHabCarbon_60yrACD, "Outputs/allHabCarbon_60yrACD_bayesian.csv", row.names = FALSE)
# 
# # Print summaries
# cat("Primary (UnLogged) posterior mean ACD:", round(prime_mean, 1), "\n")
# cat("Primary 95% CI:", round(prime_ci[1], 1), "-", round(prime_ci[2], 1), "\n")
# cat("Once-logged slope (posterior mean):", round(slope_once_mean, 4), "Mg.ha-1 per year\n")
# cat("Once-logged slope 95% CI:", round(slope_once_lwr, 4), "-", round(slope_once_upr, 4), "\n")
