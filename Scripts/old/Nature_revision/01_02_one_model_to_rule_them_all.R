# using one big fat models 

# ============================================================
# 0) LIBRARIES
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)

set.seed(123)

# ============================================================
# 1) PREPARE PLANTATION DATA
# ============================================================

SSB <- read.csv("RawData/SSB_plot_data.csv") %>%
  rename(
    month = Age..Month.,
    dbh = DBH..cm.,
    height = Height..m.,
    live_trees_ha = Live.Trees.SPH
  ) %>%
  mutate(
    live_trees_ha = as.numeric(str_remove_all(live_trees_ha, ",")),
    YEAR = month / 12
  )

# Wood densities
EP_wd <- 0.629
AB_wd <- 0.406

SSB <- SSB %>%
  mutate(
    wood_density = case_when(
      Species == "EP" ~ EP_wd,
      Species == "AF" ~ AB_wd
    ),
    tree_biomass = 0.0673 * (wood_density * dbh^2 * height)^0.976,
    biomass_ha = tree_biomass * live_trees_ha,
    AGB = biomass_ha / 1000,
    ACD = AGB * 0.47
  )

plantation <- SSB %>%
  filter(!is.na(ACD), YEAR < 15) %>%
  mutate(
    state = case_when(
      Species == "EP" ~ "plantation_EP",
      Species == "AF" ~ "plantation_AF"
    ),
    time = YEAR,
    plot_id = Block.No
  ) %>%
  select(ACD, time, state, plot_id)

# ============================================================
# 2) PREPARE LOGGED + PRIMARY DATA
# ============================================================

data_phil <- read.csv("RawData/Philipson20_PlotData2.csv")

Logged <- data_phil %>%
  filter(Forest == "Logged") %>%
  mutate(
    state = case_when(
      FACE == "ProjectScenario" ~ "restored",
      TRUE ~ "once_logged"
    ),
    time = YearsSinceLogging,
    plot_id = Plot
  ) %>%
  select(ACD, time, state, plot_id)

Primary <- data_phil %>%
  filter(Forest == "UnLogged") %>%
  mutate(
    state = "primary",
    time = 0,   # primary treated as baseline state
    plot_id = Plot
  ) %>%
  select(ACD, time, state, plot_id)

# ============================================================
# 3) COMBINE ALL DATA
# ============================================================


all_data <- bind_rows(
  plantation,
  Logged,
  Primary
)



# ============================================================
# PREP DATA
# ============================================================


# Ensure primary is reference level
all_data$state <- relevel(factor(all_data$state), ref = "primary")
all_data$plot_id <- factor(all_data$plot_id)

# Primary has no temporal variation
all_data <- all_data %>%
  mutate(
    time_state = ifelse(state == "primary", 0, time))

# ============================================================
# PRIOR PREDICTIVE CHECKS FOR THE JOINT CARBON MODEL
# ============================================================
# Goal:
# Evaluate whether the priors specified for the carbon recovery
# model generate realistic ecological outcomes BEFORE fitting
# the model to data.
#
# Why this matters:
# If priors are too tight they can strongly constrain inference.
# If priors are too wide they can produce unrealistic trajectories
# (e.g. negative carbon or >1000 MgC ha-1 forests).
#
# Prior predictive simulation samples from the priors and generates
# simulated carbon values WITHOUT using the observed data likelihood.
#
# We check:
# 1) Distribution of simulated carbon values
# 2) Carbon trajectories through time
# 3) Implied starting carbon values
# 4) Implied recovery slopes
# 5) Whether plausible ecological ranges are covered
#
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)

set.seed(123)

# ============================================================
# PREP DATA (same as main model)
# ============================================================

all_data$state <- relevel(factor(all_data$state), ref = "primary")
all_data$plot_id <- factor(all_data$plot_id)

all_data <- all_data %>%
  mutate(
    time_state = ifelse(state == "primary", 0, time)
  )

# ============================================================
# BUILD PREDICTION GRID
# ============================================================
# This grid represents the ecological space we want to explore
# under the priors.

prediction_grid <- expand.grid(
  time  = 0:75,
  state = levels(all_data$state)
) %>%
  dplyr::filter(!(state == "primary" & time != 0))

# ============================================================
# PRIOR-ONLY MODEL - test the impacts of the prior selections a priori
# ============================================================
# sample_prior = "only" tells brms:
# - Ignore the likelihood entirely
# - Draw ALL parameters from the specified priors
# - Generate implied (prior predictive) ACD values purely from assumptions
#
# CONSEQUENCE:
# Everything below defines a *generative model* of forest carbon dynamics
# BEFORE seeing any data. The outputs represent what my ecological
# assumptions alone imply about ACD levels and recovery trajectories.


# ============================================================
# PRIOR-ONLY MODEL — Ecologically constrained, weakly informative variance
# ============================================================
# PURPOSE:
# - Generate realistic prior predictive carbon trajectories
# - Encode strong ecological knowledge where justified (means/slopes)
# - Keep variance components weakly informative (not dominating signal)

prior_model <- brm(
  bf(
    ACD ~ 1 + 
      state +
      state:time - time +     # No global slope
      # Each state has its own independent recovery rate
      # CONSEQUENCE:
      # - No shrinkage toward a shared recovery trend
      # - Differences between states are fully prior-driven
      
      (1 | plot_id),          # Plot-level random intercept
    # CONSEQUENCE:
    # - Captures heterogeneity among plots within states
    # - Adds vertical spread around trajectories
    
    sigma ~ state             # State-specific residual variance
    # CONSEQUENCE:
    # - Different forest types can have different noise levels
  ),
  
  data = all_data,
  family = gaussian(),
  
  prior = c(
    
    # =========================================================
    # PRIMARY FOREST BASELINE (ANCHOR OF MODEL)
    # =========================================================
    # Mean ACD at time = 0 for primary forest
    
    prior(normal(200, 30), class = "Intercept"),
    
    # RATIONALE:
    # - Primary forests are typically ~150–250 MgC/ha
    # - Slightly relaxed from SD = 25 → allows broader realism
    #
    # CONSEQUENCE:
    # - Anchors all other states (defined relative to primary)
    # - Prevents extreme low (<100) or high (>350) values
    # - Keeps trajectories in ecologically plausible range
    
    
    # =========================================================
    # INITIAL CARBON DIFFERENCES (STATE OFFSETS)
    # =========================================================
    # These define ACD at time = 0 relative to primary
    
    # Once logged (~40 MgC/ha)
    prior(normal(-160, 20),
          class = "b",
          coef = "stateonce_logged"),
    
    # CONSEQUENCE:
    # - Most values ~0–80 MgC/ha
    # - Strong depletion relative to primary
    # - Allows near-zero but rarely negative
    
    
    # Restored (~50 MgC/ha)
    prior(normal(-150, 20),
          class = "b",
          coef = "staterestored"),
    
    # CONSEQUENCE:
    # - Slightly higher than once logged
    # - Overlap reflects uncertainty in restoration outcomes
    
    
    # Plantations (~0 MgC/ha)
    prior(normal(-200, 20),
          class = "b",
          coef = "stateplantation_AF"),
    
    prior(normal(-200, 20),
          class = "b",
          coef = "stateplantation_EP"),
    
    # RATIONALE:
    # - Tightened from SD = 30 → reduces negative ACD
    #
    # CONSEQUENCE:
    # - Most plantation values near 0–50 MgC/ha
    # - Rarely negative
    # - Maintains strong contrast with natural forests
    
    
    # =========================================================
    # RECOVERY SLOPES (KEY ECOLOGICAL PROCESS)
    # =========================================================
    # Units: MgC ha⁻¹ year⁻¹
    
    # Once logged recovery (moderate)
    prior(normal(2, 1.5),
          class = "b",
          coef = "stateonce_logged:time"),
    
    # CONSEQUENCE:
    # - Gradual recovery
    # - Allows some decline, but mostly positive
    
    
    # Restored recovery (faster)
    prior(normal(4, 2),
          class = "b",
          coef = "staterestored:time"),
    
    # CONSEQUENCE:
    # - Strong upward trajectories
    # - Encodes expectation that restoration accelerates recovery
    
    
    # Plantation slopes (uncertain)
    prior(normal(0, 2),
          class = "b",
          coef = "stateplantation_AF:time"),
    
    prior(normal(0, 2),
          class = "b",
          coef = "stateplantation_EP:time"),
    
    # RATIONALE:
    # - Reduced from SD = 3 → avoids extreme trajectories
    #
    # CONSEQUENCE:
    # - Allows:
    #   - growth
    #   - stagnation
    #   - mild decline
    # - Prevents explosive or highly negative trends
    
    
    # Primary forest slope (FORCED STABILITY)
    prior(normal(0, 0.2),
          class = "b",
          coef = "stateprimary:time"),
    
    # RATIONALE:
    # - Primary forests assumed near equilibrium
    #
    # CONSEQUENCE:
    # - Flat trajectories over time
    # - Removes unrealistic increases/decreases
    # - Critical for ecological realism
    
    
    # =========================================================
    # RANDOM EFFECT VARIATION (WEAKLY INFORMATIVE)
    # =========================================================
    
    prior(normal(0, 10), class = "sd"),
    
    # RATIONALE:
    # - Allows realistic between-plot variation (~±20 MgC/ha)
    # - Avoids extreme heavy-tailed behaviour
    #
    # CONSEQUENCE:
    # - Adds spread without dominating trajectories
    # - Prevents unrealistic extreme forests
    
    
    # =========================================================
    # RESIDUAL ERROR (WEAKLY INFORMATIVE)
    # =========================================================
    
    # Baseline residual SD
    prior(normal(0, 1), class = "sigma"),
    
    # State-specific deviations in residual variance
    prior(normal(0, 0.5), class = "b", dpar = "sigma")
    
    # RATIONALE:
    # - Keeps noise present but controlled
    # - Allows some states to be noisier than others
    #
    # CONSEQUENCE:
    # - Predictive intervals are realistic
    # - Noise does not overwhelm recovery trends
  ),
  
  
  # =========================================================
  # PRIOR SAMPLING
  # =========================================================
  
  sample_prior = "only",
  
  # CONSEQUENCE:
  # - Model ignores data entirely
  # - Outputs reflect ONLY prior assumptions
  # - Use posterior_predict() to visualise trajectories
  
  
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  backend = "cmdstanr",
  seed = 123
)# ============================================================
# 1. BASIC PRIOR PREDICTIVE DISTRIBUTION
# ============================================================

pp_check(prior_model)

# What to look for:
#
# Good signs:
# • Most values between ~0 and ~400 MgC/ha
#
# Warning signs:
# • Many negative carbon values
# • Large density >500 MgC/ha
# • Extremely wide distributions


# ============================================================
# 2. PRIOR TRAJECTORIES THROUGH TIME
# ============================================================

prior_draws <- prior_model %>%
  add_epred_draws(
    newdata = prediction_grid,
    re_formula = NA,
    ndraws = 300
  )

ggplot(prior_draws,
       aes(x = time, y = .epred, group = .draw)) +
  geom_line(alpha = 0.05) +
  facet_wrap(~ state) +
  theme_minimal() +
  ylim(-100, 500) +
  labs(
    y = "Simulated Carbon (Mg C ha⁻¹)",
    x = "Years",
    title = "Prior Predictive Carbon Trajectories"
  )

#primary forest 
prior_draws %>% filter(state == "primary") %>% 
ggplot(
       aes(x = time, y = .epred, group = .draw)) +
  geom_point(alpha = 0.05) +
  facet_wrap(~ state) +
  theme_minimal() +
  ylim(-100, 500) +
  labs(
    y = "Simulated Carbon (Mg C ha⁻¹)",
    x = "Years",
    title = "Prior Predictive Carbon Trajectories"
  )
# What to look for:
#
# Good:
# • Logged forests gradually recover
# • Plantations start near zero
# • Primary stays near ~200
#
# Bad:
# • Exploding trajectories (>600 MgC/ha)
# • Strong negative trends
# • Implausibly fast recovery (>10 MgC/ha/year)


# ============================================================
# 3. PRIOR DISTRIBUTION OF STARTING CARBON
# ============================================================

prior_draws %>%
  filter(time == 0) %>%
  ggplot(aes(x = .epred)) +
  geom_density(fill = "grey70") +
  facet_wrap(~ state, scales = "free") +
  theme_minimal() +
  labs(
    x = "Carbon at Time 0",
    title = "Prior Distribution of Initial Carbon"
  )

# Things to check:
#
# Primary:
# ~150–250
#
# Logged:
# ~20–80
#
# Plantations:
# near zero
#
# If distributions are extremely narrow,
# your priors may be overly informative.


# ============================================================
# 5. PRIOR SUMMARY TABLE
# ============================================================

prior_draws %>%
  group_by(state, time) %>%
  summarise(
    mean = mean(.epred),
    lower = quantile(.epred, 0.025),
    upper = quantile(.epred, 0.975),
    .groups = "drop"
  ) %>%
  filter(time %in% c(0, 25, 50, 75))

# This shows the range of carbon implied by your priors
# at different times.


# ============================================================
# INTERPRETING RESULTS
# ============================================================

# If priors are TOO TIGHT:
# • posterior trajectories will closely follow prior expectations
# • credible intervals are very narrow
# • prior trajectories look nearly identical
#
# If priors are TOO WIDE:
# • simulated carbon >700 MgC/ha
# • negative forests (<0 MgC/ha)
# • extremely steep slopes (>10 MgC/ha/yr)
#
# Ideal priors:
# • produce trajectories spanning realistic ecological ranges
# • but not implausible extremes.


# ============================================================
# AFTER THIS CHECK
# ============================================================
# If priors look reasonable, fit the real model:

# joint_model <- brm(...)

# ============================================================


# ============================================================
#  MODEL v.1
# ============================================================

joint_model <- brm(
  bf(
    ACD ~ 1 + 
      state +                 # state-specific intercept differences
      state:time - time +     # state-specific slopes; and remove primary slope 
      (1 | plot_id),          # plot-level random intercept
    sigma ~ 1                #  homoskedastic residual error
  ),
  data = all_data,
  family = gaussian(), #CHANGE TO STUDENT T IF NEEDED
  prior = c(
    
    # --------------------------------------------------------
    # PRIMARY FOREST (reference category)
    # --------------------------------------------------------
    # Prediction:
    #   ACD_primary ≈ 100-300 MgC ha⁻¹
    #   No temporal trend (time forced to 0)
    prior(normal(200, 50), class = "Intercept"),
    
    
    # --------------------------------------------------------
    # INTERCEPT DIFFERENCES (relative to primary)
    # --------------------------------------------------------
    #30–130 MgC starting values
    # Once-logged:
    #   200 - 150 = 50 MgC at time = 0
    prior(normal(-150, 40),
          class = "b",
          coef = "stateonce_logged"),
    #~30–130 MgC starting values
    # Restored:
    #   200 - 150 = 50 MgC at time = 0
    prior(normal(-150, 40),
          class = "b",
          coef = "staterestored"),
    
    # Plantation AF:
    #   200 - 200 = 0 MgC at time = 0
    prior(normal(-200, 30),
          class = "b",
          coef = "stateplantation_AF"),
    
    # Plantation EP:
    #   200 - 200 = 0 MgC at time = 0
    prior(normal(-200, 30),
          class = "b",
          coef = "stateplantation_EP"),
    
    
    # --------------------------------------------------------
    # STATE-SPECIFIC RECOVERY SLOPES (MgC ha⁻¹ yr⁻¹)
    # --------------------------------------------------------
    
    # Once-logged:
    #   Expected ≈ 2 MgC ha⁻¹ yr⁻¹
    #   Trajectory ≈ 40 + 2 * time
    prior(normal(2, 0.5),
          class = "b",
          coef = "stateonce_logged:time"),
    
    # Restored:
    #   Expected ≈ 4 MgC ha⁻¹ yr⁻¹
    #   Trajectory ≈ 50 + 4 * time
    prior(normal(4, 1.5),
          class = "b",
          coef = "staterestored:time"),
    
    # Plantation AF:
    #   Weak prior centered at 0
    #   Data determine recovery
    prior(normal(0, 3),
          class = "b",
          coef = "stateplantation_AF:time"),
    
    # Plantation EP:
    #   Weak prior centered at 0
    prior(normal(0, 3),
          class = "b",
          coef = "stateplantation_EP:time"),
    
    
    
    # --------------------------------------------------------
    # RANDOM EFFECTS
    # --------------------------------------------------------
    prior(student_t(3, 0, 10), class = "sd"),
    
    
    # --------------------------------------------------------
    # RESIDUAL VARIATION
    # --------------------------------------------------------
    prior(student_t(3, 0, 30), class = "Intercept", dpar = "sigma")),
    chains = 4,
  iter = 6000,
  warmup = 2000,
  cores = 4,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.99),
  seed = 123
)


saveRDS(joint_model, "Models/unified_linear_carbon_model_24.02_26.rds")

# ============================================================
#  CHECK MODEL FIT
# ============================================================

# List of states
states <- levels(all_data$state)


# Open a PDF device
pdf("pp_checks_by_state_v2.pdf", width = 8, height = 6)

for(s in states){
  state_data <- all_data %>% filter(state == s)
  
  # Default density overlay
  p1 <- pp_check(joint_model, newdata = state_data)
  print(p1 + ggtitle(paste("Density overlay -", s)))
  
  # Histogram
  p2 <- pp_check(joint_model, newdata = state_data, type = "hist")
  print(p2 + ggtitle(paste("Histogram -", s)))
  
  # Scatter observed vs predicted
  p3 <- pp_check(joint_model, newdata = state_data, type = "scatter_avg")
  print(p3 + ggtitle(paste("Scatter observed vs predicted -", s)))
  
  # Distribution of predicted means
  p4 <- pp_check(joint_model, newdata = state_data, type = "stat", stat = "mean")
  print(p4 + ggtitle(paste("Predicted means -", s)))
  
  # Distribution of predicted SDs
  p5 <- pp_check(joint_model, newdata = state_data, type = "stat", stat = "sd")
  print(p5 + ggtitle(paste("Predicted SDs -", s)))
}

# Close the PDF device
dev.off()

# # For Logged model
# pp_check(joint_model)  # default density overlay
# pp_check(joint_model, type = "hist")   # histograms of observed vs simulated
# pp_check(joint_model, type = "scatter_avg")  # scatter of observed vs avg predicted
# pp_check(joint_model, type = "stat", stat = "mean")  # distribution of predicted means
# pp_check(joint_model, type = "stat", stat = "sd")    # distribution of predicted SDs
# ============================================================
# POPULATION-LEVEL POSTERIOR TRAJECTORIES (NO RANDOM EFFECTS)
# ============================================================

#can start here
joint_model <-  readRDS("Models/unified_linear_carbon_model_24.02_26.rds")

# 1) Build prediction grid consistent with fitted model
all_data <- all_data %>% mutate(state = factor(state))
newdata <- expand.grid(
  time  = 0:75,
  state = levels(all_data$state)
) %>%
  dplyr::filter(!(state == "primary" & time != 0))

# 2) Get posterior expected values (population-level only)
draws_long <- joint_model %>%
  add_epred_draws(
    newdata = newdata,
    re_formula = NA,   # removes plot-level random intercepts
    ndraws = 500
  ) %>%
  rename(
    draw = .draw,
    ACD  = .epred
  ) 


# ------------------------------------------------------------
# SPAGHETTI PLOT
# ------------------------------------------------------------

ggplot(draws_long,
       aes(x = time, y = ACD, group = draw)) +
  geom_line(alpha = 0.05) +
  facet_wrap(~ state) +
  theme_minimal() +
  labs(
    y = "Aboveground Carbon (Mg C ha⁻¹)",
    x = "Years",
    title = "Population-Level Posterior Carbon Trajectories"
  )

# ------------------------------------------------------------
# MEAN + 95% CREDIBLE INTERVAL
# ------------------------------------------------------------

summary_df <- draws_long %>%
  group_by(state, time) %>%
  summarise(
    mean  = mean(ACD),
    lower = quantile(ACD, 0.025),
    upper = quantile(ACD, 0.975),
    .groups = "drop"
  )

ggplot(summary_df,
       aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2) +
  geom_line(linewidth = 1) +
  facet_wrap(~ state) +
  theme_minimal() +
  ylim(0,300)+
  labs(
    y = "Aboveground Carbon (Mg C ha⁻¹)",
    x = "Years",
    title = "Population-Level Predicted Carbon (Mean ± 95% CI)"
  )


#_______________________________________
#Save models and draws 
#_______________________________________
draws_long <- draws_long %>% ungroup %>%  
  select(-c( ".row"  ,     ".chain"  ,   ".iteration")) %>%  
   mutate(
    state = case_when(
      state == "plantation_EP" ~ "eucalyptus_current",
      state == "plantation_AF" ~ "albizia_current",
      TRUE ~ state
    )
  ) %>%
  rename(habitat = state)


saveRDS(draws_long, "Outputs/onemodel_ACD_draws.rds")

saveRDS(joint_model, "Models/unified_linear_carbon_model_24.02_26.rds")

