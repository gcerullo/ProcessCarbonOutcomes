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
#  MODEL v.1
# ============================================================

joint_model <- brm(
  bf(
    ACD ~ 1 + 
      state +                 # state-specific intercept differences
      state:time - time +     # state-specific slopes; and remove primary slope 
      (1 | plot_id),          # plot-level random intercept
    sigma ~ state             # state-specific residual variance
  ),
  data = all_data,
  family = gaussian(), #CHANGE TO STUDENT T IF NEEDED
  prior = c(
    
    # --------------------------------------------------------
    # PRIMARY FOREST (reference category)
    # --------------------------------------------------------
    # Prediction:
    #   ACD_primary ≈ 200 MgC ha⁻¹
    #   No temporal trend (time forced to 0)
    prior(normal(200, 25), class = "Intercept"),
    
    
    # --------------------------------------------------------
    # INTERCEPT DIFFERENCES (relative to primary)
    # --------------------------------------------------------
    
    # Once-logged:
    #   200 - 160 = 40 MgC at time = 0
    prior(normal(-160, 20),
          class = "b",
          coef = "stateonce_logged"),
    
    # Restored:
    #   200 - 150 = 50 MgC at time = 0
    prior(normal(-150, 20),
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
    
    # Primary (slope, reference level)
    #   Weak prior; no strong assumption
    prior(normal(0, 3),
          class = "b",
          coef = "stateprimary:time"),
    
    
    # --------------------------------------------------------
    # RANDOM EFFECTS
    # --------------------------------------------------------
    prior(student_t(3, 0, 10), class = "sd"),
    
    
    # --------------------------------------------------------
    # RESIDUAL VARIATION
    # --------------------------------------------------------
    prior(student_t(3, 0, 1), class = "b", dpar = "sigma")    
  ),
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

