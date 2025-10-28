#2.10.25

#plantation carbon calculation

library(tidyverse)
library(terra)
library(ggpubr)
library(mgcv)
library(ggplot2)
library(brms)
library(bayesplot)
library(tidybayes)

#read in SABAH SOFTWOOD DATA 
SSB <- read.csv("RawData/SSB_plot_data.csv") %>%  rename(
  month = "Age..Month.", 
  volume_m3_ha = "Volume..m3.ha.", 
  dbh = "DBH..cm.", 
  height = "Height..m.", 
  num_trees_ha = "Actual.Trees.SPH", 
  live_trees_ha = "Live.Trees.SPH", 
  MAI = "MAI..m3.ha.yr." )  %>%  
  mutate(num_trees_ha = str_remove_all(num_trees_ha, ","),
         live_trees_ha = str_remove_all(live_trees_ha, ",")) %>% 
  mutate(num_trees_ha = as.numeric(num_trees_ha), 
         live_trees_ha = as.numeric(live_trees_ha)) 
names(SSB)
SSB$Species %>% unique()


#calculate AGB of tree plantations using chave et al equation ####

#WOOD DENSITY VALUES ARE FROM JARAPUDIN ET AL. 2020 (for 6 year Albizia and 6 year eucalyptus)

# Alb moluccana was 560 kg m−3 (0.56 g/cm3) and 252 (0.252 g/cm3) kg m−3 at six and five years of age, respectively, 
# E. pellita had a mean basic density of 629 (0.629 g/cm3) kg m−3 at six years of age.
(0.56+0.252)/2 #0.406

#define wood density 
EP_wd <- 0.629
AB_wd <- 0.406

SSB <- SSB %>% mutate(wood_density = case_when(
  Species == "EP" ~ EP_wd,
  Species == "AF" ~ AB_wd 
))

#CHAVE AGBest = 0.0673 x (Wood density x diametre^2 x height)^0.976 - tree level biomass

SSB <- SSB %>% 
  #biomass per tree
  mutate(tree_biomass = 0.0673*(wood_density*dbh^2* height)^0.976) %>% 
  #biomass per hectare (measured in kilogrames) 
  mutate(biomass_ha = tree_biomass * live_trees_ha) %>%  
  # mutate(biomass_ha = tree_biomass * num_trees_ha) %>%  
  #MgHa (measured in tonnes)
  mutate(AGB = biomass_ha/1000) %>%  
  #covert AGB into above-ground carbon density per hectare by multiplying by 0.47 [from Chave et al.]
  mutate(ACD_Mg_ha = AGB*0.47) %>%  
  mutate(YEAR = month/12)


# determine sampling coverage 
num_trees_ALB <- SSB %>% filter(Species == "AF") %>% filter(YEAR <13)
sum(num_trees_ALB$"Area..ha..")
sum(num_trees_ALB$"live_trees_ha") #76355

num_trees_EC <- SSB %>% filter(Species == "EP") %>% filter(YEAR <7)
sum(num_trees_EC$"Area..ha..")  #2879.25+  6201.2
sum(num_trees_EC$"live_trees_ha") #295275

#determine albizia best-fitting model ####



# -------------------------
# 1) Prepare data subsets
# -------------------------
AL_filtered <- SSB %>% filter(Species == "AF", YEAR < 13, !is.na(ACD_Mg_ha)) %>%    mutate(YEAR = as.numeric(YEAR))
EC_filtered <- SSB %>% filter(Species == "EP", YEAR < 13, !is.na(ACD_Mg_ha))  %>% mutate(YEAR = as.numeric(YEAR))

# -------------------------
# 2) Priors (example weakly informative)
# -------------------------
# Compute mean outside
mean_AL <- mean(AL_filtered$ACD_Mg_ha, na.rm = TRUE)
mean_EC <- mean(EC_filtered$ACD_Mg_ha, na.rm = TRUE)

print(mean_AL)
print(mean_EC)

priors_basic <- c(
  prior(normal(0, 10), class = "b"),
  prior(student_t(3, 0, 10), class = "sigma")
)

priors_basic_EC <- c(
  prior(normal(0, 10), class = "b"),
  prior(student_t(3, 0, 10), class = "sigma")
)
# -------------------------
# 3) Fit Bayesian models
# -------------------------
# Linear (Bayesian lm) for Albizia or Eucalyptus as needed
fit_AL_lm <- brm(
  formula = ACD_Mg_ha ~ YEAR +(1|Block.No),
  data = AL_filtered,
  family = gaussian(),
  prior = priors_basic,
  iter = 6000, warmup = 2000, adapt_delta = 0.95,
  backend = "cmdstanr"           # optional: recommended if installed
)

fit_EC_lm <- brm(
  formula = ACD_Mg_ha ~  YEAR +(1|Block.No),
  data = EC_filtered,
  family = gaussian(),
  prior = priors_basic_EC,
  iter = 6000, warmup = 2000, adapt_delta = 0.95,
  backend = "cmdstanr"
)

#After checking effects of spline model - linear is better
# #Spline model instead 
# # Fit spline model for Albizia
# fit_AL_spline <- brm(
#   formula = ACD_Mg_ha ~0 +  s(YEAR) + (1|Block.No),
#   data = AL_filtered,
#   family = gaussian(),
#   prior = prior(student_t(3, 0, 10), class = "sigma"),
#   chains = 4, iter = 4000, warmup = 1000,
#   backend = "cmdstanr", 
#   adapt_delta = 0.9
# )
# 
# unique(AL_filtered$YEAR)
# # Fit spline model for Eucalyptus
# fit_EC_spline <- brm(
#   formula = ACD_Mg_ha ~0 +  s(YEAR),
#   data = EC_filtered,
#   family = gaussian(),
#   prior = prior(student_t(3, 0, 10), class = "sigma"),
#   chains = 4, iter = 3000, warmup = 1000,
#   backend = "cmdstanr",
#   adapt_delta = 0.9
#   
# )
# 
# # Compute LOO for both models
# loo_AL_lm <- loo(fit_AL_lm) #AL linear model is preferred
# loo_AL_spline <- loo(fit_AL_spline)
# 
# loo_EC_lm <- loo(fit_EC_lm)
# loo_EC_spline <- loo(fit_EC_spline)
# 
# # Compare models
# loo_compare(loo_AL_lm, loo_AL_spline)
# loo_compare(loo_EC_lm, loo_EC_spline)

#run some basic pp checks

# Basic graphical checks

pp_check(fit_EC_lm, type = "dens_overlay", nsamples = 50) +
  ggtitle("Albizia: Posterior Predictive Check")
pp_check(fit_AL_lm, type = "dens_overlay", nsamples = 50) +
  ggtitle("Eucalyptus: Posterior Predictive Check")


# Check means and variances (all look good)
pp_check(fit_EC_lm, type = "stat", stat = "mean")
pp_check(fit_EC_lm, type = "stat", stat = "sd")

pp_check(fit_AL_lm, type = "stat", stat = "mean")
pp_check(fit_AL_lm, type = "stat", stat = "sd")


# -------------------------
# Make model predictions
# -------------------------
new_AL <- data.frame(YEAR = seq(0, 12, by = 0.5))
new_EC <- data.frame(YEAR = seq(0, 10, by = 0.5))

# Use posterior_epred to get draws on response scale (includes sigma noise)
ep_AL_lm     <- posterior_predict(fit_AL_lm,    newdata = new_AL, re_formula = NA)
ep_EC_lm     <- posterior_predict(fit_EC_lm,     newdata = new_EC,re_formula = NA)

# Summarise draws -> mean and 95% CI
summarise_draws <- function(ep_matrix, newdata) {
  sm <- t(apply(ep_matrix, 2, function(x) c(mean = mean(x),
                                            lwr = quantile(x, 0.025),
                                            upr = quantile(x, 0.975))))
  sm_df <- as.data.frame(sm)
  bind_cols(newdata, sm_df)
}

AL_lm_pred     <- summarise_draws(ep_AL_lm, new_AL) %>% 
  rename(lwr ='lwr.2.5%', 
         upr = 'upr.97.5%')

EC_lm_pred  <- summarise_draws(ep_EC_lm, new_EC)%>% 
  rename(lwr ='lwr.2.5%', 
         upr = 'upr.97.5%')


# -------------------------
#Plots: mean + 95% credible ribbon + observed points
# -------------------------

p_AL_lm <- ggplot(AL_lm_pred, aes(x = YEAR, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey80", alpha = 0.4) +
  geom_line(color = "black") +
  geom_point(data = AL_filtered, aes(x = YEAR, y = ACD_Mg_ha), alpha = 0.2) +
  labs(title = "Albizia falcataria", x = "Year", y = "ACD (Mg ha^-1)") +
  theme_minimal(base_size = 14)

p_EC_lm <- ggplot(EC_lm_pred, aes(x = YEAR, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey80", alpha = 0.4) +
  geom_line(color = "black") +
  xlim(0,9)+
  geom_point(data = AL_filtered, aes(x = YEAR, y = ACD_Mg_ha), alpha = 0.2) +
  labs(title = "Eucalyptus pellita", x = "Year", y = "ACD (Mg ha^-1)") +
  theme_minimal(base_size = 14)

plantation_plot <- cowplot::plot_grid(p_AL_lm,p_EC_lm)     

#get 500 posterior draws ###

# --- Define ages for prediction ---
ages_EC <- 0:6
ages_AL <- 0:12

# --- Get posterior predictions (fitted values) ---
# Albizia
AL_draws <- fit_AL_lm %>%
  add_epred_draws(
    newdata   = tibble(YEAR = ages_AL),
    re_formula = NA,     # exclude random effects
    ndraws     = 500     # request 500 draws
  ) %>% ungroup %>% 
  select(.draw, YEAR, .epred) %>%
  mutate(species = "albizia_current") %>%  
  rename(draw = .draw,
         plantationAge = YEAR, 
         ACD = .epred)

# Eucalyptus
EC_draws <- fit_EC_lm %>%
  add_epred_draws(
    newdata   = tibble(YEAR = ages_EC),
    re_formula = NA,     # exclude random effects
    ndraws     = 500     # request 500 draws
  ) %>% ungroup %>% 
  select(.draw, YEAR, .epred) %>%
  mutate(species = "eucalyptus_current") %>%  
  rename(draw = .draw,
         plantationAge = YEAR, 
         ACD = .epred)

# combine
plantation_draws <- bind_rows(AL_draws, EC_draws)

# inspect
plantation_draws %>% group_by(species) %>% summarise(unique_draws = n_distinct(draw))



#OUTPUT #### 
#plantation carbon draws 
saveRDS(plantation_draws, "Outputs/plantation_carbon_draws.rds")

#output figure
ggsave(
  filename = "Figures/plantation_ACD_growth.png",
  plot = plantation_plot,
  width = 7.5,      # width in inches — good for one-column Nature/Science layout
  height = 4,       # adjust height for aspect ratio
  dpi = 600,        # high resolution for publication
  bg = "white",     # ensures white background
  units = "in"      # specify units for clarity
)
