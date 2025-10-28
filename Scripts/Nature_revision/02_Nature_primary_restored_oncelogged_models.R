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


#OUTPUTS ####
#save BRMS models 
saveRDS(log_brm, "Models/logged_restored_model.rds")
saveRDS(primary_brm, "Models/primary_model.rds")

