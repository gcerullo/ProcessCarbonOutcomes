#04.06.24 This code reads 

library(tidyverse)
library(terra)
library(ggpubr)
library(mgcv)
library(ggplot2)

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

#Come back and change this equations!!!! 1/10/2025
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
# 
# #BAYES FRAMEWORK
# # subset species
# AL_filtered <- SSB %>%
#   filter(Species == "AF", YEAR < 13)
# 
# EC_filtered <- SSB %>%
#   filter(Species == "EP", YEAR < 13, !is.na(ACD_Mg_ha))
# 
# # ---- Example 1: linear model ----
# fit_lm <- brm(
#   ACD_Mg_ha ~ YEAR,
#   data = AL_filtered,
#   family = gaussian(),
#   chains = 4, cores = 4, iter = 2000
# )
# 
# # ---- Example 2: spline model (Bayesian GAM) ----
# fit_spline <- brm(
#   ACD_Mg_ha ~ s(YEAR),
#   data = EC_filtered,
#   family = gaussian(),
#   chains = 4, cores = 4, iter = 4000
# )
# 
# # ---- Predictions with uncertainty ----
# newdat <- data.frame(YEAR = seq(0, 13, length.out = 100))
# 
# pred_spline <- posterior_epred(fit_spline, newdata = newdat) # posterior draws
# pred_spline_mean <- apply(pred_spline, 2, mean)
# pred_spline_ci <- apply(pred_spline, 2, quantile, probs = c(0.025, 0.975))
# 
# plot_df <- cbind(newdat, 
#                  mean = pred_spline_mean,
#                  lwr = pred_spline_ci[1,],
#                  upr = pred_spline_ci[2,])
# 
# ggplot(plot_df, aes(x = YEAR, y = mean)) +
#   geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey80", alpha = 0.4) +
#   geom_line(color = "black") +
#   geom_point(data = EC_filtered, aes(x = YEAR, y = ACD_Mg_ha), alpha = 0.2) +
#   labs(x = "Plantation age", 
#        y = expression("Aboveground Carbon Density (Mg ha"^-1*")")) +
#   theme_pubr(base_size = 16)

#
# # Fit the lm model
#
# AL_filtered <- SSB %>% filter(Species == "AF") %>% filter(YEAR <13)
# model_lm <- lm(ACD_Mg_ha ~ YEAR, data = AL_filtered)
# predictions_lm <- predict(model_lm)
# AL_filtered %>% select(Block.No) %>% unique() %>% count()
#
# # Fit the loess model
# model_loess <- loess(ACD_Mg_ha ~ YEAR, data = AL_filtered)
# predictions_loess <- predict(model_loess)
#
# # Fit the spline model
# model_spline <- gam(ACD_Mg_ha ~ s(YEAR), data = AL_filtered)
# predictions_spline <- predict(model_spline)
#
# # Calculate MSE for lm, loess, and spline
# mse_lm <- mean((AL_filtered$ACD_Mg_ha %>% na.omit - predictions_lm)^2)
# mse_loess <- mean((AL_filtered$ACD_Mg_ha %>% na.omit - predictions_loess)^2)
# mse_spline <- mean((AL_filtered$ACD_Mg_ha %>% na.omit - predictions_spline)^2)
#
# # Compare the MSE values
# lowest_mse <- min(mse_lm, mse_loess, mse_spline)
# best_model <- switch(match(lowest_mse, c(mse_lm, mse_loess, mse_spline)), "lm", "loess", "spline")
#
# # Print the MSE values and the best model   #LOESS IS BEST
# cat("MSE (lm):", mse_lm, "\n")
# cat("MSE (loess):", mse_loess, "\n")
# cat("MSE (spline):", mse_spline, "\n")
# cat("Best Model:", best_model, "\n")
# cat("Lowest MSE:", lowest_mse, "\n")
#
# # determine eucalyptus best-fitting model ####
#
# # Fit the lm model
# EC_filtered <- SSB %>% filter(Species == "EP") %>% filter(YEAR <13) %>% filter(!is.na(ACD_Mg_ha))
# model_lm <- lm(ACD_Mg_ha ~ YEAR, data = EC_filtered)
# predictions_lm <- predict(model_lm)
#
# # Fit the loess model
# model_loess <- loess(ACD_Mg_ha ~ YEAR, data = EC_filtered)
# predictions_loess <- predict(model_loess)
#
# # Fit the spline model
# model_spline <- gam(ACD_Mg_ha ~ s(YEAR), data = EC_filtered)
# predictions_spline <- predict(model_spline)
# confidence_interval <- predict(model_spline, se.fit = TRUE)
# upper_bound <- predictions_spline + 1.96 * confidence_interval$se.fit
# lower_bound <- predictions_spline - 1.96 * confidence_interval$se.fit
#
#
# # Calculate MSE for lm, loess, and spline
# mse_lm <- mean((EC_filtered$ACD_Mg_ha %>% na.omit() - predictions_lm)^2)
# mse_loess <- mean((EC_filtered$ACD_Mg_ha %>% na.omit() - predictions_loess)^2)
# mse_spline <- mean((EC_filtered$ACD_Mg_ha %>% na.omit - predictions_spline)^2)
#
# # Compare the MSE values
# lowest_mse <- min(mse_lm, mse_loess, mse_spline)
# best_model <- switch(match(lowest_mse, c(mse_lm, mse_loess, mse_spline)), "lm", "loess", "spline")
#
# # Print the MSE values and the best model   #LOESS IS BEST
# cat("MSE (lm):", mse_lm, "\n")
# cat("MSE (loess):", mse_loess, "\n")
# cat("MSE (spline):", mse_spline, "\n")
# cat("Best Model:", best_model, "\n")
# cat("Lowest MSE:", lowest_mse, "\n")
#
#
# #Plot data ####
#
# #automatically plots 95% confidence intervals
# sample_coverage <- AL_filtered %>% rbind(EC_filtered) %>% rename(area = "Area..ha..")
# sample_coverage %>% summarise(meanArea = mean(area),
#                               sd(area))
#
# names(sample_coverage)
# ec <- EC_filtered %>%
#   ggplot(aes(x = YEAR, y = ACD_Mg_ha)) +
#   geom_point(alpha = 0.1) +
#   geom_smooth(method = "gam", formula = y ~ s(x), colour = "black") +
#   geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.3, fill = "grey80") +
#   labs(x = "Plantation age", y = "Aboveground Carbon Density (Mg ha\u207B\u00B9)") +
#   ggtitle("Eucalyptus pellita") +
#   theme_pubr(base_size = 16)
#
#
#
#
# #albizia (loess)
# al <- AL_filtered %>%
#   filter(YEAR < 13.5) %>%
#   ggplot(aes(x = YEAR, y = ACD_Mg_ha)) +
#   geom_point(alpha = 0.1) +
#   geom_smooth(method = "loess", colour = "black") +
#   labs(x = "Plantation age", y = "Aboveground Carbon Density (Mg ha\u207B\u00B9)") +
#   ggtitle("  Albizia falcataria")+
#   theme_pubr(base_size =16)
#
# plant_figures <- cowplot::plot_grid(ec,al)


#----export figures ------

#set figure export path
path = "Figures/"
# Set the dimensions for A4 size in inches
width <- 8.27
height <- 11.69

ggsave(plant_figures,
       filename = paste0(path, "//PlantationACD.pdf"),
       width =  width, #in pixels
       height = height/2,
       units = "in")


#EXTRACT plantation values 
#(mean, and 5/95 confidence intervals) for a given age of plantation 


#define the plantation age you want to extract mean and error for
age_to_extract_AF <- seq(1,12, by = 1)
age_to_extract_EC <- seq(1,6, by =1)

#1. For eucalytpus, linear relationship is fine becuase plantations are harvested during the linear growth phase. 
#2 For albizia, carbon over time varies; so should I fit loess or linear or spline relationship? 


#FIT FINAL MODELS AND PLOT FIGURES ####


#fit loess for albizia
model_AF <- loess(ACD_Mg_ha ~ YEAR, data = AL_filtered) 
AF_extract <- data.frame(YEAR = age_to_extract_AF)
AF_prediction <- AF_extract %>%
  mutate(ACD = predict(model_AF, newdata = AF_extract),
         ACD_se = predict(model_AF, newdata = AF_extract, se = TRUE)$se) %>%
  cbind(habitat = "albizia") %>%
  # multiply SE by 1.96 to get the 5th and 95th percentile
  mutate(upr_ACD = ACD + ACD_se * 1.96,
         lwr_ACD = ACD - ACD_se * 1.96) %>%
  select(habitat, ACD, lwr_ACD, upr_ACD) %>%
  cbind(plantationAge = age_to_extract_AF)


#fit lm for albizia
model_EC <- lm(ACD_Mg_ha ~ YEAR, data = EC_filtered)
EC_extract <- data.frame(YEAR = age_to_extract_EC)
EC_prediction <- as.data.frame(predict(model_EC, newdata = EC_extract, type = 'response', se.fit = TRUE))  %>% 
  cbind(habitat = "eucalyptus")  %>%  
  rename(ACD = fit, 
         ACD_se = se.fit) %>% 
  #multiply SE by 1.96 to get the 5 and 95th percentile
  mutate(upr_ACD=  ACD + ACD_se * 1.96,
         lwr_ACD= ACD - ACD_se * 1.96) %>%  
  select(habitat, ACD,lwr_ACD, upr_ACD) %>% 
  cbind(plantationAge = age_to_extract_EC )


# final plantation carbon results ####
plantation_carbon <- rbind(AF_prediction, EC_prediction)

write.csv(plantation_carbon, "Outputs/plantation_carbon.csv")

