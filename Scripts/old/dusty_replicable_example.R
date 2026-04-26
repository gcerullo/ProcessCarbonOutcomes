#replicable code for dusty 
library(tidyverse)
library(data.table)

#function for making sure we can visualise the order of scenarios correctly
scenario_row_order_visualisation_fun <- function(x){
  x %>%  
    group_by(index, production_target, habitat, original_habitat, harvest_delay) %>%  
    arrange(true_year, .by_group = TRUE) %>%
    ungroup()
}

#read in a two single posterior draws of scenarios
draw1 <- readRDS("outputs/dusty_single_post_draw.rds") %>% scenario_row_order_visualisation_fun()
draw2 <- readRDS("outputs/dusty_single_post_draw2.rds") %>% scenario_row_order_visualisation_fun()

#combine two draws 
df <- draw1 %>% rbind(draw2)%>% scenario_row_order_visualisation_fun()

setDT(df)   # ensure data.table format
names(df)

#This is the maths and steps for 1 posterior draw. 

#==============================================================
# STEP 1 — Compute carbon per 10 km² parcel
#==============================================================
# ACD is Mg/ha; a 10 km² parcel = 1000 ha
df[, ACD_per_parcel := ACD * 1000]


#==============================================================
# STEP 2 — Handle staggered parcel allocation
#==============================================================

df[, ACD_10km2_stag := ACD_per_parcel * parcels_per_delay]


#==============================================================
# STEP 3 — Summarise to habitat-transition × year
#==============================================================
# Group by: index, production_target, original_habitat, habitat, true_year

habitat_year_summary <- df[, .(
  
  #sum across the staggered schedule for each harvest transition
  hab_ACD_year = sum(ACD_10km2_stag, na.rm = TRUE),
  scenarioName = first(scenarioName),
  scenarioStart = first(scenarioStart)
  
), by = .(
  index,
  production_target,
  original_habitat,
  habitat,
  true_year,
  draw
)]


#==============================================================
# STEP 4 — Summarise to scenario × year
#==============================================================
#sum across transitions to calculate yearly carbon across landscape
scenario_year_summary <- habitat_year_summary[, .(
  
  scen_ACD_year = sum(hab_ACD_year, na.rm = TRUE),
  scenarioName = first(scenarioName),
  scenarioStart = first(scenarioStart)
  
), by = .(
  index,
  production_target,
  true_year,
  draw
)]


#==============================================================
# STEP 5 — Summarise cumulative carbon across all years
#==============================================================
# This gives total cumulative ACD over the full simulation period.

carbon_summary <- scenario_year_summary[, .(
  
  cumulative_stock_year = sum(scen_ACD_year, na.rm = TRUE),
  scenarioName = first(scenarioName),
  scenarioStart = first(scenarioStart)
  
), by = .(
  production_target,
  index,
  draw
)]


#summarise across the two posterior draws 
carbon_stock_years_summary <- carbon_summary[, .(
  
  # Mean of ACD per across draws
  mean_cum_stock_year = mean(cumulative_stock_year, na.rm = TRUE),
  
  # 95% credible interval
  lwr_cum_stock_year_95 = quantile(cumulative_stock_year, probs = 0.025, na.rm = TRUE),
  upr_cum_stock_year_95 = quantile(cumulative_stock_year, probs = 0.975, na.rm = TRUE),
  
  # 80% credible interval
  lwr_cum_stock_year_80 = quantile(cumulative_stock_year, probs = 0.10, na.rm = TRUE),
  upr_cum_stock_year_80 = quantile(cumulative_stock_year, probs = 0.90, na.rm = TRUE),
  
  lwr_cum_stock_year_50 = quantile(cumulative_stock_year, probs = 0.25, na.rm = TRUE),
  upr_cum_stock_year_50 = quantile(cumulative_stock_year, probs = 0.75, na.rm = TRUE),
  
  # Optional metadata (if constant across draws)
  scenarioName = first(scenarioName),
  scenarioStart = first(scenarioStart)
  
), by = .(index, production_target)]


#plot 

plot <- carbon_stock_years_summary %>%
  ggplot(aes(x = production_target,
             y = mean_cum_stock_year)) +
  
  # Error bars (clearer, thicker, slightly transparent)
  geom_errorbar(
    aes(ymin = lwr_cum_stock_year_95,
        ymax = upr_cum_stock_year_95),
    width = 0.015,
    linewidth = 0.6,
    alpha = 0.5
  ) +
  
  # Points (larger, semi-transparent, with slight jitter)
  geom_point(
    position = position_jitter(width = 0.015, height = 0),
    size = 2.2,
    alpha = 0.8
  ) +
  

  # Axes and labels
  xlim(0, 1) +
  xlab("Production target") +
  ylab("Total aboveground carbon (Mg C)") +
  
  facet_wrap(~scenarioName, ncol = 4) +
  
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    strip.text = element_blank()
  )

plot

