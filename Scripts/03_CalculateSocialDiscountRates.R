#09.06.24
#this code genereates estimates of the social cost of carbon through time, which I can subsequently use 
# to calculate carbon impact of scenarios by using discounted social cost of carbon values. 

library(tidyverse)
library(magrittr)
library(gridExtra)

#define params ####
rm(list = ls())

discount_rates = c(0.02,0.04, 0.06)
n_years<-60
rcp<-'RCP2.6'   
marginal_damage_exponent<-1

#read in data ####

rcps<-read_csv("Inputs/Venmans_Value of an offset 25_1evening.csv", skip = 2, col_types = 'd') %>%
  slice(-(1)) %>%
  select(starts_with('RCP'))

# either 2.6, 4.5, 6.0, 8.5
temp_rise<-rcps %>% pull(rcp) %>% as.numeric()

scc_gen<-function(temp_rise, r_discount = 0.03, r_GDP_growth = 0.017, timeframe = 1500, n_years = length(temp_rise),
                  marginal_damage_exponent = 1, GDP_start = 85000){
 
   start_year= 2020
  # len_tr<-length(temp_rise)
  # if(n_years<len_tr){
  #   temp_rise <- temp_rise[1:n_years]
  # }
  # else
  #   n_years<-len_tr
  
  # if(timeframe>(n_years-500))
  #   stop('timeframe is longer than the number of years of temperature data provided; check temp_rise or n_years')
  
  zeta<-0.0006       # zeta: Transient Climate Response to cumulative Emissions TCRE 	0.0006	°C/GtCO2
  kappa<-0.0077      #  % output loss for 1°C (9xlarger for 3°C)
  # Total damages=exp(-kappa*T²)=exp(-kappa*zeta²*S²)
  gamma <- kappa * 2 # Total damages=exp(-gamma/2*T²)=exp(-gamma/2*zeta²S)
  
  years<-start_year:(start_year+n_years-1)
  time<-years-start_year
  
  gdp_world <- GDP_start * (1+r_GDP_growth)^time
  marginal_damage<-gamma * temp_rise^marginal_damage_exponent * zeta * gdp_world
  discount<-exp(-r_discount*time)
  marginal_scc<-marginal_damage * discount
  
  end<-timeframe
  total_scc<-sapply(1:timeframe, function(i){
    sum(marginal_scc[i:timeframe])/discount[i]
  }) %>% unlist()
  
  #   total_scc<-sapply(0:end, function(i){
  #   timeframe_i<-ifelse((i+timeframe)<=end, i+timeframe, end)
  #   sum(marginal_scc[(i+1):timeframe_i])/discount[i+1]
  # }) %>% unlist()
  
  years<-years[1:length(total_scc)]
  names(total_scc)<-as.character(years)
  return(total_scc)
}

scc<-scc_gen(temp_rise, timeframe = 1500, marginal_damage_exponent = 1)


#Initialize a list to store results
scc_list <- list()

# Loop through each discount rate
for (r_discount in discount_rates) {
  scc <- scc_gen(temp_rise, r_discount = r_discount, timeframe = 1500, marginal_damage_exponent = marginal_damage_exponent)
  scc_discounted <- scc[1:n_years] / (1 + r_discount)^(0:(n_years - 1))
  
  # Create a data frame for the current discount rate
  df_current <- as.data.frame(scc_discounted) %>%
    mutate(year = row_number(), discount_rate = paste0(r_discount * 100, "%"))
  
  # Append to the list
  scc_list[[length(scc_list) + 1]] <- df_current
}

# Combine all data frames into one
df_scc <- bind_rows(scc_list)

# Plot the results
ggplot(df_scc, aes(x = year, y = scc_discounted, color = discount_rate)) +
  geom_line() +
  labs(title = "Discounted SCC over time for different discount rates",
       x = "Year",
       y = "SCC Discounted",
       color = "Discount Rate") +
  theme_minimal()

# Write CSV
write.csv(df_scc, "Outputs/SocialDiscountRates_2_4_6pc_DR.csv", row.names = FALSE)

