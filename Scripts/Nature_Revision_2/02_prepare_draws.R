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

#params
years <- 0:75

# ============================================================
# OUTPUT LOCATIONS (simple, deterministic)
# ============================================================

source(file.path("Scripts", "Nature_Revision_2", "_config.R"))
paths <- nr2_step_paths("02_draws")
nr2_ensure_dirs(paths)

# -----------------------------------------
# Inputs: outputs of step 01
# -----------------------------------------
joint_model_path <- file.path(nr2_out_root, "01_one_model", "models", "unified_linear_carbon_model.rds")
stopifnot(file.exists(joint_model_path))
joint_model <- readRDS(joint_model_path)

# -----------------------------------------
# Outputs
# -----------------------------------------
prep_out_root <- paths$root
prep_fig_dir <- paths$figures
prep_tab_dir <- paths$tables
prep_rds_dir <- paths$rds

# -----------------------------------------
# Generate draws_long from fitted model output
# -----------------------------------------
states <- levels(joint_model$data$state)
newdata <- expand.grid(
  time = years,
  state = states
) %>%
  dplyr::filter(!(state == "primary" & time != 0))

draws_long <- joint_model %>%
  tidybayes::add_epred_draws(
    newdata = newdata,
    re_formula = NA,
    ndraws = 500
  ) %>%
  transmute(
    draw = .draw,
    habitat = dplyr::case_when(
      state == "plantation_EP" ~ "eucalyptus_current",
      state == "plantation_AF" ~ "albizia_current",
      TRUE ~ as.character(state)
    ),
    functionalhabAge = time,
    ACD = pmax(.epred, 0)
  )

draws_long %>% group_by(draw, habitat) %>% count()

# get slopes for once-logged forest recovery for each draw - COME BACK TO AMD CHECK
delta <- 1  # 1-year slope
slopes_once_logged <- draws_long %>%
  filter(habitat == "once_logged") %>%
  arrange(functionalhabAge, draw) %>%
  group_by(draw) %>%
  summarise(
    slope_once_logged = (ACD[functionalhabAge == delta] - ACD[functionalhabAge == 0]) / delta,
    .groups = "drop"
  ) 
# ------------------------------------------------
# enforce plateau of ACD to primary mean and CIs
# ------------------------------------------------
# for each posterior draw, if the ACD estimate surpasses the primary estimate, replace with the primary estimate to enforce plateau
# Convert to data.table
dt <- as.data.table(draws_long)
primary_dt <- dt[habitat == "primary", .(draw, primary_ACD = ACD)]

# Merge primary draws on draw
dt <- merge(dt, primary_dt, by = "draw", all.x = TRUE)

# Apply plateau per draw
dt[ACD > primary_ACD, ACD := primary_ACD]

draws_long <- as_tibble(dt) %>%
  select(-primary_ACD)

## ============================================================
## CUSTOM: build `restored_start` functional habitat per draw
## ============================================================
##
## Concept:
## - Cohort A (28%): 28% per ha is subjected to liana cutting and planting and so follows restored ACD for ages 0–30, then is harvested and
##   resets to once-logged recovery starting at age 0 for true_year 31–75.
## - Cohort B (72%): follows once-logged ACD for ages 0–75.
## - `restored_start` ACD = 0.28 * cohortA + 0.72 * cohortB (per draw, per year).
##
## NOTE on the harvest "reset":
## - We map once-logged age 0..44 onto true_year 31..75 (age0 at year31).

weight_restored_cohort <- 0.285714286 #*percentage harvestable near roads (200/700) - the area of restored forest harvested
weight_once_cohort <- 1-weight_restored_cohort # the area of un-enrichment planted forest that is left recovering as 1L
age_restored_end <- 30
reset_year <- 31

hab_vals <- sort(unique(draws_long$habitat))
once_name <- if ("once_logged" %in% hab_vals) "once_logged" else if ("once-logged" %in% hab_vals) "once-logged" else NA_character_
rest_name <- if ("restored" %in% hab_vals) "restored" else NA_character_

if (is.na(once_name)) stop("Could not find once-logged habitat in draws_long$habitat. Found: ", paste(hab_vals, collapse = ", "))
if (is.na(rest_name)) stop("Could not find restored habitat in draws_long$habitat. Found: ", paste(hab_vals, collapse = ", "))

rest_dt <- draws_long %>%
  filter(habitat == rest_name) %>%
  select(draw, functionalhabAge, ACD_restored = ACD)

once_dt <- draws_long %>%
  filter(habitat == once_name) %>%
  select(draw, functionalhabAge, ACD_once = ACD)

# Cohort A: restored ages 0..30
cohortA_rest <- rest_dt %>%
  filter(functionalhabAge %in% 0:age_restored_end) %>%
  transmute(draw, functionalhabAge, ACD_A = ACD_restored)

# Cohort A continuation: once-logged ages 0..44 mapped to years 31..75
cohortA_once <- once_dt %>%
  filter(functionalhabAge %in% 0:(max(years) - reset_year)) %>%
  transmute(draw, functionalhabAge = functionalhabAge + reset_year, ACD_A = ACD_once)

cohortA <- bind_rows(cohortA_rest, cohortA_once) %>%
  filter(functionalhabAge %in% years)

# Cohort B: once-logged ages 0..75
cohortB <- once_dt %>%
  filter(functionalhabAge %in% years) %>%
  transmute(draw, functionalhabAge, ACD_B = ACD_once)

# Combine cohorts
restored_start <- cohortA %>%
  left_join(cohortB, by = c("draw", "functionalhabAge")) %>%
  mutate(
    ACD = weight_restored_cohort * ACD_A + weight_once_cohort * ACD_B,
    habitat = "restored_start"
  ) %>%
  select(draw, habitat, functionalhabAge, ACD)

# Enforce the same per-draw primary plateau for restored_start
primary_by_draw <- draws_long %>%
  filter(habitat == "primary") %>%
  group_by(draw) %>%
  summarise(primary_ACD = mean(ACD, na.rm = TRUE), .groups = "drop")

restored_start <- restored_start %>%
  left_join(primary_by_draw, by = "draw") %>%
  mutate(ACD = pmin(ACD, primary_ACD)) %>%
  select(-primary_ACD)

# Diagnostic plot (one draw): cohort A vs cohort B vs combined
draw_plot_id <- sort(unique(draws_long$draw))[1]
plot_df <- cohortA %>%
  filter(draw == draw_plot_id) %>%
  left_join(cohortB %>% filter(draw == draw_plot_id), by = c("draw", "functionalhabAge")) %>%
  left_join(restored_start %>% filter(draw == draw_plot_id) %>% select(draw, functionalhabAge, ACD_comb = ACD),
            by = c("draw", "functionalhabAge")) %>%
  pivot_longer(
    cols = c(ACD_A, ACD_B, ACD_comb),
    names_to = "series",
    values_to = "ACD"
  ) %>%
  mutate(series = recode(
    series,
    ACD_A = paste0("Cohort A: restored 0–", age_restored_end, " then once-logged"),
    ACD_B = "Cohort B: once-logged 0–75",
    ACD_comb = paste0("Combined: 0.28*A + 0.72*B (restored_start)")
  ))

p_restored_start <- ggplot(plot_df, aes(x = functionalhabAge, y = ACD, colour = series)) +
  geom_line(linewidth = 1.0, alpha = 0.9, na.rm = TRUE) +
  geom_vline(xintercept = age_restored_end, linetype = 2, colour = "grey40") +
  geom_vline(xintercept = reset_year, linetype = 3, colour = "grey40") +
  theme_minimal(base_size = 13) +
  labs(
    title = paste0("restored_start construction (one posterior draw: ", draw_plot_id, ")"),
    subtitle = paste0("Cohort A uses restored 0–", age_restored_end, ", then resets to once-logged at year ", reset_year,
                      "; Combined = 0.28*A + 0.72*B"),
    x = "functionalhabAge (years)",
    y = "ACD (Mg/ha)",
    colour = NULL
  )

ggsave(
  filename = file.path(prep_fig_dir, "restored_start__cohorts_and_combined__one_draw.png"),
  plot = p_restored_start,
  width = 9.5,
  height = 5.5,
  dpi = 220
)




# ------------------------------------------------
# Summary  plot to visualise predictions
# ------------------------------------------------
ACD_summary <- draws_long %>%
  group_by(habitat, functionalhabAge) %>%
  summarise(
    mean = mean(ACD),
    lwr  = quantile(ACD, 0.025),
    upr  = quantile(ACD, 0.975),
    .groups = "drop"
  )

#might need to go back and look at albizia...
quick_plot <- ggplot(ACD_summary, aes(x = functionalhabAge, y = mean,
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
quick_plot

#______________________________________
#INFER twice-logged dynamics ####
#______________________________________

#Longhand for what is going on below 
#define amount of ACD loss from second rotation (ACD loss from getting 31.2m3 of timber in second harvest)

#RATIONALE - BEING USED for 2L -> 2L: assume forest that starts out as twice logged was logged fifteen years 
# #after first logging (e.g. first harvest = t= -15, second harvest t = 0; tracks reality)
# #1. For model of ACD recovery after once-logging,  it looks like the first harvesting rotation led to a decline of 148.8714 mg/ha   
# ACD_l1_yr0 <-  ACD_summary %>% filter(functionalhabAge == 0 &habitat == "once-logged") %>% select(mean) %>% pull()
# primary_ACD <- ACD_summary %>%  filter(habitat == "primary") %>% select(mean) %>% pull()
# decline_from_first_logging <- primary_ACD - l1_yr0
# #2 Thus 148.8714/112.96 =  1.317912 decline in ACD per m3 harvest
# carbon_loss_pr_m3 <- decline_from_first_logging/112.96 
# #3. Again, from models, once-logged forest has an ACD of 75.08815  after 15 years of recovery.  
# ACD_l1_yr15 <- ACD_summary %>% filter(functionalhabAge == 15 &habitat == "once-logged") %>% select(mean) %>% pull()
# ACD_l1_yr30 <- ACD_summary %>% filter(functionalhabAge == 30 &habitat == "once-logged") %>% select(mean) %>% pull()
# #4. Twice-logging removes 31.24 m3 of timber (+/- 10.4) 
# # If I assume the same ACD loss per m3 harvested as once-logged 
# ACD_loss_from_2nd_harvest <-  31.24 * carbon_loss_pr_m3 #= 41.17158 = ACD loss from second harvest
# #7. Thus given a 15yr once-logged ACD before logging, then:
# ACD_2L_starting2L_15yrAfter1L <-  ACD_l1_yr15 - ACD_loss_from_2nd_harvest #= 40.62932 in twice-logged in yr 0   
# ACD_2L_yr0_30yrAfter1L_30yrAfter1L <- ACD_l1_yr30 - ACD_loss_from_2nd_harvest #= 40.62932 in twice-logged in yr 0   

# --- constants ---
vol_first_rotation <- 112.96   # m³ removed in first harvest
vol_second_rotation <- 31.24   # m³ removed in second harvest

# --- extract per-draw ACD values ---
once_draws <- draws_long %>% filter(habitat == "once_logged")
primary_draws <- draws_long %>%
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

#across posterior draws, quite substantial variation in loss per m3
plot(loss_per_m3_df$loss_per_m3)
mean(loss_per_m3_df$loss_per_m3)

# --- define 2L alternative slopes as factors of once-logged slope ---
slope_factors <- c(0.8, 1, 1.2)  # slower = 80%, same = 100%, faster = 120% than once-logged recovery


# second-rotation ACD loss per draw
loss_second_df <- loss_per_m3_df %>%
  mutate(
    loss_second = vol_second_rotation * loss_per_m3
  ) %>%
  select(draw, loss_second)


twice_starting_ACD <- once_age15 %>%
  left_join(once_age30, by = "draw") %>%
  left_join(loss_second_df, by = "draw") %>%
  mutate(
    ACD_2L_start_15yrAfter1L_raw = once_ACD_age15 - loss_second,
    ACD_2L_start_30yrAfter1L_raw = once_ACD_age30 - loss_second,
    
    # clamp to minimum plausible ACD
    ACD_2L_start_15yrAfter1L = pmax(ACD_2L_start_15yrAfter1L_raw, 15),
    ACD_2L_start_30yrAfter1L = pmax(ACD_2L_start_30yrAfter1L_raw, 15)
  ) %>%
  select(draw,
         ACD_2L_start_15yrAfter1L,
         ACD_2L_start_30yrAfter1L)

plot(twice_starting_ACD$ACD_2L_start_15yrAfter1L)
plot(twice_starting_ACD$ACD_2L_start_30yrAfter1L)

# --- expand twice-logged draws across slope scenarios]

ACD_twice_draws <- twice_starting_ACD %>%
  left_join(slopes_once_logged, by = "draw") %>%
  tidyr::expand_grid(
    functionalhabAge = years,
    slope_factor = slope_factors
  ) %>% unique() %>% 
  # apply scaled slope to each draw
  mutate(
    slope_scaled = slope_once_logged * slope_factor,
    # raw ACD trajectories for different slopes
    ACD_twice_logged_15yrStart_raw = ACD_2L_start_15yrAfter1L + (slope_scaled * functionalhabAge),
    ACD_twice_logged_30yrStart_raw = ACD_2L_start_30yrAfter1L + (slope_scaled * functionalhabAge)
  ) %>%
  select(draw, functionalhabAge, slope_factor,
         ACD_twice_logged_15yrStart_raw, ACD_twice_logged_30yrStart_raw)


#ensure correct plateau on a per-draw basis
dt <- as.data.table(ACD_twice_draws)
dt <- merge(dt, primary_draws, by = "draw", all.x = TRUE)
dt[, ACD_twice_logged_15yrStart := pmin(ACD_twice_logged_15yrStart_raw, primary_ACD)]
dt[, ACD_twice_logged_30yrStart := pmin(ACD_twice_logged_30yrStart_raw, primary_ACD)]
dt[, primary_ACD := NULL]
ACD_twice_draws <- as_tibble(dt)

#PLOT

# --- reshape for plotting ---
ACD_twice_long <- ACD_twice_draws %>%
  pivot_longer(
    cols = c(ACD_twice_logged_15yrStart, ACD_twice_logged_30yrStart),
    names_to = "start_age",
    values_to = "ACD"
  ) %>%
  mutate(
    start_age = case_when(
      start_age == "ACD_twice_logged_15yrStart" ~ "15yrAfter1L - e.g if parcel starts scenario 2L",
      start_age == "ACD_twice_logged_30yrStart" ~ "30yrAfter1L - e.g if primary goes to 2L during scenario"
    ),
    slope_factor = factor(slope_factor, levels = slope_factors)
  )  %>%  
  select(-c(ACD_twice_logged_15yrStart_raw, ACD_twice_logged_30yrStart_raw))

# --- summarize posterior draws ---
ACD_twice_summary <- ACD_twice_long %>%
  group_by(functionalhabAge, slope_factor, start_age) %>%
  summarise(
    mean_ACD = mean(ACD),
    lwr_ACD  = quantile(ACD, 0.025),
    upr_ACD  = quantile(ACD, 0.975),
    .groups = "drop"
  )

#quick comparison of 0-year intercepts - nb 2L mean intercept for 15y after once-logged is moderately higher than after once-logging. I think this makes sense - removing much less wood/biomass
print(ACD_summary)
print(ACD_twice_summary)

# --- plot ---
twice_logged_plot <-ACD_twice_summary %>% 
  # filter(slope_factor == 1.2) %>% 
  ggplot(aes(x = functionalhabAge, y = mean_ACD,
             color = slope_factor, fill = slope_factor)) +
  geom_ribbon(aes(ymin = lwr_ACD, ymax = upr_ACD), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  facet_grid(start_age ~ slope_factor) +
  
  # facet_wrap(~start_age, ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Years Since 2nd Logging",
    y = "Aboveground Carbon Density (ACD)",
    color = "Slope factor",
    fill = "Slope factor",
    title = "Twice-logged recovery trajectories with different slopes"
  )

cowplot::plot_grid(quick_plot, twice_logged_plot)

#_______________________________________________________________________________
#transition rules 
#_______________________________________________________________________________
#primary - no rules 

#once logged
#a. P -> 1L  If parcel goes to from primary to once-logged, we start at 1L yr 0
#b. 1L -> 1L if parcel starts as once-logged and stays 1L, assume logging happened 15 yrs befor scnanario start (once-logged_start). Likewise, restoration happens immediately after once-logging at t-15
#c. 1L -> 2L if parcel starts as once logged before re-harvest, we assume 1L recovery (already explicit in scenarios) that resets to yr 0 once-logged values.

#for b. 
# oncelogged <- L1_R %>%
#   filter(original_habitat == "once-logged") %>% 
#   #once-logged forest starting in our scenarios is actually 15 yrs old already 
#   #so this corrects for this. 
#   mutate(functionalhabAge = functionalhabAge - 15) %>% 
#   #filter delay to always be >15 yrs for any forest parcel starting with once-logged habitat
#   filter(harvest_delay>15) 

# 
# #if parcel starts as twice-logged, use the 15yr after 1L curves - with different slope factors, tracking the assumption that all 2nd harvest occur at yr 0 for forest that was harvested for the first time 15 yrs prior 
# #if parcel starts as primary and goes to twice-logged, use the 30yrafter 1L, with different slope factors 
# unique(ACD_draws$habitat)
# ACD_twice_long <- ACD_twice_long %>% select(-ACD) %>%   pivot_longer(
#     cols = c(ACD_twice_logged_15yrStart_raw, ACD_twice_logged_30yrStart_raw),
#     names_to = "hab_trans_rules",
#     values_to = "ACD"
#   ) %>%
#   mutate(
#     hab_trans_rules = case_when(
#       hab_trans_rules == "ACD_twice_logged_15yrStart_raw" ~ "15yrAfter1L",
#       hab_trans_rules == "ACD_twice_logged_30yrStart_raw" ~ "30yrAfter1L",
#       TRUE ~ hab_trans_rules
#     )
#   ) %>%
#   mutate(habitat = "twice-logged") %>%
#   select(-hab_trans_rules)
# 
# ACD_draws
draws_long <- draws_long %>%
  mutate(
    slope_factor = 1,
    start_age = NA
  )


#make sure habitat starting as twice-logged has a distrinct trajectory
ACD_twice_long <- ACD_twice_long %>%  
  mutate(habitat = "twice-logged" ) %>% 
  mutate(
    habitat = case_when(
      str_starts(start_age, "15yrAfter") ~ "twice_logged_start",
      TRUE ~ habitat
    )
  )


# combine once-logged, primary, twice-logged and restored and plantation
draws_long_bind <- draws_long %>%
  mutate(slope_factor = factor(slope_factor, levels = slope_factors)) %>%
  select(draw, habitat, functionalhabAge, ACD, slope_factor, start_age)

ACD_twice_bind <- ACD_twice_long %>%
  mutate(slope_factor = factor(as.character(slope_factor), levels = as.character(slope_factors))) %>%
  select(draw, habitat, functionalhabAge, ACD, slope_factor, start_age)

final_ACD_draws <- bind_rows(draws_long_bind, ACD_twice_bind)

head(final_ACD_draws)

#Build a plot of above ground carbon dynamics through times
unique(final_ACD_draws$habitat)
head(final_ACD_draws)
final_ACD_summary <- final_ACD_draws %>% group_by(habitat, slope_factor,functionalhabAge, start_age) %>% 
  summarise(
    mean_ACD = mean(ACD),
    lwr_ACD  = quantile(ACD, 0.025),
    upr_ACD  = quantile(ACD, 0.975),
    .groups = "drop"
  ) %>%  
  filter(habitat != "twice_logged_start")

#make carbon through time plot #####
p_1_r <- final_ACD_summary %>% 
  filter(habitat %in% c("primary", "once_logged", "restored"))
unique(p_1_r$habitat)

# Get full age range used in other habitats
age_range <- range(p_1_r$functionalhabAge, na.rm = TRUE)

primary_expanded <- p_1_r %>%
  filter(habitat == "primary") %>%
  select(-functionalhabAge) %>%   # remove the existing column
  crossing(functionalhabAge = seq(age_range[1], age_range[2]))

p_1_r <- bind_rows(
  p_1_r %>% filter(habitat != "primary"),
  primary_expanded
)


tl <-final_ACD_summary %>% filter(habitat == "twice-logged")
pl <- final_ACD_summary %>% filter(habitat%in% c("eucalyptus_current", "albizia_current"))


# p_1_plot %>% group_by(habitat) %>%   ggplot(aes(x = functionalhabAge, y = mean_ACD,
#            color = slope_factor, fill = habitat)) +
#   geom_ribbon(aes(ymin = lwr_ACD, ymax = upr_ACD), alpha = 0.2, color = NA) +
#   geom_line(linewidth = 1.2) +
#   facet_grid(start_age ~ slope_factor) +
#   
#   # facet_wrap(~start_age, ncol = 1, scales = "free_y") +
#   theme_minimal(base_size = 14) +
#   labs(
#     x = "Years Since Logging",
#     y = "Aboveground Carbon Density (ACD)",
#     color = "Habitat",
#     title = "Primary, Once-logged and Restored trajectories"
#   ) + 
#   xlim(0,60)

p_1_r_plot <- p_1_r %>% 
  ggplot(aes(x = functionalhabAge,
             y = mean_ACD,
             colour = habitat,
             linetype = habitat,
             group = habitat)) +
  
  geom_ribbon(aes(ymin = lwr_ACD, ymax = upr_ACD),
              fill = "grey70",
              alpha = 0.25,
              colour = NA) +
  
  geom_line(linewidth = 1.1) +
  
  scale_colour_grey(start = 0.2, end = 0.6) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank()
    
  ) +
  
  labs(
    x = "Years Since Logging",
    y = "Aboveground Carbon Density (ACD)",
    colour = "Habitat",
    linetype = "Habitat",
    title = "Primary, Once-logged and Restored trajectories"
  ) +
  
  coord_cartesian(xlim = c(0,60))


pl_plot <- pl %>%  
  ggplot(aes(x = functionalhabAge,
             y = mean_ACD,
             colour = habitat,
             linetype = habitat,
             group = habitat)) +
  
  geom_ribbon(aes(ymin = lwr_ACD, ymax = upr_ACD),
              fill = "grey70",
              alpha = 0.25,
              colour = NA) +
  
  geom_line(linewidth = 1.1) +
  
  scale_colour_grey(start = 0.2, end = 0.6) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "top", 
    legend.title = element_blank()
  ) +
  
  labs(
    x = "Years Since Logging",
    y = "Aboveground Carbon Density (ACD)",
    colour = "Habitat",
    linetype = "Habitat",
    title = "Plantation trajectories"
  ) +
  xlim(0,12)+
  coord_cartesian(xlim = c(0,12))+
  ylim(0,250)

tl_plot <- tl %>%  
  ggplot(aes(x = functionalhabAge,
             y = mean_ACD,
             colour = slope_factor,
             linetype = habitat,
             group = habitat)) +
  
  geom_ribbon(aes(ymin = lwr_ACD, ymax = upr_ACD),
              fill = "grey70",
              alpha = 0.25,
              colour = NA) +
  
  geom_line(linewidth = 1.1) +
  
  scale_colour_grey(start = 0.2, end = 0.6) +
  
  theme_minimal(base_size = 14) +
  facet_wrap(~slope_factor)+
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  
  labs(
    x = "Years Since Logging",
    y = "Aboveground Carbon Density (ACD)",
    colour = "Habitat",
    linetype = "Habitat",
    title = "Interpolated twice-logged trajectories
    with different slope assumptions"
  ) +
  xlim(0,60)+
  coord_cartesian(xlim = c(0,60))

log_plant <- cowplot::plot_grid(p_1_r_plot,pl_plot )
final_carbon_curve_fig <- cowplot::plot_grid(log_plant, tl_plot, ncol = 1)
ggsave(
  filename = file.path(prep_fig_dir, "all_carbon_curves.png"),
  plot = final_carbon_curve_fig,
  units = "mm",
  height = 297,
  width = 220
)


# add *_start values to ACD long

# Append restored_start into the final draws table (so downstream code treats it like any habitat)
# IMPORTANT: give it slope_factor = 1 so it is retained when later scripts
# filter to a single slope_factor trajectory.
restored_start <- restored_start %>%
  mutate(
    slope_factor = factor("1", levels = as.character(slope_factors)),
    start_age = NA
  ) %>%
  select(draw, habitat, functionalhabAge, ACD, slope_factor, start_age)


#all restored start was planted at the same time as the once-logging (ie at t-15)
restored_start <- restored_start %>%  
  mutate(functionalhabAge = functionalhabAge - 15) %>%  
  filter(functionalhabAge >-1)
final_ACD_draws <- bind_rows(final_ACD_draws, restored_start)


# Build once-logged_start from once-logged predictions:
# same ACD draws, but shifted back by 15 years so scenario year 0 maps to
# once-logged functional age 15.
#since scenarios that start once-logged were harvested 15 years prior
once_logged_start <- draws_long_bind %>%
  filter(habitat %in% c("once_logged", "once-logged")) %>%
  mutate(
    habitat = "once-logged_start",
    functionalhabAge = functionalhabAge - 15
  ) %>%
  filter(functionalhabAge > -1) %>%
  select(draw, habitat, functionalhabAge, ACD, slope_factor, start_age)

final_ACD_draws <- bind_rows(final_ACD_draws, once_logged_start)

#_________________________
#Outputs ####
#_________________________

# save ACD draw (no below ground, and assuming different slope trajectories for twice logged forest)
saveRDS(final_ACD_draws, file.path(prep_rds_dir, "acdraws_aboveground.rds"))


