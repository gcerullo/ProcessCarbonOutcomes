# ----------------------------------------------------------------------------
# Nature Revision 2 — joint habitat carbon model (fits)
#
# This is where I fit the shared state-space / recovery model across habitats and push posterior summaries into the NR2 folder structure.
# Inputs: cleaned carbon plot/habitat RDS and covariate tables referenced in the CONFIG section; sources _config.R.
# Outputs: models/ and rds/ under the active NR2 step root (posterior draws, summaries for step 02).
# ----------------------------------------------------------------------------

# Nature Revision 2 version of: 01_02_one_model_to_rule_them_all.R
# Purpose:
# - Data prep → prior predictive checks → fit model → posterior predictive checks → simple trajectory plots
# - Weakly informative priors calibrated from ecological assumptions
# - Write all outputs under Outputs/Nature_Revision_Outputs/NR2/current/ (no CursorOutputs/CursorModels)

# ============================================================
# 0) LIBRARIES + OUTPUT LOCATIONS (Nature Revision Outputs)
# ============================================================

library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)

set.seed(123)

# ============================================================
# 0b) RUN CONTROLS (so you can iterate on priors fast)
# ============================================================

# Suggested workflow:
# - Set RUN_JOINT_MODEL = FALSE while tuning priors
# - Iterate until prior predictive plots/diagnostics look sensible
# - Then set RUN_JOINT_MODEL = TRUE and RUN_EXPORT_DRAWS = TRUE
# RUN_PRIOR_ONY <- TRUE
# RUN_JOINT_MODEL <- TRUE
# RUN_EXPORT_DRAWS <- FALSE
# RUN_PLOT_TRAJECTORIES <- TRUE

# Faster defaults for prior-only runs (enough for stable plots)
# CHAINS_PRIOR <- 2
# ITER_PRIOR <- 1500
# WARMUP_PRIOR <- 500
# CORES_PRIOR <- 2

# More expensive settings for the actual fitted model
CHAINS_FIT <- 4
ITER_FIT <- 6000
WARMUP_FIT <- 2000
CORES_FIT <- 4

# ============================================================
# 0c) OUTPUT LOCATIONS (simple, deterministic)
# ============================================================

source(file.path("Scripts", "Nature_Revision_2", "_config.R"))
paths <- nr2_step_paths("01_one_model", include_models = TRUE)
nr2_ensure_dirs(paths)

step_root <- paths$root
nr_model_root <- paths$models
nr_fig_dir <- paths$figures
nr_tab_dir <- paths$tables
nr_rds_dir <- paths$rds

writeLines(capture.output(sessionInfo()), con = file.path(step_root, "sessionInfo.txt"))

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
    time = 0,
    plot_id = Plot
  ) %>%
  select(ACD, time, state, plot_id)

# ============================================================
# 3) COMBINE + FINAL PREP
# ============================================================

all_data <- bind_rows(plantation, Logged, Primary)
all_data$state <- relevel(factor(all_data$state), ref = "primary")
all_data$plot_id <- factor(all_data$plot_id)

# Plot identifiers come from different studies; keep them separate per state
all_data <- all_data %>%
  mutate(plot_state = factor(interaction(state, plot_id, drop = TRUE)))

# Primary has no temporal variation; keep a helper column to make intent explicit
all_data <- all_data %>%
  mutate(time_state = ifelse(state == "primary", 0, time))

# Prediction grid for prior predictive trajectories
prediction_grid <- expand.grid(
  time = 0:75,
  state = levels(all_data$state)
) %>%
  dplyr::filter(!(state == "primary" & time != 0))

# ============================================================
# 4) PRIOR CALIBRATION HELPERS (turn assumptions into priors)
# ============================================================

sd_from_95 <- function(lower, upper) {
  (upper - lower) / (2 * qnorm(0.975))
}

normal_from_95 <- function(mean, lower, upper) {
  list(mean = mean, sd = sd_from_95(lower, upper))
}

# We use a truncated Gaussian likelihood (lb = 0) to enforce a hard zero on ACD,
# while keeping a linear-in-time mean structure on the original scale.
assumptions <- list(
  # Absolute ACD at time = 0 by state
  primary_acd_t0 = list(mean = 200, lower95 = 100, upper95 = 300),
  once_logged_acd_t0 = list(mean = 70, lower95 = 20, upper95 = 140),
  restored_acd_t0 = list(mean = 80, lower95 = 30, upper95 = 160),
  plantation_AF_acd_t0 = list(mean = 15, lower95 = 0, upper95 = 60),
  plantation_EP_acd_t0 = list(mean = 15, lower95 = 0, upper95 = 60),

  # Slopes are on the original ACD scale (Mg C ha^-1 yr^-1)
  primary_slope = list(mean = 0.0, lower95 = -0.5, upper95 = 0.5),
  once_logged_slope = list(mean = 2.0, lower95 = 0.0, upper95 = 5.0),
  restored_slope = list(mean = 2.5, lower95 = 0.5, upper95 = 6.0),
  plantation_AF_slope = list(mean = 6.0, lower95 = 0.0, upper95 = 12.0),
  plantation_EP_slope = list(mean = 8.0, lower95 = 0.0, upper95 = 16.0),

  # variability terms
  plot_sd = list(df = 3, scale = 20),
  # sigma model is on log(sigma) scale in brms (sigma has log link)
  sigma_log = list(mean = log(30), sd = 0.7)
)

pri_primary_intercept <- normal_from_95(
  mean = assumptions$primary_acd_t0$mean,
  lower = assumptions$primary_acd_t0$lower95,
  upper = assumptions$primary_acd_t0$upper95
)

pri_once_logged_intercept <- normal_from_95(
  mean = assumptions$once_logged_acd_t0$mean,
  lower = assumptions$once_logged_acd_t0$lower95,
  upper = assumptions$once_logged_acd_t0$upper95
)

pri_restored_intercept <- normal_from_95(
  mean = assumptions$restored_acd_t0$mean,
  lower = assumptions$restored_acd_t0$lower95,
  upper = assumptions$restored_acd_t0$upper95
)

pri_plantation_AF_intercept <- normal_from_95(
  mean = assumptions$plantation_AF_acd_t0$mean,
  lower = assumptions$plantation_AF_acd_t0$lower95,
  upper = assumptions$plantation_AF_acd_t0$upper95
)

pri_plantation_EP_intercept <- normal_from_95(
  mean = assumptions$plantation_EP_acd_t0$mean,
  lower = assumptions$plantation_EP_acd_t0$lower95,
  upper = assumptions$plantation_EP_acd_t0$upper95
)

slope_primary <- normal_from_95(
  mean = assumptions$primary_slope$mean,
  lower = assumptions$primary_slope$lower95,
  upper = assumptions$primary_slope$upper95
)

slope_once_logged <- normal_from_95(
  mean = assumptions$once_logged_slope$mean,
  lower = assumptions$once_logged_slope$lower95,
  upper = assumptions$once_logged_slope$upper95
)

slope_restored <- normal_from_95(
  mean = assumptions$restored_slope$mean,
  lower = assumptions$restored_slope$lower95,
  upper = assumptions$restored_slope$upper95
)

slope_plantation_AF <- normal_from_95(
  mean = assumptions$plantation_AF_slope$mean,
  lower = assumptions$plantation_AF_slope$lower95,
  upper = assumptions$plantation_AF_slope$upper95
)

slope_plantation_EP <- normal_from_95(
  mean = assumptions$plantation_EP_slope$mean,
  lower = assumptions$plantation_EP_slope$lower95,
  upper = assumptions$plantation_EP_slope$upper95
)

priors_calibrated <- c(
  # Absolute intercepts at time = 0 (because model uses `0 + state`)
  set_prior(paste0("normal(", pri_primary_intercept$mean, ", ", pri_primary_intercept$sd, ")"),
            class = "b", coef = "stateprimary"),
  set_prior(paste0("normal(", pri_once_logged_intercept$mean, ", ", pri_once_logged_intercept$sd, ")"),
            class = "b", coef = "stateonce_logged"),
  set_prior(paste0("normal(", pri_restored_intercept$mean, ", ", pri_restored_intercept$sd, ")"),
            class = "b", coef = "staterestored"),
  set_prior(paste0("normal(", pri_plantation_AF_intercept$mean, ", ", pri_plantation_AF_intercept$sd, ")"),
            class = "b", coef = "stateplantation_AF"),
  set_prior(paste0("normal(", pri_plantation_EP_intercept$mean, ", ", pri_plantation_EP_intercept$sd, ")"),
            class = "b", coef = "stateplantation_EP"),

  set_prior(paste0("normal(", slope_primary$mean, ", ", slope_primary$sd, ")"),
            class = "b", coef = "stateprimary:time"),
  set_prior(paste0("normal(", slope_once_logged$mean, ", ", slope_once_logged$sd, ")"),
            class = "b", coef = "stateonce_logged:time"),
  set_prior(paste0("normal(", slope_restored$mean, ", ", slope_restored$sd, ")"),
            class = "b", coef = "staterestored:time"),
  set_prior(paste0("normal(", slope_plantation_AF$mean, ", ", slope_plantation_AF$sd, ")"),
            class = "b", coef = "stateplantation_AF:time"),
  set_prior(paste0("normal(", slope_plantation_EP$mean, ", ", slope_plantation_EP$sd, ")"),
            class = "b", coef = "stateplantation_EP:time"),

  set_prior(paste0("student_t(", assumptions$plot_sd$df, ", 0, ", assumptions$plot_sd$scale, ")"), class = "sd"),

  # state-specific residual SD on log scale (sigma has log link)
  set_prior(paste0("normal(", assumptions$sigma_log$mean, ", ", assumptions$sigma_log$sd, ")"),
            class = "b", dpar = "sigma")
)

write.csv(
  as.data.frame(priors_calibrated),
  file = file.path(nr_tab_dir, "priors_calibrated.csv"),
  row.names = FALSE
)

# ============================================================
# 5) PRIOR PREDICTIVE CHECKS (calibrated priors)
# ============================================================
# 
# #if (RUN_PRIOR_ONLY) {
#   prior_only_model <- brm(
#     bf(
#       ACD | trunc(lb = 0) ~ 0 +
#         state +
#         0 + state:time +
#         (1 | plot_state),
#       sigma ~ 0 + state,
#       center = FALSE
#     ),
#     data = all_data,
#     family = gaussian(),
#     prior = priors_calibrated,
#     sample_prior = "only",
#     chains = CHAINS_PRIOR,
#     iter = ITER_PRIOR,
#     warmup = WARMUP_PRIOR,
#     cores = CORES_PRIOR,
#     backend = "cmdstanr",
#     seed = 123
#   )
# 
#   saveRDS(prior_only_model, file.path(nr_model_root, "prior_only_model.rds"))
# 
#   # Basic prior predictive distribution
#   yrep <- posterior_predict(prior_only_model, ndraws = 200)
#   keep <- which(rowSums(is.na(yrep)) == 0)
#   yrep_ok <- yrep[keep, , drop = FALSE]
#   if (nrow(yrep_ok) >= 5) {
#     pdf(file.path(nr_fig_dir, "prior_pp_check_density.pdf"), width = 8, height = 6)
#     show_n <- min(50, nrow(yrep_ok))
#     print(bayesplot::ppc_dens_overlay(y = all_data$ACD, yrep = yrep_ok[seq_len(show_n), ]))
#     dev.off()
#   }
#   write.csv(
#     tibble::tibble(
#       ndraws_requested = 200,
#       ndraws_kept = nrow(yrep_ok),
#       ndraws_dropped_due_to_NA = nrow(yrep) - nrow(yrep_ok)
#     ),
#     file = file.path(nr_tab_dir, "prior_predictive_yrep_NA_diagnostics.csv"),
#     row.names = FALSE
#   )
# 
#   prior_draws <- prior_only_model %>%
#     add_epred_draws(
#       newdata = prediction_grid,
#       re_formula = NA,
#       ndraws = 500
#     ) %>%
#     filter(is.finite(.epred))
# 
#   # Trajectories plot
#   p_prior_traj <- ggplot(prior_draws, aes(x = time, y = .epred, group = .draw)) +
#     geom_line(alpha = 0.04) +
#     facet_wrap(~ state) +
#     theme_minimal() +
#     coord_cartesian(ylim = c(-50, 450)) +
#     labs(
#       y = "Simulated ACD (Mg C ha^-1)",
#       x = "Years",
#       title = "Prior predictive trajectories (calibrated priors)"
#     )
# 
#   ggsave(
#     filename = file.path(nr_fig_dir, "prior_predictive_trajectories.png"),
#     plot = p_prior_traj,
#     width = 12,
#     height = 7,
#     units = "in",
#     dpi = 200
#   )
# 
#   # Starting carbon distributions
#   p_prior_t0 <- prior_draws %>%
#     filter(time == 0) %>%
#     ggplot(aes(x = .epred)) +
#     geom_density(fill = "grey80") +
#     facet_wrap(~ state, scales = "free") +
#     theme_minimal() +
#     labs(
#       x = "ACD at time = 0",
#       title = "Prior predictive initial ACD by state"
#     )
# 
#   ggsave(
#     filename = file.path(nr_fig_dir, "prior_predictive_initial_ACD.png"),
#     plot = p_prior_t0,
#     width = 12,
#     height = 7,
#     units = "in",
#     dpi = 200
#   )
# 
#   # Quantitative prior diagnostics (helps detect dominating or nonsensical priors)
#   prior_diag <- prior_draws %>%
#     mutate(flag_negative = .epred < 0, flag_too_high = .epred > 400) %>%
#     group_by(state, time) %>%
#     summarise(
#       mean = mean(.epred, na.rm = TRUE),
#       p02_5 = quantile(.epred, 0.025, na.rm = TRUE),
#       p50 = quantile(.epred, 0.5, na.rm = TRUE),
#       p97_5 = quantile(.epred, 0.975, na.rm = TRUE),
#       pr_negative = mean(flag_negative, na.rm = TRUE),
#       pr_gt400 = mean(flag_too_high, na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     filter(time %in% c(0, 10, 25, 50, 75))
# 
#   write.csv(prior_diag, file.path(nr_tab_dir, "prior_predictive_diagnostics.csv"), row.names = FALSE)
# #}

# ============================================================
# 6) FIT THE MODEL
# ============================================================

#if (RUN_JOINT_MODEL) {
  joint_model <- brm(
    bf(
      ACD | trunc(lb = 0) ~ 0 +
        state +
        0 + state:time +
        (1 | plot_state),
      sigma ~ 0 + state,
      center = FALSE
    ),
    data = all_data,
    family = gaussian(),
    prior = priors_calibrated,
    chains = CHAINS_FIT,
    iter = ITER_FIT,
    warmup = WARMUP_FIT,
    cores = CORES_FIT,
    backend = "cmdstanr",
    control = list(adapt_delta = 0.999, max_treedepth = 12),
    seed = 123
  )

  saveRDS(joint_model, file.path(nr_model_root, "unified_linear_carbon_model.rds"))

  # Posterior predictive checks by state
  states <- levels(all_data$state)
  pdf(file.path(nr_fig_dir, "posterior_pp_checks_by_state.pdf"), width = 8, height = 6)
  for (s in states) {
    state_data <- all_data %>% filter(state == s)

    yrep <- posterior_predict(joint_model, newdata = state_data, ndraws = 200)
    keep <- which(rowSums(is.na(yrep)) == 0)
    if (length(keep) < 5) next
    yrep_ok <- yrep[keep, , drop = FALSE]
    show_n <- min(50, nrow(yrep_ok))
    yrep_show <- yrep_ok[seq_len(show_n), , drop = FALSE]

    p1 <- bayesplot::ppc_dens_overlay(y = state_data$ACD, yrep = yrep_show)
    print(p1 + ggtitle(paste("Density overlay -", s)))

    p2 <- bayesplot::ppc_hist(y = state_data$ACD, yrep = yrep_show)
    print(p2 + ggtitle(paste("Histogram -", s)))

    p3 <- bayesplot::ppc_scatter_avg(y = state_data$ACD, yrep = yrep_show)
    print(p3 + ggtitle(paste("Scatter observed vs predicted -", s)))

    p4 <- bayesplot::ppc_stat(y = state_data$ACD, yrep = yrep_show, stat = "mean")
    print(p4 + ggtitle(paste("Predicted means -", s)))

    p5 <- bayesplot::ppc_stat(y = state_data$ACD, yrep = yrep_show, stat = "sd")
    print(p5 + ggtitle(paste("Predicted SDs -", s)))
  }
  dev.off()
#}

# ============================================================
# 6b) STRAIGHTFORWARD TRAJECTORY PLOTS (posterior epred vs time)
# ============================================================

#if (RUN_JOINT_MODEL && RUN_PLOT_TRAJECTORIES) {
  # For logged/restored, you project to 75y; for plantations you care about 0–12y.
  states <- levels(all_data$state)
  time_max_by_state <- tibble::tibble(
    state = states,
    time_max = dplyr::case_when(
      state %in% c("plantation_AF", "plantation_EP") ~ 12,
      state == "primary" ~ 0,
      TRUE ~ 75
    )
  )

  newdata_pred <- time_max_by_state %>%
    tidyr::uncount(weights = time_max + 1, .remove = FALSE) %>%
    group_by(state) %>%
    mutate(time_plot = 0:dplyr::first(time_max)) %>%
    ungroup() %>%
    mutate(time = ifelse(state == "primary", 0, time_plot)) %>%
    select(state, time, time_plot)

  # Population-level (no random effects) trajectories
  draws_pred <- joint_model %>%
    tidybayes::add_epred_draws(
      newdata = newdata_pred,
      re_formula = NA,
      ndraws = 800
    ) %>%
    rename(draw = .draw, ACD = .epred)

  # Credible intervals (median + 95% CI)
  traj_summary <- draws_pred %>%
    group_by(state, time_plot) %>%
    summarise(
      med = median(ACD),
      lwr = quantile(ACD, 0.025),
      upr = quantile(ACD, 0.975),
      .groups = "drop"
    )

  write.csv(
    traj_summary,
    file.path(nr_tab_dir, "posterior_trajectory_summary.csv"),
    row.names = FALSE
  )

  p_ci <- ggplot(traj_summary, aes(x = time_plot, y = med)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_line(linewidth = 1) +
    facet_wrap(~ state, scales = "free_x") +
    theme_minimal() +
    labs(
      x = "Years",
      y = "Predicted ACD (Mg C ha^-1)",
      title = "Posterior predicted trajectories (median ± 95% credible interval)"
    )

  ggsave(
    filename = file.path(nr_fig_dir, "posterior_pred_trajectories_CI.png"),
    plot = p_ci,
    width = 12,
    height = 7,
    units = "in",
    dpi = 220
  )

  # Spaghetti plot (subset of draws)
  set.seed(123)
  keep_draws <- sample(unique(draws_pred$draw), size = min(120, length(unique(draws_pred$draw))))
  draws_spag <- draws_pred %>% filter(draw %in% keep_draws)

  p_spag <- ggplot(draws_spag, aes(x = time_plot, y = ACD, group = draw)) +
    geom_line(alpha = 0.08) +
    facet_wrap(~ state, scales = "free_x") +
    theme_minimal() +
    labs(
      x = "Years",
      y = "Predicted ACD (Mg C ha^-1)",
      title = "Posterior predicted trajectories (spaghetti; population-level)"
    )

  ggsave(
    filename = file.path(nr_fig_dir, "posterior_pred_trajectories_spaghetti.png"),
    plot = p_spag,
    width = 12,
    height = 7,
    units = "in",
    dpi = 220
  )
#}

# ============================================================
# 7) EXPORT DRAWS FOR DOWNSTREAM SCRIPTS
# ============================================================

#if (RUN_EXPORT_DRAWS) {
#  if (!RUN_JOINT_MODEL) {
#    stop("RUN_EXPORT_DRAWS = TRUE requires RUN_JOINT_MODEL = TRUE.")
#  }

  newdata <- expand.grid(
    time = 0:75,
    state = levels(all_data$state)
  ) %>%
    dplyr::filter(!(state == "primary" & time != 0))

  draws_long <- joint_model %>%
    add_epred_draws(
      newdata = newdata,
      re_formula = NA,
      ndraws = 500
    ) %>%
    rename(
      draw = .draw,
      ACD = .epred
    )

  draws_long <- draws_long %>%
    ungroup() %>%
    select(-c(".row", ".chain", ".iteration")) %>%
    mutate(
      state = case_when(
        state == "plantation_EP" ~ "eucalyptus_current",
        state == "plantation_AF" ~ "albizia_current",
        TRUE ~ state
      )
    ) %>%
    rename(habitat = state)

  saveRDS(draws_long, file.path(nr_rds_dir, "onemodel_ACD_draws.rds"))
  write.csv(draws_long, file.path(nr_rds_dir, "onemodel_ACD_draws.csv"), row.names = FALSE)
#}

