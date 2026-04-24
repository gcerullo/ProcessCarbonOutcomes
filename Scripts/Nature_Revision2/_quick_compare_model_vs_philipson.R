# ----------------------------------------------------------------------------
# Nature Revision 2 — quick model vs Philipson recovery check
#
# I use this as a sanity check when I change priors or recovery curvature and want a fast plot/table against the published Philipson-style expectations.
# Inputs: RDS produced by the main carbon model steps (paths defined at top of file).
# Outputs: comparison figures/tables next to the carbon NR2 figure folders (see ggsave/write paths in-script).
# ----------------------------------------------------------------------------

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(posterior)
})

source(file.path("Scripts", "Nature_Revision_2", "_config.R"))
paths <- nr2_step_paths("carbon_recovery_comparison")
nr2_ensure_dirs(paths)

fig_dir <- paths$figures
tab_dir <- paths$tables

joint_model_path <- file.path(
  nr2_out_root,
  "01_one_model",
  "models",
  "unified_linear_carbon_model.rds"
)
stopifnot(file.exists(joint_model_path))

joint_model <- readRDS(joint_model_path)
draws_df <- as.data.frame(as_draws_df(joint_model))

summarise_draw <- function(x) {
  tibble(
    estimate = mean(x, na.rm = TRUE),
    lwr95 = as.numeric(quantile(x, 0.025, na.rm = TRUE)),
    upr95 = as.numeric(quantile(x, 0.975, na.rm = TRUE))
  )
}

model_tbl <- bind_rows(
  summarise_draw(draws_df[["b_stateprimary:time"]]) %>%
    mutate(habitat = "primary", metric = "annual_slope"),
  summarise_draw(draws_df[["b_stateonce_logged:time"]]) %>%
    mutate(habitat = "once-logged", metric = "annual_slope"),
  summarise_draw(draws_df[["b_staterestored:time"]]) %>%
    mutate(habitat = "restored", metric = "annual_slope"),
  summarise_draw(draws_df[["b_stateprimary"]]) %>%
    mutate(habitat = "primary", metric = "acd_at_35"),
  summarise_draw(draws_df[["b_stateonce_logged"]] + 35 * draws_df[["b_stateonce_logged:time"]]) %>%
    mutate(habitat = "once-logged", metric = "acd_at_35"),
  summarise_draw(draws_df[["b_staterestored"]] + 35 * draws_df[["b_staterestored:time"]]) %>%
    mutate(habitat = "restored", metric = "acd_at_35")
) %>%
  mutate(source = "my_model")

# Philipson-equivalent refit from the same raw data and formula family
# used in the original workflow (no need for a saved Philipson model object).
data_phil <- read.csv("RawData/Philipson20_PlotData2.csv")
Logged <- subset(data_phil, Forest == "Logged")
UnLogged <- subset(data_phil, Forest == "UnLogged")

Logged$FACE <- factor(Logged$FACE)
if ("Baseline" %in% levels(Logged$FACE)) {
  Logged$FACE <- relevel(Logged$FACE, ref = "Baseline")
}

m1 <- lmer(
  ACD ~ YearsSinceLogging * FACE + (1 | Plot) + (1 | LoggingMethod:Coupe),
  REML = TRUE,
  data = Logged,
  na.action = na.fail
)

fixed_form <- delete.response(terms(lme4::nobars(formula(m1))))
beta <- fixef(m1)
vc <- as.matrix(vcov(m1))

lincomb_est_ci <- function(x_vec) {
  if (is.matrix(x_vec)) {
    x_named <- as.numeric(x_vec[1, ])
    names(x_named) <- colnames(x_vec)
  } else {
    x_named <- as.numeric(x_vec)
  }
  v <- rep(0, length(beta))
  names(v) <- names(beta)
  common_terms <- intersect(names(v), names(x_named))
  v[common_terms] <- x_named[common_terms]
  est <- sum(v * beta)
  se <- sqrt(as.numeric(t(v) %*% vc %*% v))
  tibble(
    estimate = est,
    lwr95 = est - 1.96 * se,
    upr95 = est + 1.96 * se
  )
}

new_once_34 <- data.frame(
  YearsSinceLogging = 34,
  FACE = factor("Baseline", levels = levels(Logged$FACE))
)
new_once_35 <- data.frame(
  YearsSinceLogging = 35,
  FACE = factor("Baseline", levels = levels(Logged$FACE))
)
new_rest_34 <- data.frame(
  YearsSinceLogging = 34,
  FACE = factor("ProjectScenario", levels = levels(Logged$FACE))
)
new_rest_35 <- data.frame(
  YearsSinceLogging = 35,
  FACE = factor("ProjectScenario", levels = levels(Logged$FACE))
)

x_once_34 <- model.matrix(fixed_form, new_once_34)
x_once_35 <- model.matrix(fixed_form, new_once_35)
x_rest_34 <- model.matrix(fixed_form, new_rest_34)
x_rest_35 <- model.matrix(fixed_form, new_rest_35)

prime_mod <- lmer(ACD ~ 1 + (1 | MeasureTime), data = UnLogged)
prime_est <- as.numeric(fixef(prime_mod)[1])
prime_se <- sqrt(as.numeric(vcov(prime_mod)[1, 1]))

philipson_refit_tbl <- bind_rows(
  tibble(
    habitat = "primary",
    metric = "annual_slope",
    estimate = 0,
    lwr95 = NA_real_,
    upr95 = NA_real_
  ),
  lincomb_est_ci(x_once_35 - x_once_34) %>%
    mutate(habitat = "once-logged", metric = "annual_slope"),
  lincomb_est_ci(x_rest_35 - x_rest_34) %>%
    mutate(habitat = "restored", metric = "annual_slope"),
  tibble(
    habitat = "primary",
    metric = "acd_at_35",
    estimate = prime_est,
    lwr95 = prime_est - 1.96 * prime_se,
    upr95 = prime_est + 1.96 * prime_se
  ),
  lincomb_est_ci(x_once_35) %>%
    mutate(habitat = "once-logged", metric = "acd_at_35"),
  lincomb_est_ci(x_rest_35) %>%
    mutate(habitat = "restored", metric = "acd_at_35")
) %>%
  mutate(source = "philipson_refit")

# LIDAR-derived sensitivity: implied annual delta ACD (Mg C ha^-1 yr^-1) by plot class.
# Plotted on the *same* panel as annual slopes (green), not a separate facet.
# Mapping to x-axis habitat:
#   liana_cutting -> restored
#   enrichment, sbe_once_logged_control, sbe_once_logged_control_plus_enrichment -> once-logged
#     (each with a short text label above the point)
#   primary -> primary; twice_logged -> twice-logged
habitat_levels <- c("primary", "once-logged", "restored", "twice-logged")

lidar_path <- file.path("Inputs", "lidar_SBE_class_summary.csv")
stopifnot(file.exists(lidar_path))
lidar_raw <- read.csv(lidar_path, stringsAsFactors = FALSE)

lidar_tbl <- lidar_raw %>%
  mutate(
    habitat = case_when(
      plot_class == "primary" ~ "primary",
      plot_class == "enrichment" ~ "once-logged",
      plot_class == "liana_cutting" ~ "restored",
      plot_class == "sbe_once_logged_control" ~ "once-logged",
      plot_class == "twice_logged" ~ "twice-logged",
      plot_class == "sbe_once_logged_control_plus_enrichment" ~ "once-logged",
      TRUE ~ NA_character_
    ),
    lidar_label = case_when(
      plot_class == "enrichment" ~ "enrichment",
      plot_class == "sbe_once_logged_control" ~ "SBE control",
      plot_class == "sbe_once_logged_control_plus_enrichment" ~ "SBE control + enrichment",
      TRUE ~ NA_character_
    ),
    metric = "annual_slope",
    source = "lidar",
    estimate = mean_implied_delta_carbon_MgC_ha_yr,
    lwr95 = lwr_95_MgC_ha_yr,
    upr95 = upr_95_MgC_ha_yr
  ) %>%
  filter(!is.na(habitat)) %>%
  mutate(
    habitat = factor(habitat, levels = habitat_levels),
    plot_class = factor(plot_class, levels = unique(plot_class)[order(plot_class)])
  ) %>%
  select(habitat, metric, source, estimate, lwr95, upr95, plot_class, lidar_label)

plot_tbl <- bind_rows(
  bind_rows(model_tbl, philipson_refit_tbl) %>%
    mutate(
      plot_class = NA_character_,
      lidar_label = NA_character_,
      habitat = factor(habitat, levels = habitat_levels),
      source = factor(source, levels = c("my_model", "philipson_refit"))
    ),
  lidar_tbl %>% mutate(plot_class = as.character(plot_class))
) %>%
  mutate(
    source = factor(source, levels = c("my_model", "philipson_refit", "lidar")),
    metric = factor(
      metric,
      levels = c("annual_slope", "acd_at_35"),
      labels = c(
        "Annual ACD increase (Mg C ha^-1 yr^-1)",
        "Predicted ACD at 35 years since logging (Mg C ha^-1)"
      )
    ),
    dodge_grp = ifelse(
      is.na(plot_class) | plot_class == "",
      paste0(as.character(source), "__", as.character(habitat)),
      paste0("lidar__", plot_class)
    )
  )

out_csv <- file.path(tab_dir, "carbon_recovery__quick_compare_model_vs_philipson.csv")
out_png <- file.path(fig_dir, "carbon_recovery__quick_compare_model_vs_philipson.png")
out_pdf <- file.path(fig_dir, "carbon_recovery__quick_compare_model_vs_philipson.pdf")

write.csv(as.data.frame(plot_tbl), out_csv, row.names = FALSE)

dodge_w <- 0.78
p <- ggplot(plot_tbl, aes(x = habitat, y = estimate, color = source, group = dodge_grp)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(
    aes(ymin = lwr95, ymax = upr95),
    width = 0.15,
    position = position_dodge(width = dodge_w),
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  geom_point(
    position = position_dodge(width = dodge_w),
    size = 2.5
  ) +
  geom_text(
    data = plot_tbl %>%
      dplyr::filter(!is.na(lidar_label)) %>%
      dplyr::mutate(label_y = dplyr::coalesce(upr95, estimate)),
    aes(
      x = habitat,
      y = label_y,
      group = dodge_grp,
      label = lidar_label
    ),
    inherit.aes = FALSE,
    position = position_dodge(width = dodge_w),
    vjust = -0.35,
    size = 2.4,
    color = "#238b45",
    lineheight = 0.95
  ) +
  scale_color_manual(
    values = c(
      "my_model" = "#1f78b4",
      "philipson_refit" = "#d95f02",
      "lidar" = "#2ca02c"
    ),
    labels = c(
      "my_model" = "My model",
      "philipson_refit" = "Philipson-equivalent refit",
      "lidar" = "LIDAR SBE (Inputs/lidar_SBE_class_summary)"
    ),
    drop = FALSE
  ) +
  facet_wrap(~ metric, scales = "free", ncol = 2) +
  labs(
    x = NULL,
    y = "Estimate",
    color = NULL,
    title = "My model vs Philipson-equivalent refit (LIDAR annual on left panel)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, size = 9)
  )

ggsave(out_png, p, width = 9.5, height = 5.2, dpi = 220)
ggsave(out_pdf, p, width = 9.5, height = 5.2)

message("Saved: ", out_csv)
message("Saved: ", out_png)
message("Saved: ", out_pdf)
