# ----------------------------------------------------------------------------
# Unified script: allometry comparison + volume-based validation
# ----------------------------------------------------------------------------

library(tidyverse)
library(stringr)

set.seed(123)

# ============================================================
# USER INPUT
# ============================================================

albizia_models_to_use <- c("pow")   # "pow", "ln", or c("pow","ln")
albizia_model_choice  <- "pow"      # used for validation

# BEF: stand inventory volume × wood density is stem (commercial) wood mass; values
# >1 expand to whole-tree *aboveground* biomass to account for branches/bark/ foliage
# not in the volume table. We also compare BEF=1 (no expansion) as a sensitivity plot.
BEF_values <- c(low = 1.1, mid = 1.25, high = 1.4)

# ============================================================
# 1) LOAD + PREP DATA
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
    YEAR = month / 12,
    species_group = case_when(
      Species == "AF" ~ "Albizia",
      Species == "EP" ~ "Eucalyptus",
      TRUE ~ NA_character_
    )
  ) %>%
  
  # 🔑 Correct conditional filtering (once only)
  filter(
    (species_group == "Eucalyptus" & YEAR < 7) |
      (species_group == "Albizia"   & YEAR < 12)
  ) %>%  
  rename(volume_m3_ha  ="Volume..m3.ha." )

# ============================================================
# 2) PARAMETERS
# ============================================================

wd_albizia    <- 0.406
wd_eucalyptus <- 0.629

# ============================================================
# 3) BUILD ALL MODELS (ONCE)
# ============================================================

SSB <- SSB %>%
  mutate(
    wood_density = case_when(
      species_group == "Albizia" ~ wd_albizia,
      species_group == "Eucalyptus" ~ wd_eucalyptus
    ),
    
    # --------------------------------------------------------
    # CHAVE (pantropical)
    # --------------------------------------------------------
    tree_biomass_chave =
      0.0673 * (wood_density * dbh^2 * height)^0.976,
    
    AGB_chave =
      (tree_biomass_chave * live_trees_ha) / 1000,
    
    # --------------------------------------------------------
    # EUCALYPTUS MODEL
    # https://onlinelibrary.wiley.com/doi/epdf/10.1111/gcb.13201
    # --------------------------------------------------------
    tree_biomass_euc =
      if_else(species_group == "Eucalyptus",
              exp(-2.016 + 2.375 * log(dbh)) * 1.067,
              NA_real_),
    
    AGB_euc =
      (tree_biomass_euc * live_trees_ha) / 1000
  )

# ---- Albizia models (conditional) ----
if ("pow" %in% albizia_models_to_use) {
  SSB <- SSB %>%
    mutate(
      # https://www.tandfonline.com/doi/full/10.1080/21580103.2023.2256355
      tree_biomass_albizia_pow =
        if_else(species_group == "Albizia",
                0.08062 * dbh^2.36816,
                NA_real_),
      
      AGB_albizia_pow =
        (tree_biomass_albizia_pow * live_trees_ha) / 1000
    )
}

if ("ln" %in% albizia_models_to_use) {
  SSB <- SSB %>%
    mutate(
      # https://doi.org/10.1088/1755-1315/1506/1/012018
      tree_biomass_albizia_ln =
        if_else(species_group == "Albizia",
                exp(1.4926 + 0.068 * dbh + 0.0456 * height),
                NA_real_),
      
      AGB_albizia_ln =
        (tree_biomass_albizia_ln * live_trees_ha) / 1000
    )
}

# dplyr::case_when() evaluates all RHS expressions; keep missing columns so
# validation does not reference AGB_albizia_ln when only "pow" was fitted.
if (!"ln" %in% albizia_models_to_use) {
  SSB <- SSB %>% mutate(
    tree_biomass_albizia_ln = NA_real_,
    AGB_albizia_ln = NA_real_
  )
}
if (!"pow" %in% albizia_models_to_use) {
  SSB <- SSB %>% mutate(
    tree_biomass_albizia_pow = NA_real_,
    AGB_albizia_pow = NA_real_
  )
}

# ============================================================
# 4) ALLOMETRY COMPARISON (MODEL vs MODEL)
# ============================================================

model_cols <- c(
  if ("pow" %in% albizia_models_to_use) "AGB_albizia_pow",
  if ("ln"  %in% albizia_models_to_use) "AGB_albizia_ln",
  "AGB_euc"
) 

comparison_table <- SSB %>%
  select(species_group, AGB_chave, all_of(model_cols)) %>%
  pivot_longer(
    cols = all_of(model_cols),
    names_to = "model",
    values_to = "AGB"
  ) %>%
  filter(!is.na(AGB)) %>%
  filter(
    (species_group == "Albizia"   & grepl("albizia", model)) |
      (species_group == "Eucalyptus" & model == "AGB_euc")
  ) %>%
  group_by(species_group, model) %>%
  summarise(
    mean_ratio = mean(AGB / AGB_chave),
    correlation = cor(AGB, AGB_chave),
    n = n(),
    .groups = "drop"
  )

print(comparison_table)

# ---- Plot ----
fig_allometry <- SSB %>%
  select(species_group, AGB_chave, all_of(model_cols)) %>%
  pivot_longer(cols = all_of(model_cols), names_to = "model", values_to = "AGB") %>%
  filter(!is.na(AGB)) %>%
  mutate(
    model = recode(model,
                   AGB_albizia_pow = "Albizia power",
                   AGB_albizia_ln  = "Albizia log-linear",
                   AGB_euc         = "Eucalyptus-specific")
  ) %>%
  ggplot(aes(AGB_chave, AGB)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ paste(species_group, model, sep = " — "), scales = "free") +
  theme_bw()

print(fig_allometry)

# ============================================================
# 5) VOLUME-BASED VALIDATION
# ============================================================

SSB_val <- SSB %>%
  mutate(
    AGB_volume_noBEF = volume_m3_ha * wood_density,
    AGB_volume_low  = volume_m3_ha * wood_density * BEF_values["low"],
    AGB_volume_mid  = volume_m3_ha * wood_density * BEF_values["mid"],
    AGB_volume_high = volume_m3_ha * wood_density * BEF_values["high"],
    
    AGB_species = case_when(
      species_group == "Albizia"   & albizia_model_choice == "pow" ~ AGB_albizia_pow,
      species_group == "Albizia"   & albizia_model_choice == "ln"  ~ AGB_albizia_ln,
      species_group == "Eucalyptus" ~ AGB_euc
    )
  )

comparison_df <- SSB_val %>%
  select(species_group, AGB_chave, AGB_species,
         AGB_volume_low, AGB_volume_mid, AGB_volume_high) %>%
  pivot_longer(
    cols = starts_with("AGB_volume"),
    names_to = "scenario",
    values_to = "AGB_volume"
  ) %>%
  mutate(
    scenario = recode(
      scenario,
      AGB_volume_low = paste0("BEF_", unname(BEF_values["low"])),
      AGB_volume_mid = paste0("BEF_", unname(BEF_values["mid"])),
      AGB_volume_high = paste0("BEF_", unname(BEF_values["high"]))
    )
  ) %>%
  filter(!is.na(AGB_species), !is.na(AGB_volume))

# ---- Metrics ----
metrics <- comparison_df %>%
  group_by(species_group, scenario) %>%
  summarise(
    bias_chave   = mean(AGB_chave - AGB_volume),
    bias_species = mean(AGB_species - AGB_volume),
    rmse_chave   = sqrt(mean((AGB_chave - AGB_volume)^2)),
    rmse_species = sqrt(mean((AGB_species - AGB_volume)^2)),
    .groups = "drop"
  )

print(metrics)

# ---- Validation plot (BEF = 1; volume × WD only, no biomass expansion) ----
comparison_no_bef <- SSB_val %>%
  select(species_group, AGB_chave, AGB_species, AGB_volume_noBEF) %>%
  filter(!is.na(AGB_species), !is.na(AGB_volume_noBEF))

metrics_no_bef <- comparison_no_bef %>%
  group_by(species_group) %>%
  summarise(
    bias_chave = mean(AGB_chave - AGB_volume_noBEF, na.rm = TRUE),
    bias_species = mean(AGB_species - AGB_volume_noBEF, na.rm = TRUE),
    rmse_chave = sqrt(mean((AGB_chave - AGB_volume_noBEF) ^ 2, na.rm = TRUE)),
    rmse_species = sqrt(mean((AGB_species - AGB_volume_noBEF) ^ 2, na.rm = TRUE)),
    .groups = "drop"
  )

message("Validation vs volume (BEF = 1: no expansion, stem-wood from V×WD only):")
print(metrics_no_bef)

fig_validation_no_bef <- comparison_no_bef %>%
  pivot_longer(
    cols = c(AGB_chave, AGB_species),
    names_to = "model",
    values_to = "AGB_model"
  ) %>%
  mutate(
    model = recode(
      model,
      AGB_chave = "Chave (pantropical)",
      AGB_species = "Species allometry"
    )
  ) %>%
  ggplot(aes(AGB_volume_noBEF, AGB_model)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_grid(species_group ~ model, scales = "free") +
  labs(
    x = "Reference from inventory: volume x wood density, no BEF (Mg dry mass ha-1)",
    y = "Allometric AGB (Mg ha-1)",
    title = "Validation: allometry vs V x wood density (BEF = 1, no expansion)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10))

print(fig_validation_no_bef)

# ---- Validation plot ----
fig_validation <- comparison_df %>%
  pivot_longer(
    cols = c(AGB_chave, AGB_species),
    names_to = "model",
    values_to = "AGB_model"
  ) %>%
  ggplot(aes(AGB_volume, AGB_model)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 200)) +
  facet_grid(species_group ~ model + scenario, scales = "free") +
  theme_bw()

print(fig_validation)

# ============================================================
# SAVE
# ============================================================

out_dir <- file.path("Outputs", "plantation_allometrics_compare")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(out_dir, "fig_allometry_1to1.png"),
  fig_allometry,
  width = 8,
  height = 4,
  dpi = 200
)
ggsave(
  file.path(out_dir, "fig_validation.png"),
  fig_validation,
  width = 10,
  height = 6,
  dpi = 200
)
ggsave(
  file.path(out_dir, "fig_validation_noBEF.png"),
  fig_validation_no_bef,
  width = 8,
  height = 5,
  dpi = 200
)
message("Saved figures (PNG) to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))