# ----------------------------------------------------------------------------
# Unified script: allometry comparison + volume-based validation
# ----------------------------------------------------------------------------

library(tidyverse)
library(stringr)
library(cowplot)

set.seed(123)

# ============================================================
#  INPUTS
# ============================================================

albizia_models_to_use <- c("pow")   # "pow", "ln", or c("pow","ln")
albizia_model_choice  <- "pow"      # used for validation

# BEF: stand inventory volume × wood density is stem (commercial) wood mass; a
# BEF > 1 expands this to whole-tree aboveground biomass so it is comparable to
# Chave-style AGB estimates. We use one supplement-ready reference value only.
BEF_value <- 1.2

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
#from Jarapudin et al 2020
#https://research.usc.edu.au/esploro/outputs/journalArticle/Growth-performance-of-selected-taxa-as/99450828302621
wd_albizia    <- 0.406
wd_eucalyptus <- 0.629
#wd_eucalyptus <- 0.500

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
# Comparison of Chave vs specific-specific allometries

# ============================================================
# 5) VOLUME-BASED VALIDATION
# ============================================================
#compare allometry performance vs measured volume of merchantable bole 
SSB_val <- SSB %>%
  mutate(
    AGB_volume_noBEF = volume_m3_ha * wood_density,
    
    AGB_species = case_when(
      species_group == "Albizia"   & albizia_model_choice == "pow" ~ AGB_albizia_pow,
      species_group == "Albizia"   & albizia_model_choice == "ln"  ~ AGB_albizia_ln,
      species_group == "Eucalyptus" ~ AGB_euc
    )
  )

comparison_df <- SSB_val %>%
  mutate(
    AGB_volume = volume_m3_ha * wood_density * BEF_value,
    scenario = paste0("BEF_", BEF_value)
  ) %>%
  select(species_group, AGB_chave, AGB_species, AGB_volume, scenario) %>%
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

# ---- Validation plot ----
fig_validation <- comparison_df %>%
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
  ggplot(aes(AGB_volume, AGB_model)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 200)) +
  facet_grid(species_group ~ model, scales = "free") +
  labs(
    x = paste0("Inventory-derived AGB using merchantable volume, wood density, and BEF = ", BEF_value, " (Mg ha-1)"),
    y = "Allometric AGB (Mg ha-1)"
  ) +
  theme_bw()

print(fig_validation)

# ============================================================
# 6) SINGLE SUPPLEMENT FIGURE
# ============================================================

panel_a <- fig_allometry +
  labs(
    x = "Chave AGB (Mg ha-1)",
    y = "Species-specific AGB (Mg ha-1)"
  ) +
  theme(
    plot.title = element_blank(),
    strip.text = element_text(size = 9)
  )

panel_b <- fig_validation +
  theme(
    plot.title = element_blank(),
    strip.text = element_text(size = 9)
  )

supp_figure <- plot_grid(
  panel_a,
  panel_b,
  labels = c("A", "B"),
  ncol = 1,
  align = "v",
  rel_heights = c(0.9, 1.1)
)

# ============================================================
# SAVE
# ============================================================

out_dir <- file.path("Outputs", "plantation_allometrics_compare")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(out_dir, "fig_nature_supp__plantation_allometrics.png"),
  supp_figure,
  width = 10,
  height = 10,
  dpi = 200
)
message("Saved figures (PNG) to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
