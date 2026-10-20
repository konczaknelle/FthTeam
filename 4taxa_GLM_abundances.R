# ============================================================
# MACROFAUNA ANALYSIS
# Taxa: Formicidae, Isoptera,
#       Lumbricina_oligochaetes, Coleoptera_sum
# ============================================================


# 1. DATA PREPARATION -----------------------------------------

library(ggplot2)
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)

data <- macrofaune_counts

data[, 1:6] <- lapply(data[, 1:6], as.factor)

data <- data %>%
  mutate(
    soil_type = factor(soil_type),
    layer = factor(layer, levels = c("L", "10", "20")),
    monolith_ID = factor(monolith_ID)
  )

str(data)


# 2. SELECT THE FOUR TAXA ------------------------------------

taxa_cols <- c(
  "Formicidae",
  "Isoptera",
  "Lumbricina_oligochaetes",
  "Coleoptera_sum"
)


# 3. COLOURBLIND-FRIENDLY PALETTE -----------------------------

cb_colors <- c(
  "#0072B2", # blue
  "#E69F00", # orange
  "#009E73", # green
  "#D55E00"  # vermillion
)

names(cb_colors) <- taxa_cols


# 4. MEAN DENSITY PER M2 -------------------------------------

data_long_mean <- data %>%
  pivot_longer(
    cols = all_of(taxa_cols),
    names_to = "taxon",
    values_to = "count"
  ) %>%
  group_by(soil_type, layer, taxon) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    density = mean_count * 16,
    taxon = factor(taxon, levels = taxa_cols)
  )

p_density <- ggplot(
  data_long_mean,
  aes(
    y = layer,
    x = density,
    fill = taxon
  )
) +
  geom_col(position = position_stack(reverse = TRUE)) +
  facet_wrap(~ soil_type) +
  scale_fill_manual(
    values = cb_colors,
    breaks = taxa_cols
  ) +
  scale_y_discrete(limits = c("20", "10", "L")) +
  labs(
    y = "Layer",
    x = expression("Mean density per m"^2),
    fill = "Taxon"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.title = element_text(face = "bold")
  )

print(p_density)


# 5. MEAN ln(COUNT + 1) --------------------------------------

data_long_log_mean <- data %>%
  pivot_longer(
    cols = all_of(taxa_cols),
    names_to = "taxon",
    values_to = "count"
  ) %>%
  mutate(
    ln_count = log(count + 1)
  ) %>%
  group_by(soil_type, layer, taxon) %>%
  summarise(
    mean_ln_count = mean(ln_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    taxon = factor(taxon, levels = taxa_cols)
  )

p_log <- ggplot(
  data_long_log_mean,
  aes(
    y = layer,
    x = mean_ln_count,
    fill = taxon
  )
) +
  geom_col(position = position_stack(reverse = TRUE)) +
  facet_wrap(~ soil_type) +
  scale_fill_manual(
    values = cb_colors,
    breaks = taxa_cols
  ) +
  scale_y_discrete(limits = c("20", "10", "L")) +
  labs(
    y = "Layer",
    x = "Mean ln(count + 1)",
    fill = "Taxon"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.title = element_text(face = "bold")
  )

print(p_log)


# 6. LONG-FORMAT DATA FOR STATISTICAL ANALYSIS ---------------

data_long <- data %>%
  pivot_longer(
    cols = all_of(taxa_cols),
    names_to = "taxon",
    values_to = "count"
  ) %>%
  mutate(
    taxon = factor(taxon, levels = taxa_cols),
    soil_type = factor(soil_type),
    layer = factor(layer, levels = c("L", "10", "20")),
    monolith_ID = factor(monolith_ID)
  )

str(data_long)


# 7. NEGATIVE-BINOMIAL MIXED MODELS --------------------------
#
# Separate model for each taxon:
# count ~ soil_type * layer + (1 | monolith_ID)
#
# The interaction tests whether layer effects differ
# between soil types.


models_nb <- list()
diagnostics_nb <- list()
anova_results <- list()
emm_layers <- list()
pairwise_layers <- list()


for (tx in taxa_cols) {
  
  cat("\n============================================\n")
  cat("TAXON:", tx, "\n")
  cat("============================================\n")
  
  dat_tx <- data_long %>%
    filter(
      taxon == tx,
      !is.na(count),
      !is.na(soil_type),
      !is.na(layer),
      !is.na(monolith_ID)
    )
  
  if (any(dat_tx$count < 0) ||
      any(dat_tx$count != floor(dat_tx$count))) {
    stop(
      "Invalid count values for taxon: ", tx
    )
  }
  
  # Fit negative-binomial GLMM
  m_nb <- glmmTMB(
    count ~ soil_type * layer + (1 | monolith_ID),
    family = nbinom2(link = "log"),
    data = dat_tx
  )
  
  models_nb[[tx]] <- m_nb
  
  # Model summary
  print(summary(m_nb))
  
  # DHARMa residual diagnostics
  res_nb <- simulateResiduals(
    fittedModel = m_nb,
    n = 1000
  )
  
  diagnostics_nb[[tx]] <- res_nb
  
  plot(
    res_nb,
    main = paste("DHARMa diagnostics:", tx)
  )
  
  print(testDispersion(res_nb))
  print(testZeroInflation(res_nb))
  
  # Type II tests of fixed effects
  anova_results[[tx]] <- car::Anova(
    m_nb,
    type = 2
  )
  
  print(anova_results[[tx]])
  
  # Estimated marginal means for layers,
  # averaged across soil types
  emm_layers[[tx]] <- emmeans(
    m_nb,
    ~ layer
  )
  
  print(emm_layers[[tx]])
  
  # Tukey-adjusted pairwise layer comparisons
  pairwise_layers[[tx]] <- as.data.frame(
    pairs(
      emm_layers[[tx]],
      adjust = "tukey"
    )
  )
  
  print(pairwise_layers[[tx]])
}


# 8. COMBINE PAIRWISE COMPARISONS -----------------------------
pairwise_layers_df <- bind_rows(
  lapply(names(pairwise_layers), function(tx) {
    pairwise_layers[[tx]] %>%
      mutate(taxon = tx, .before = 1)
  })
)

print(pairwise_layers_df)


# 9. COMBINE ANOVA RESULTS ------------------------------------

anova_results_df <- bind_rows(
  lapply(names(anova_results), function(tx) {
    as.data.frame(anova_results[[tx]]) %>%
      rownames_to_column("term") %>%
      mutate(taxon = tx, .before = 1)
  })
)

print(anova_results_df)


# ============================================================
# END OF ANALYSIS
# ============================================================


