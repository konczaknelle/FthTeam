data <- macrofaune_counts
data[, 1:6] <- lapply(data[, 1:6], as.factor)
str(data)



# Colorblind-friendly palette
# Extended colourblind-friendly palette (28 colours)
cb_colors <- c(
  "#0072B2", # blue
  "#E69F00", # orange
  "#009E73", # green
  "#D55E00", # vermillion
  "#CC79A7", # purple
  "#56B4E9", # sky blue
  "#F0E442", # yellow
  "#000000", # black
  "#117733", # dark green
  "#332288", # dark blue-purple
  "#88CCEE", # light blue
  "#DDCC77", # tan
  "#AA4499", # magenta-purple
  "#44AA99", # teal
  "#CC6677", # rose
  "#999933", # olive
  "#882255", # dark magenta
  "#661100", # dark brown
  "#6699CC", # muted blue
  "#117744", # forest green
  "#AA5500", # dark orange
  "#7755AA", # violet
  "#228833", # green
  "#AA3377", # pink-purple
  "#4477AA", # steel blue
  "#EE6677", # coral
  "#66CCEE", # cyan
  "#BBBB44"  # yellow-green
)
#caluclate density by m^2
library(ggplot2)
library(tidyverse)

taxa_cols <- names(data)[7:34]

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
    density = mean_count * 16
  )


# Alphabetical order for the legend
data_long_mean <- data_long_mean %>%
  mutate(
    taxon = factor(taxon, levels = sort(unique(taxon)))
  )

ggplot(
  data_long_mean,
  aes(
    y = factor(layer, levels = c("20", "10", "L")),
    x = density,
    fill = taxon
  )
) +
  geom_col(position = position_stack(reverse = TRUE)) +
  facet_wrap(~ soil_type) +
  scale_fill_manual(values = cb_colors) +
  labs(
    y = "Layer",
    x = "Mean density per m^2",
    fill = "Taxon"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.title = element_text(face = "bold")
  )

######################
# Long format, ln(x + 1) transformation,
# then mean by soil type, layer and taxon
data_long_mean <- data %>%
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
  )

# Alphabetical order for the legend
data_long_mean <- data_long_mean %>%
  mutate(
    taxon = factor(taxon, levels = sort(unique(taxon)))
  )

ggplot(
  data_long_mean,
  aes(
    y = factor(layer, levels = c("20", "10", "L")),
    x = mean_ln_count,
    fill = taxon
  )
) +
  geom_col(position = position_stack(reverse = TRUE)) +
  facet_wrap(~ soil_type) +
  scale_fill_manual(values = cb_colors) +
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

##########################


###GLM densities----
library(dplyr)
library(tidyr)
library(glmmTMB)
library(DHARMa)

taxa_cols <- names(data)[7:34]

data_long <- data %>%
  pivot_longer(
    cols = all_of(taxa_cols),
    names_to = "taxon",
    values_to = "count"
  )
library(dplyr)
library(tidyr)
library(glmmTMB)
library(DHARMa)

taxa_cols <- names(data)[7:34]

data_long <- data %>%
  pivot_longer(
    cols = all_of(taxa_cols),
    names_to = "taxon",
    values_to = "count"
  ) %>%
  mutate(
    taxon = factor(taxon)
  )
str(data_long)

data_total <- data %>%
  mutate(
    total_count = rowSums(
      across(all_of(taxa_cols)),
      na.rm = FALSE
    ),
    soil_type = factor(soil_type),
    layer = factor(layer, levels = c("L", "10", "20")),
    monolith_ID = factor(monolith_ID)
  ) %>%
  mutate(
    density = total_count * 16
  )


library(glmmTMB)

m_poisson <- glmmTMB(
  total_count ~ soil_type * layer +
    (1 | monolith_ID),
  family = poisson(link = "log"),
  data = data_total
)

summary(m_poisson)

res <- simulateResiduals(m_poisson)
plot(res)
testDispersion(res)




m_nb <- glmmTMB(
  total_count ~ soil_type + layer,
  family = poisson,
  data = data_total
)

res_nb <- simulateResiduals(m_nb, n = 1000)
plot(res_nb)

summary(m_nb)

testDispersion(res_nb)
testZeroInflation(res_nb)










m_nb_int <- glmmTMB(
  total_count ~ soil_type * layer + (1|monolith_ID),
  family = nbinom2(link = "log"),
  data = data_total
)

summary(m_nb_int)

res_nb_int <- simulateResiduals(m_nb_int, n = 1000)
plot(res_nb_int)

summary(m_nb_int)

testDispersion(res_nb_int)
testZeroInflation(res_nb_int)
library(car)
Anova(m_nb_int)





install.packages("emmeans")
library(emmeans)

# Estimated marginal means for layers
emm_layers <- emmeans(m_nb, ~ layer)

# Pairwise comparisons between layers
pairs(emm_layers, adjust = "tukey")
#To get the results as a data frame:


pairwise_layers <- as.data.frame(
  pairs(emm_layers, adjust = "tukey")
)

pairwise_layers
