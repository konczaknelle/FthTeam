################libraries----

library(vegan)
library(permute)
library(ggplot2)
library(ggrepel)
library(grid)
library(dplyr)

layer_colors <- c(
  "L"  = "#009E73",  
  "10" = "#E69F00",  
  "20" = "#5678A8"   
) # set at start of global script for all our plots

############
# import data once globally for complete script
###########
df <- macrofaune_counts

# defining taxonomic abundance columns
community_cols <- names(df)[7:ncol(df)]

# remove zero-abundance layers (decision following vegan handbook, since zero-abundance layer does not convey community composition information)
keep <- rowSums(df[, community_cols], na.rm = TRUE) > 0
df <- df[keep, , drop = FALSE]

# create the community matrix
comm <- as.matrix(df[, community_cols])
storage.mode(comm) <- "numeric" # making sure its num

# check dimensions (should be less for comm, should be 59 rows (60 - one zero-obs))
dim(df)
dim(comm)

# make sure qualitative variables are factors
df$monolith_ID <- factor(df$monolith_ID)
df$soil_type <- factor(df$soil_type, levels = c("alt", "dls"))
df$layer <- factor(df$layer, levels = c("L", "10", "20"))


######## analysis 1: PERMANOVA layer effects----

meta_layer <- df[, c("monolith_ID", "soil_type", "layer")] 

# Bray-Curtis PERMANOVA
# set seed for reproducible permutation chain
set.seed(123)

# define permutation scheme, limit WITHIN monoliths (bc spatial dependence)
perm_layer <- how(
  nperm = 9999,
  blocks = meta_layer$monolith_ID
) # chose large number of perms for stabilizing p value, increase if borderline

# run based on defined scheme
permanova_layer <- adonis2(
  comm ~ layer,
  data = meta_layer,
  method = "bray",
  permutations = perm_layer
)

print(permanova_layer)
# => significant, R2 16% 

# check dispersion usinf betadisper (if significant, permanova result might only be to unequal dispersion; a bit like homogeneity of variance in anova)
dist_layer <- vegdist(comm, method = "bray")

bd_layer <- betadisper(
  dist_layer,
  group = meta_layer$layer,
  type = "median",
  bias.adjust = TRUE
)

# set seed again for reproducibility, same permutation scheme
set.seed(123)

dispersion_layer <- permutest(
  bd_layer,
  permutations = perm_layer
)

print(dispersion_layer)
# not significant, permanova significance due to actual differences 

# quick visualisation (PCoA)
plot(bd_layer)

# pairwise PERMANOVA testing for layers
layer_pairs <- list(
  c("L", "10"),
  c("L", "20"),
  c("10", "20")
)

# function based to apply to all pairs
pairwise_results <- lapply(layer_pairs, function(pair) {
  
  # Keep only the two layers being compared
  keep_pair <- meta_layer$layer %in% pair
  
  comm_pair <- comm[keep_pair, , drop = FALSE]
  meta_pair <- droplevels(meta_layer[keep_pair, , drop = FALSE])
  
  # restrict perms to monoliths, like in global PERMANOVA
  perm_pair <- how(
    nperm = 9999,
    blocks = meta_pair$monolith_ID
  )
  
  # perform PERMANOVA
  set.seed(123)
  res <- adonis2(
    comm_pair ~ layer,
    data = meta_pair,
    method = "bray",
    permutations = perm_pair
  )
  #for nice table notation
  data.frame(
    comparison = paste(pair, collapse = " vs "),
    F = res$F[1],
    R2 = res$R2[1],
    p_value = res$`Pr(>F)`[1]
  )
})

pairwise_results <- do.call(rbind, pairwise_results)

# Correct for the three pairwise tests
pairwise_results$p_adjusted <- p.adjust(
  pairwise_results$p_value,
  method = "bonferroni" ) # BH correction also tested, in both corrections all pairs significant, more conservative bonferroni kept


print(pairwise_results)



######## analysis 2: PERMANOVA soil effects ----
# sum across all layers within each monolith => to test independent of layer effects, since dependent of monolith
comm_mono <- rowsum(
  comm,
  group = df$monolith_ID,
  reorder = FALSE
)

# one metadata row per monolith
meta_mono <- df[
  !duplicated(df$monolith_ID),
  c("monolith_ID", "soil_type")
]

# aligning matrix with metadata
comm_mono <- comm_mono[
  as.character(meta_mono$monolith_ID),
  ,
  drop = FALSE
]

# Bray-Curtis PERMANOVA
#set seed again
set.seed(123)

permanova_soil <- adonis2(
  comm_mono ~ soil_type,
  data = meta_mono,
  method = "bray",
  permutations = 9999
) 

print(permanova_soil)
# not significant, no differences by soil type

# betadisper
dist_soil <- vegdist(comm_mono, method = "bray")

bd_soil <- betadisper(
  dist_soil,
  group = meta_mono$soil_type,
  type = "median",
  bias.adjust = TRUE
)

set.seed(123)

dispersion_soil <- permutest(
  bd_soil,
  permutations = 9999
)

print(dispersion_soil)

# simple plot
plot(bd_soil)


########## PCoA visualization for layers ----
# calculate, using bray like in PERMANOVA
dist_layer <- vegdist(comm, method = "bray")
pcoa_layer <- cmdscale(
  dist_layer, k = 2, eig = TRUE, add = TRUE
) # set k for number of dimensions +> k =2 two dimensional representation

# sample scores
pcoa_df <- data.frame(
  PCoA1 = pcoa_layer$points[, 1],
  PCoA2 = pcoa_layer$points[, 2],
  layer = factor(df$layer, levels = c("L", "10", "20"))
)

# calculate hulls around samples
hull_df <- pcoa_df %>%
  group_by(layer) %>%
  slice(chull(PCoA1, PCoA2))

# fit taxon abundances to the PCoA (using same permutation scheme as in PERMANOVA)
fit <- envfit(
  pcoa_layer$points,
  comm,
  permutations = perm_layer
)

# extract taxon results and adjust p-values (since we do multiple testing)
taxa_results <- data.frame(
  taxon = rownames(fit$vectors$arrows),
  R2 = fit$vectors$r,
  p_value = fit$vectors$pvals
)

taxa_results$p_adjusted <- p.adjust(
  taxa_results$p_value,
  method = "bonferroni"
) # BH correction also tested, since no different result in taxa with p < 0.05, more conservative bonferroni kept

# extract taxon vectors 
vec <- as.data.frame(scores(fit, display = "vectors"))
vec$taxon <- rownames(vec)
names(vec)[1:2] <- c("PCoA1", "PCoA2")

# match adjusted p-values by taxon
vec$p_adjusted <- taxa_results$p_adjusted[
  match(vec$taxon, taxa_results$taxon)
]

# keep those with adjusted p < 0.05
vec <- vec[
  !is.na(vec$p_adjusted) &
    vec$p_adjusted < 0.05,
]

# calculate dim1 and dim2 percentages
eig <- pcoa_layer$eig
positive_eig <- eig[eig > 0]

pc1_pct <- 100 * eig[1] / sum(positive_eig)
pc2_pct <- 100 * eig[2] / sum(positive_eig)


# plot PCoA
p <- ggplot(pcoa_df, aes(PCoA1, PCoA2, color = layer)) +
  geom_polygon(
    data = hull_df,
    aes(x = PCoA1, y = PCoA2, fill = layer, group = layer),
    alpha = 0.12,
    color = NA,
    inherit.aes = FALSE
  ) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = layer_colors) +
  scale_fill_manual(values = layer_colors) +
  theme_classic(base_size = 12) +
  labs(
    title = "PCoA of macrofauna communities by layer",
    color = "Layer",
    fill = "Layer",
    x = sprintf("PCoA1 (%.1f%%)", pc1_pct),
    y = sprintf("PCoA2 (%.1f%%)", pc2_pct)
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )


# taxon arrows and labels
  p <- p +
    geom_segment(
      data = vec,
      aes(
        x = 0, y = 0,
        xend = PCoA1, yend = PCoA2
      ),
      inherit.aes = FALSE,
      color = "grey25",
      linewidth = 0.9,
      arrow = arrow(
        length = unit(0.2, "cm"),
        type = "closed"
      )
    ) +
    geom_label_repel(
      data = vec,
      aes(
        x = PCoA1, y = PCoA2,
        label = taxon
      ),
      inherit.aes = FALSE,
      color = "grey15",
      fill = alpha("white", 0.6),
      label.size = NA,
      size = 3.5,
      fontface = "italic",
      box.padding = 0.6,
      point.padding = 0.5,
      force = 3,
      max.overlaps = Inf,
      min.segment.length = Inf,
      segment.color = NA,
      seed = 123
    )

print(p)

##################### FIN ----