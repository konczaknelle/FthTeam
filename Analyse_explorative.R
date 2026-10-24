# =============================================================================
# MACROFAUNE DU SOL - PARACOU - Alt (drainage vertical) vs DLS (drainage lateral)
# 01 - EXPLORATION : tableaux recapitulatifs, densites, boxplots, ACP
#
# Dispositif : 2 sites Alt + 2 sites DLS
#              x 5 monolithes TSBF par site
#              x 3 couches (litiere, 0-10 cm, 10-20 cm)
#              = 4 x 5 x 3 = 60 echantillons
#
# Ce script ne fait AUCUN test statistique : il sert a regarder les donnees
# avant de modeliser (c'est l'etape "analyse exploratoire").
# =============================================================================


# -----------------------------------------------------------------------------
# 0. PACKAGES ET PARAMETRES
# -----------------------------------------------------------------------------
# Si un package manque : install.packages("nom_du_package")
library(dplyr)     # manipuler les tableaux (group_by, summarise, mutate...)
library(tidyr)     # changer la forme des tableaux (pivot_longer, pivot_wider)
library(ggplot2)   # faire les figures

# Nom du fichier (il doit etre dans le "working directory" : getwd() pour le voir)
fichier <- "data_brut.csv"

# Taille d'un monolithe TSBF : 25 cm x 25 cm = 0.0625 m2
# (a verifier : c'est la valeur du protocole TSBF standard)
aire_monolithe <- 0.25 * 0.25

# Taxons qu'on ne veut PAS analyser.
# "indet" = individus non identifies. On peut ajouter ici d'autres colonnes,
# par exemple c("indet", "Acarina", "Collembola") si decide que la
# meso-faune (trop petite pour un tri a la main) ne fait pas partie de l'etude.
taxons_a_exclure <- c("indet")

# Les 3 "ingenieurs de l'ecosysteme" (pour la partie ingenieurs).
# Pour ajouter les coleopteres plus tard : c(..., "Coleoptera_sum")
taxons_ingenieurs <- c("Formicidae", "Isoptera", "Lumbricina_oligochaetes")

# Couleurs des deux types de sol
couleurs_type <- c("Alt" = "#2a9d8f", "DLS" = "#e76f51")

# Dossier de sortie pour les tableaux et figures
dir.create("output2", showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 1. IMPORT ET NETTOYAGE
# -----------------------------------------------------------------------------
# Le CSV utilise le point-virgule comme separateur -> read.csv2()
# "UTF-8-BOM" evite que le premier nom de colonne soit mal lu
donnees <- read.csv2(fichier, fileEncoding = "UTF-8-BOM")

# Les colonnes "descriptives" (qui decrivent l'echantillon)
colonnes_info <- c("monolith_ID", "sample_ID", "soil_type", "layer",
                   "date", "identification")

# Toutes les autres colonnes sont des taxons (comptages d'individus)
taxons <- setdiff(names(donnees), colonnes_info)

# Une case vide = taxon absent = 0
donnees[taxons][is.na(donnees[taxons])] <- 0

# On cree des variables "propres" avec des facteurs ordonnes :
#   - type    : Alt ou DLS
#   - couche  : litiere / 0-10 / 10-20 (dans cet ordre, de haut en bas)
#   - site    : Alt_1, Alt_2, DLS_1, DLS_2 (= monolith_ID sans le dernier numero)
donnees <- donnees %>%
  mutate(type   = factor(soil_type, levels = c("alt", "dls"),
                         labels = c("Alt", "DLS")),
         couche = factor(layer, levels = c("L", "10", "20"),
                         labels = c("Litiere", "0-10 cm", "10-20 cm")),
         site   = sub("_[0-9]+$", "", monolith_ID))

# --- Verifications : le dispositif doit etre exactement celui prevu ---
cat("\n--- Verification du dispositif ---\n")
cat("Nombre d'echantillons (attendu 60) :", nrow(donnees), "\n")
cat("Nombre de monolithes  (attendu 20) :", n_distinct(donnees$monolith_ID), "\n")
cat("Echantillons par site x couche (attendu 5 partout) :\n")
print(table(donnees$site, donnees$couche))

# Individus non identifies / exclus : combien on en perd ?
cat("\nIndividus exclus de l'analyse :\n")
print(colSums(donnees[taxons_a_exclure]))


# -----------------------------------------------------------------------------
# 2. LE TABLEAU "LONG" : UNE LIGNE = UN TAXON DANS UN ECHANTILLON
# -----------------------------------------------------------------------------
# Le fichier d'origine est "large" (une colonne par taxon, 60 lignes).
# Pour les figures et les modeles, on a besoin du format "long" :
# une ligne par combinaison echantillon x taxon (60 x ~28 taxons).
# C'est dans ce tableau que "taxon" devient une variable a part entiere.

long <- donnees %>%
  pivot_longer(cols = all_of(taxons),
               names_to = "taxon",
               values_to = "abondance") %>%          # abondance = nb d'individus
  filter(!taxon %in% taxons_a_exclure) %>%           # on retire les taxons exclus
  mutate(densite = abondance / aire_monolithe)       # ind / m2

# Taxons conserves (sans les exclus)
taxons_gardes <- setdiff(taxons, taxons_a_exclure)


# -----------------------------------------------------------------------------
# 3. TABLEAUX RECAPITULATIFS
# -----------------------------------------------------------------------------

# 3a. Un tableau par taxon : combien d'individus, quel pourcentage, et dans
#     combien d'echantillons il est present (= prevalence).
tab_taxons <- long %>%
  group_by(taxon) %>%
  summarise(total_individus = sum(abondance),
            nb_echantillons_presents = sum(abondance > 0),   # sur 60
            densite_moyenne_m2 = mean(densite)) %>%          # moyenne des 60 echantillons
  mutate(pourcent_abondance = 100 * total_individus / sum(total_individus),
         pourcent_presence  = 100 * nb_echantillons_presents / nrow(donnees)) %>%
  arrange(desc(total_individus))

print(as.data.frame(tab_taxons), digits = 3)
write.csv(tab_taxons, "output2/tab_taxons.csv", row.names = FALSE)

# 3b. Pourcentage d'abondance de chaque taxon DANS chaque type de sol
tab_type <- long %>%
  group_by(type, taxon) %>%
  summarise(total = sum(abondance), .groups = "drop") %>%
  group_by(type) %>%
  mutate(pourcent = 100 * total / sum(total)) %>%    # % au sein de Alt, puis au sein de DLS
  ungroup() %>%
  pivot_wider(names_from = type, values_from = c(total, pourcent)) %>%
  arrange(desc(total_Alt + total_DLS))

print(as.data.frame(tab_type), digits = 3)
write.csv(tab_type, "output2/tab_pourcent_par_type.csv", row.names = FALSE)

# 3c. Densite moyenne (ind/m2) par type x couche x taxon
tab_type_couche <- long %>%
  group_by(type, couche, taxon) %>%
  summarise(densite_moyenne_m2 = mean(densite),
            ecart_type = sd(densite), .groups = "drop")
print(as.data.frame(tab_type_couche), digits = 3)
write.csv(tab_type_couche, "output2/tab_densite_type_couche.csv", row.names = FALSE)

# 3d. Un tableau "par echantillon" avec les variables de synthese :
#     abondance totale, densite totale, richesse (= nb de taxons differents)
par_echantillon <- long %>%
  group_by(sample_ID, monolith_ID, site, type, couche) %>%
  summarise(abondance_totale = sum(abondance),
            richesse = sum(abondance > 0), .groups = "drop") %>%
  mutate(densite_totale = abondance_totale / aire_monolithe)

# 3e. Densite d'un monolithe ENTIER (litiere + 0-10 + 10-20 additionnees),
#     c'est-a-dire le nombre d'individus sous 1 m2 de sol (profil 0-20 cm).
par_monolithe <- long %>%
  group_by(monolith_ID, site, type) %>%
  summarise(abondance_totale = sum(abondance), .groups = "drop") %>%
  mutate(densite_profil_m2 = abondance_totale / aire_monolithe)

cat("\nDensite totale du profil (litiere + 0-20 cm), ind/m2, par type de sol :\n")
par_monolithe %>%
  group_by(type) %>%
  summarise(moyenne = mean(densite_profil_m2), mediane = median(densite_profil_m2),
            min = min(densite_profil_m2), max = max(densite_profil_m2)) %>%
  print()


# -----------------------------------------------------------------------------
# 4. DIAGNOSTIC DES DONNEES (pour justifier le choix de la loi des modeles)
# -----------------------------------------------------------------------------
# Une loi de Poisson suppose variance = moyenne. Si la variance est tres
# superieure a la moyenne ("surdispersion"), on prend une binomiale negative.
# Ici on calcule le rapport variance / moyenne (Poisson => proche de 1)
# et le % de zeros pour chaque taxon ingenieur, par type x couche.
diagnostic <- long %>%
  filter(taxon %in% taxons_ingenieurs) %>%
  group_by(taxon, type, couche) %>%
  summarise(moyenne = mean(abondance),
            variance = var(abondance),
            variance_sur_moyenne = variance / moyenne,
            pourcent_zeros = 100 * mean(abondance == 0), .groups = "drop")

print(as.data.frame(diagnostic), digits = 3)


# -----------------------------------------------------------------------------
# 5. BOXPLOTS
# -----------------------------------------------------------------------------
# Un boxplot resume la distribution : trait epais = mediane, boite = 50 % des
# valeurs centrales (quartiles), moustaches = etendue "normale" des donnees.
# On ajoute les points reels (formes differentes selon le site) pour voir
# TOUTES les donnees, pas seulement le resume.
#
# Echelle "pseudo-log" : les abondances vont de 0 a plusieurs milliers
# d'ind/m2. Sur une echelle normale, tout serait ecrase pres de zero.
# La pseudo-log garde les zeros (contrairement au log) et etire le bas.

graduations <- c(0, 10, 100, 1000, 5000)    # valeurs affichees sur l'axe y
formes_sites <- c(16, 17, 15, 18)           # rond, triangle, carre, losange

# Un petit "modele de figure" qu'on reutilise : on lui donne un tableau,
# le nom de la colonne a representer, et le titre de l'axe y.
faire_boxplot <- function(tableau, y, titre_y) {
  ggplot(tableau, aes(x = couche, y = .data[[y]], fill = type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6,
                 position = position_dodge(0.8)) +
    geom_point(aes(shape = site),
               position = position_jitterdodge(jitter.width = 0.15,
                                               dodge.width = 0.8),
               size = 1.8, alpha = 0.8) +
    scale_fill_manual(values = couleurs_type, name = "Type de sol") +
    scale_shape_manual(values = formes_sites, name = "Site") +
    labs(x = NULL, y = titre_y) +
    theme_bw()
}

# 5a. Densite totale de macrofaune (tous taxons confondus)
fig_total <- faire_boxplot(par_echantillon, "densite_totale",
                           expression("Densite totale (ind."~m^-2~")")) +
  scale_y_continuous(trans = "pseudo_log", breaks = graduations)
print(fig_total)
ggsave("output2/fig1_densite_totale.png", fig_total, width = 7, height = 5, dpi = 300)

# 5b. Richesse taxonomique (nombre de taxons differents par echantillon)
fig_richesse <- faire_boxplot(par_echantillon, "richesse",
                              "Richesse (nb de taxons par echantillon)")
print(fig_richesse)
ggsave("output2/fig2_richesse.png", fig_richesse, width = 7, height = 5, dpi = 300)

# 5c. Les ingenieurs : un graphique par taxon (facet_wrap)
ingenieurs <- long %>%
  filter(taxon %in% taxons_ingenieurs) %>%
  mutate(taxon = recode(taxon, Formicidae = "Fourmis", Isoptera = "Termites",
                        Lumbricina_oligochaetes = "Vers de terre"))

fig_ingenieurs <- faire_boxplot(ingenieurs, "densite",
                                expression("Densite (ind."~m^-2~")")) +
  facet_wrap(~ taxon, scales = "free_y") +    # un panneau par taxon, axe y propre a chacun
  scale_y_continuous(trans = "pseudo_log", breaks = graduations)
print(fig_ingenieurs)
ggsave("output2/fig3_ingenieurs.png", fig_ingenieurs, width = 11, height = 5, dpi = 300)

# 5d. Composition : part de chaque taxon dans l'abondance, par type x couche.
#     Les taxons qui font moins de 2 % du total sont regroupes en "Autres".
taxons_majeurs <- tab_taxons$taxon[tab_taxons$pourcent_abondance >= 2]

composition <- long %>%
  mutate(taxon_groupe = ifelse(taxon %in% taxons_majeurs, taxon, "Autres")) %>%
  group_by(type, couche, taxon_groupe) %>%
  summarise(total = sum(abondance), .groups = "drop") %>%
  group_by(type, couche) %>%
  mutate(proportion = total / sum(total))

fig_compo <- ggplot(composition, aes(x = type, y = proportion, fill = taxon_groupe)) +
  geom_col(colour = "white", width = 0.7) +
  facet_wrap(~ couche) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Type de sol", y = "Part des individus", fill = "Taxon") +
  theme_bw()
print(fig_compo)
ggsave("output2/fig4_composition.png", fig_compo, width = 9, height = 5, dpi = 300)


# -----------------------------------------------------------------------------
# 6. ACP (Analyse en Composantes Principales)
# -----------------------------------------------------------------------------
# But : resumer 25+ taxons en 2 axes pour voir si les echantillons Alt et DLS
# se separent. Deux echantillons proches sur le graphique ont une faune
# de composition semblable.
#
# Precaution 1 - transformation de Hellinger :
#   On remplace chaque comptage par la racine carree de sa PROPORTION dans
#   l'echantillon. Sans cela, les fourmis (>300 ind. dans un echantillon)
#   ecraseraient tous les autres taxons dans l'ACP. Hellinger est la
#   transformation standard pour des donnees d'abondance avec beaucoup de zeros.
# Precaution 2 - on retire les echantillons vides (impossible de calculer
#   une proportion : 0/0) et les taxons presents dans moins de 3 echantillons
#   (trop rares pour structurer quoi que ce soit).

faire_acp <- function(matrice, info, titre) {
  # matrice : une ligne par echantillon (ou monolithe), une colonne par taxon
  # info    : tableau avec type, site (et couche) dans le MEME ordre de lignes
  
  vide <- rowSums(matrice) == 0
  matrice <- matrice[!vide, ]
  info <- info[!vide, ]
  matrice <- matrice[, colSums(matrice > 0) >= 3]
  
  hellinger <- sqrt(matrice / rowSums(matrice))        # transformation de Hellinger
  acp <- prcomp(hellinger, center = TRUE, scale. = FALSE)
  
  # % de variance expliquee par chaque axe
  variance <- 100 * acp$sdev^2 / sum(acp$sdev^2)
  
  # Coordonnees des echantillons (scores) et des taxons (fleches)
  scores <- cbind(info, as.data.frame(acp$x[, 1:2]))
  fleches <- as.data.frame(acp$rotation[, 1:2])
  fleches$taxon <- rownames(fleches)
  # On ne garde que les 8 taxons qui contribuent le plus aux 2 axes
  fleches <- fleches %>%
    mutate(longueur = sqrt(PC1^2 + PC2^2)) %>%
    slice_max(longueur, n = 8)
  # Mise a l'echelle des fleches pour qu'elles tiennent dans le graphique
  echelle <- 0.8 * max(abs(scores[, c("PC1", "PC2")])) / max(abs(fleches[, c("PC1", "PC2")]))
  
  figure <- ggplot(scores, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = 0, colour = "grey70") +
    stat_ellipse(aes(colour = type), level = 0.8, linetype = 2) +
    geom_point(aes(colour = type, shape = couche), size = 3) +
    geom_segment(data = fleches,
                 aes(x = 0, y = 0, xend = PC1 * echelle, yend = PC2 * echelle),
                 arrow = arrow(length = unit(0.2, "cm")), colour = "grey30") +
    geom_text(data = fleches,
              aes(x = PC1 * echelle * 1.1, y = PC2 * echelle * 1.1, label = taxon),
              size = 3) +
    scale_colour_manual(values = couleurs_type, name = "Type de sol") +
    labs(x = sprintf("Axe 1 (%.1f %%)", variance[1]),
         y = sprintf("Axe 2 (%.1f %%)", variance[2]),
         shape = "Couche", title = titre) +
    theme_bw()
  
  cat("\n", titre, "- variance expliquee par les 2 premiers axes :",
      round(variance[1] + variance[2], 1), "%\n")
  return(figure)
}

# 6a. ACP sur les 60 echantillons
#     Attendu : l'axe 1 separe surtout les COUCHES (litiere vs sol), pas les types.
fig_acp1 <- faire_acp(matrice = donnees[taxons_gardes],
                      info = donnees[c("type", "site", "couche")],
                      titre = "ACP - 60 echantillons")
print(fig_acp1)
ggsave("output2/fig5_acp_echantillons.png", fig_acp1, width = 8, height = 6.5, dpi = 300)

# 6b. ACP sur les 20 monolithes (les 3 couches additionnees).
#     On supprime l'effet "couche" pour comparer directement Alt et DLS.
monolithes <- donnees %>%
  group_by(monolith_ID, type, site) %>%
  summarise(across(all_of(taxons_gardes), sum), .groups = "drop") %>%
  mutate(couche = "toutes")          # colonne factice pour que la fonction marche

fig_acp2 <- faire_acp(matrice = monolithes[taxons_gardes],
                      info = monolithes[c("type", "site", "couche")],
                      titre = "ACP - 20 monolithes (couches cumulees)")
print(fig_acp2)
ggsave("output2/fig6_acp_monolithes.png", fig_acp2, width = 8, height = 6.5, dpi = 300)

cat("\nTermine. Tableaux et figures dans le dossier 'output2'.\n")

