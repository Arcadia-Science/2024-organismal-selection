setwd("~/Documents/Research/github/2024-organismal-selection/")
source("code/utils.R")

###################
##### Load data#####
###################
##### Phylogeny#####
# Load phylogeny
tree <- read.newick("data/congruified_spprax_species_tree.newick")

# Calculate cophenetic distance (will be used for downstream analyses)
cophen <- cophenetic.phylo(tree)

# Extract human row (used as a gauge of distance from humans per-species)
cophen <- cophen[grep(
  "Homo-sapiens",
  rownames(cophen)
), ]

# Create tree without human (will be used for downstream analyses)
tree_nh <- drop.tip(tree, "Homo-sapiens")

##### Metadata#####
# Load
metadata <- read.csv("data/metadata.csv")
# Match colors with taxa
taxa_colors <- all_colors[1:length(unique(metadata$taxogroup2_unieuk))]
names(taxa_colors) <- unique(metadata$taxogroup2_unieuk)

# Match with species
species_colors <- taxa_colors[match(
  metadata$taxogroup2_unieuk,
  names(taxa_colors)
)]
names(species_colors) <- metadata$genus.species

###### Conservation#####
# Load
conservation_table <-
  readRDS("data/conservation_table.RDS")
conservation_species <-
  readRDS("data/conservation_by_spp.RDS")
conservation_species_best <-
  readRDS("data/conservation_best_protein_per_spp.RDS")

# Split conservation_table by gene family
conservation <- split(
  conservation_table,
  conservation_table$gene_family
)

###############################################################################
##### Plot phylogeny ###########################################################
##### (Related to figure 2; organism icons added later via illustrator)#########
###############################################################################
# Create metadata for tree (colors with be based on 'taxogroup2...')
meta <- data.frame(
  name = tree$tip.label,
  taxa = metadata$taxogroup2_unieuk[
    match(
      gsub("_", "-", tree$tip.label),
      metadata$genus.species
    )
  ],
  row.names = tree$tip.label
)

# Combine with tree via 'gettreedata'
tree_annot <- gettreedata(
  tree,
  meta
)

# Create initial plot
p <- ggtree(tree_annot,
  layout = "fan",
  open.angle = 45
)

cols <- c("taxa")
col <- cols[1]

# Reorder metadata by taxa
df <- meta[tree$tip.label, ][col]

# Match colors to tree order
colors <- species_colors[match(
  tree$tip.label,
  names(species_colors)
)]
# Use taxa as color names
names(colors) <- df$taxa

# Plot all
p + new_scale_fill() +
  geom_tippoint(
    mapping = aes(fill = .data[[col]]),
    size = 4,
    shape = 21,
    stroke = 0
  ) +
  geom_tiplab(offset = 0.5) +
  scale_fill_manual(
    values = colors,
    na.value = "white"
  )

####################################################################
##### Plot phylogeny with n gene families per species################
##### (Related to figure 3)##########################################
####################################################################
##### Plot phylogeney with gene family n#####
# Calculate n ogs per species
n_ogs_per_species <- unlist(lapply(
  conservation_species,
  function(x) {
    length(unique((x$gene_family)))
  }
))

# Combine tree and gene family number via 'phylo4d'
p4d <- phylo4d(
  tree,
  unique_species
)

# Match taxa colors to tree
cols <- species_colors[match(
  tree$tip.label,
  names(species_colors)
)]

# Plot using 'barplot.phylo4d' from the 'phylosignal' package
barplot.phylo4d(p4d,
  center = FALSE,
  scale = FALSE,
  bar.col = cols,
  tree.ladderize = TRUE,
  trait.labels = "n ogs",
  use.edge.length = FALSE
)

##### Compare vertebrates to non-vertebrates#####
## Compare gene family n of vertebrates to all other species
# Vertebrate gene family n
verts <- n_ogs_per_species[names(n_ogs_per_species) %in%
  metadata$genus.species[
    metadata$taxogroup2_unieuk == "Vertebrata"
  ]]

# Non-vertebrate gene family n
non_verts <- n_ogs_per_species[!names(n_ogs_per_species) %in%
  metadata$genus.species[
    metadata$taxogroup2_unieuk == "Vertebrata"
  ]]

# Compare with Kruskal-Wallis test
verts_ks <- kruskal.test(list(
  verts,
  non_verts
))
verts_ks

#####################################################################
##### Compare gene family n with phylogenetic distance################
##### (Related to figure 4)###########################################
#####################################################################
###### Figure 4A: Plot cophenetic distance and gene family number#####
# Regress gene family number and phylogenetic distance
plot(cophen,
  n_ogs_per_species[match(
    names(cophen),
    names(n_ogs_per_species)
  )],
  ylab = "n gene families",
  xlab = "Distance from humans",
  cex.lab = 1.5,
  cex.axis = 1.5,
  pch = 20,
  cex = 2,
  col = species_colors[match(
    names(cophen),
    names(species_colors)
  )]
)

# Add species names
text(
  cophen,
  n_ogs_per_species[match(
    names(cophen),
    names(n_ogs_per_species)
  )],
  names(cophen)
)

# Add regression line
mod <- lm(n_ogs_per_species[match(
  names(cophen),
  names(n_ogs_per_species)
)] ~ cophen)

abline(mod,
  lty = "dashed"
)

##### Compare the observed gene family distribution#####
##### to expected using permutations (Figure 4B)########
# Calculate n proteins per gene family
n_proteins_per_og <- unlist(lapply(
  conservation,
  function(x) {
    nrow(x)
  }
))

# Calculate n species per OG
n_species_per_og <- unlist(lapply(
  conservation,
  function(x) {
    length(unique((x$nonref_species)))
  }
))

# Extract fitted values from regression comparing gene family n and distance
fitted <- mod$fitted.values

# Normalize by sum of all fitted values; converts to a proportion
fitted <- fitted / sum(fitted)

# Create n observations per species (i.e. proportion of 10000 observations
# determined by 'fitted'); each permutation will randomly sample from
# this distribution
n <- 100000 * fitted

pool <- c()
for (i in 1:length(n)) {
  pool <- c(pool, rep(names(n)[i], n[i]))
}

# Create empty list for saving permutation results
perms <- list()

# Counter
pb <- txtProgressBar(
  min = 1,
  max = length(seq(
    min(n_proteins_per_og),
    max(n_proteins_per_og), 5
  )),
  style = 3,
  width = 100,
  char = "."
)

# Set seed
set.seed(666)

# Set up n permutations
n_perms <- 100

# Loop over range encompassing the minimum and maximum number of proteins
# per gene family
for (i in seq(
  min(n_proteins_per_og),
  max(n_proteins_per_og), 5
)) {
  # Update counter
  setTxtProgressBar(pb, i)

  # Create empty vector to save individual permutations
  lengths <- c()

  # Permutate
  for (j in 1:n_perms) {
    # Calculate the number of species present per sample
    lengths <- c(
      lengths,
      length(table(pool[sample(1:length(pool), i)]))
    )
  }

  # Add to results list
  perms[[as.character(i)]] <- lengths
}

# Compare saturation point of permuted and real data
# (i.e. which gene family size contains proteins from all 63 species)
perm_saturation <- unlist(lapply(
  perms,
  function(x) length(grep(63, x))
))
obs_saturation <- split(
  n_species_per_og,
  n_proteins_per_og
)
obs_saturation <- unlist(lapply(
  obs_saturation,
  function(x) length(grep(63, x))
))

# Filter to gene families with at least permutation containing all species
perm_saturation <- perm_saturation[perm_saturation > 0]
obs_saturation <- obs_saturation[obs_saturation > 0]

# Find smallest gene family size that contains all species
min(as.numeric(names(perm_saturation)))
# 264
min(as.numeric(names(obs_saturation)))
# 70

##### Plot Figures 4B + 4C#####
# Compare n proteins per og and n species
smoothScatter(log(n_proteins_per_og),
  n_species_per_og,
  nrpoints = 0,
  xlab = "n proteins per gene family (log)",
  ylab = "n species per gene family",
  cex.lab = 1.5,
  cex.axis = 1.5,
  bty = "n",
  nbin = 100,
  colramp = colorRampPalette(c(
    "white",
    rev(
      unlist(
        arcadia_magma$color_dict
      )
    )
  ))
)
lines(log(as.numeric(names(perms))),
  unlist(lapply(perms, function(x) mean(x))),
  lwd = 1.5
)

## Figure 4C: compare the number of proteins within each gene family
## to evolutionary distance
# Calculate phylogenetic dispersion per OG
phylo_dist_per_og <- list()
for (i in 1:length(conservation)) {
  # Get species in OG
  species_og <- unique(conservation[[i]]$nonref_species)

  # Get MRCA
  mrca_og <- getMRCA(
    tree,
    species_og[species_og %in% tree$tip.label]
  )

  # Subset tree
  tree_og <- extract.clade(
    tree,
    mrca_og
  )

  # Calculate "age" (max branch length)
  phylo_dist_per_og[[names(conservation)[i]]] <- max(tree_og$edge.length)
}
phylo_dist_per_og <- unlist(phylo_dist_per_og)

# Plot
smoothScatter(n_species_per_og,
  phylo_dist_per_og,
  nrpoints = 0,
  xlab = "n species per OG",
  ylab = "Time span per OG (billion years)",
  cex.lab = 1.5,
  cex.axis = 1.5,
  bty = "n",
  nbin = 100,
  colramp = colorRampPalette(c(
    "white",
    rev(
      unlist(
        arcadia_magma$color_dict
      )
    )
  ))
)

#########################################################################
##### Plot molecular conservation distributions for all gene families#####
##### (Related to figure 5)###############################################
#########################################################################
# Calculate distributions of molecular conservation (via histograms)
# across all gene families
conservation_hists <- lapply(
  conservation_species_best,
  function(x) {
    hist(x$trait_dist,
      plot = FALSE,
      breaks = seq(0, 50, 1)
    )$density
  }
)

# Normalize histograms by maximum value (to make comparable across all families)
conservation_hists <- lapply(
  conservation_hists,
  function(x) x / max(x)
)

# Combine
conservation_hists <- do.call(
  rbind,
  conservation_hists
)

# Hierarchically cluster gene family distributions
hcl <- hclust(dist(conservation_hists))

# Plot
plot(hcl, labels = FALSE)

## Plot gene family conservation distributions as points
# Calculate median conservation (for initiating plot)
stats <- lapply(
  conservation_species_best,
  function(x) median(x$trait_dist)
)
stats <- stats[hcl$order]

# Order gene families based on hierarchical clustering
hcl_conservation <- conservation_species_best[hcl$order]

# Plot
par(bg = NA)
plot(unlist(lapply(stats, function(x) x$stats[3])),
  ylim = c(0, 50),
  pch = 20,
  ylab = "Distance from human",
  xlab = "",
  cex.axis = 1.5,
  cex.lab = 1.5
)
for (i in 1:length(stats)) {
  points(rep(i, length(x[[i]]$trait_dist)),
    hcl_conservation[[i]]$trait_dist,
    pch = 20,
    cex = 2,
    col = scales::alpha("black", 0.1)
  )
}

####################################
##### Plot example gene families#####
##### (Related to figure 6)##########
####################################
##### Calculate slopes for each gene family#####
# Calculate slopes (will be used to categorize gene families)
pb <- txtProgressBar(
  min = 1,
  max = length(conservation_species_best),
  style = 3,
  width = 100,
  char = "."
)

og_slopes <- list()
for (i in 1:length(conservation_species_best)) {
  # Update counter
  setTxtProgressBar(pb, i)

  # Split on species
  s <- split(
    conservation_species_best[[i]],
    conservation_species_best[[i]]$nonref_species
  )

  # Select most conserved protein per species
  s <- lapply(s, function(x) x[which.min(x$trait_dist), ])

  # Recombine
  s <- do.call(rbind, s)

  # Get species phylo distances
  dists <- cophen[match(
    s$nonref_species,
    names(cophen)
  )]

  # Regression
  mod <- lm(s$trait_dist ~ dists)
  slope <- coef(mod)[2]
  r2 <- summary(mod)[8]

  # Add to list
  og_slopes[[names(conservation_species_best)[i]]] <- list(
    data = s,
    r2 = r2,
    slope = slope
  )
}

# Extract model fit, slope, and conservation per gene family
r2s <- unlist(lapply(og_slopes, function(x) x$r2))[
  unlist(lapply(og_slopes, function(x) nrow(x$data)) > 40)
]

slopes <- unlist(lapply(og_slopes, function(x) x$slope))[
  unlist(lapply(og_slopes, function(x) nrow(x$data)) > 40)
]

cons <- unlist(lapply(og_slopes, function(x) mean(x$data$trait_dist)))[
  unlist(lapply(og_slopes, function(x) nrow(x$data)) > 40)
]

cons_var <- unlist(lapply(og_slopes, function(x) {
  sd(x$data$trait_dist) / mean(x$data$trait_dist)
}))[
  unlist(lapply(og_slopes, function(x) nrow(x$data)) > 40)
]

# Combine in a data frame for filtering
res <- data.frame(
  protein = names(og_slopes)
  [unlist(lapply(og_slopes, function(x) nrow(x$data)) > 40)],
  r2 = r2s,
  cons = cons,
  cons_var = cons_var,
  slope = slopes
)

##### Positive relationship#####
# Plot examples
par(mfrow = c(2, 2))
# Positive relationship
x <- res[res$slope > 0, ]
x <- x[x$cons < 10, ]
x <- x[which.max(x$r2), ]
n <- x$protein

dists <- cophen[
  match(
    og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
    names(cophen)
  )
]

plot(dists,
  og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist,
  xlab = "Cophenetic distance",
  ylab = "Conservation with humans",
  cex.axis = 1.5,
  cex.lab = 1.5,
  pch = 20,
  cex = 2,
  ylim = c(0, 3),
  col = species_colors[
    match(
      og_slopes[[grep(
        n,
        names(og_slopes)
      )]]$data$nonref_species,
      names(species_colors)
    )
  ]
)
abline(lm(og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist ~ dists),
  lty = "dashed",
  lwd = 1.5
)
text(
  3000,
  3,
  paste("r2 =", signif(x$r2, 2), sep = "")
)
title(
  main = paste(n, "(PTN4)"),
  font.main = 1,
  cex.main = 1.5
)

##### Positive relationship, human specific#####
x <- res[res$slope > 0, ]
x <- x[which.max(x$r2), ]
n <- x$protein

dists <- cophen[match(
  og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
  names(cophen)
)]

plot(dists,
  og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist,
  xlab = "Cophenetic distance",
  ylab = "Conservation with humans",
  ylim = c(0, 20),
  cex.axis = 1.5,
  cex.lab = 1.5,
  pch = 20,
  cex = 2,
  col = species_colors[
    match(
      og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
      names(species_colors)
    )
  ]
)

abline(lm(og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist ~ dists),
  lty = "dashed",
  lwd = 1.5
)
text(
  3000,
  20,
  paste("r2 =", signif(x$r2, 2), sep = "")
)
title(
  main = paste(n, "(FOXA1)"),
  font.main = 1,
  cex.main = 1.5
)

###### Deep conservation#####
x <- res[res$cons < 1, ]
x <- x[which.min(x$cons_var), ]
n <- x$protein

dists <- cophen[match(
  og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
  names(cophen)
)]

plot(dists,
  og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist,
  xlab = "Cophenetic distance",
  ylab = "Conservation with humans",
  cex.axis = 1.5,
  cex.lab = 1.5,
  pch = 20,
  ylim = c(0, 2),
  cex = 2,
  col = species_colors[
    match(
      og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
      names(species_colors)
    )
  ]
)
abline(lm(og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist ~ dists),
  lty = "dashed",
  lwd = 1.5
)
text(
  3000,
  2,
  paste("r2 =", signif(x$r2, 2), sep = "")
)
title(
  main = paste(n, "(ARF3)"),
  font.main = 1,
  cex.main = 1.5
)

##### Negative relationship#####
x <- res[res$slope < 0, ]
x <- x[which.max(x$r2), ]
n <- x$protein

dists <- cophen[match(
  og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
  names(cophen)
)]

plot(dists,
  og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist,
  xlab = "Cophenetic distance",
  ylab = "Conservation with humans",
  cex.axis = 1.5,
  cex.lab = 1.5,
  pch = 20,
  cex = 2,
  ylim = c(0, 5),
  col = species_colors[
    match(
      og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
      names(species_colors)
    )
  ]
)

abline(lm(og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist ~ dists),
  lty = "dashed",
  lwd = 1.5
)
text(
  3000,
  5,
  paste("r2 =", signif(x$r2, 2), sep = "")
)
title(
  main = paste(n, "(3HIDH)"),
  font.main = 1,
  cex.main = 1.5
)

############################################
##### Generate and plot phylomorphospace#####
##### (Related to Figure 7A)#################
############################################
##### Create and plot phylomorphospace#####
# Create matrix of conservation values for each gene family
conservation_matrix <- as.data.frame(matrix(
  ncol =
    length(conservation_species_best),
  nrow =
    63
))

# Column and rownames
colnames(conservation_matrix) <- names(conservation_species_best)
rownames(conservation_matrix) <- tree_nh$tip.label

# Add in conservation values
for (i in 1:length(conservation_species_best)) {
  x <- conservation_species_best[[i]]$trait_dist[
    match(
      rownames(conservation_matrix),
      conservation_species_best[[i]]$nonref_species
    )
  ]
  conservation_matrix[, i] <- x
}

# Replace NAs with maximum conservation distance
conservation_matrix[is.na(conservation_matrix)] <-
  max(conservation_table$trait_dist)

# PCA of conservation matrix
pca <- prcomp(conservation_matrix)

# Match PC order to tree tip labels
x <- pca$x[match(
  tree_nh$tip.label,
  rownames(pca$x)
), 1:2]

# Get species colors
cols <- species_colors[match(
  rownames(x),
  names(species_colors)
)]
cols <- c(cols, rep("grey90", tree_nh$Nnode))
names(cols) <- 1:(length(tree_nh$tip) + tree_nh$Nnode)

# Plot phylomorphospace
phylomorphospace(tree_nh,
  x,
  ftype = "off",
  node.by.map = TRUE,
  node.size = 1.5,
  control = list(col.node = cols),
  xlab = paste("PC1", "(45.89%)"),
  ylab = paste("PC2", ("7.71%")),
  cex.axis = 1.5,
  cex.lab = 1.5
)

##### Analyze PCs#####
# Calculate correlation between gene family number and PC1
res <- apply(
  pca$x,
  2,
  function(y) {
    cor.test(
      n_ogs_per_species,
      y[match(
        names(n_ogs_per_species),
        rownames(pca$x)
      )]
    )
  }
)

# ID which PCs are significantly correlated
# Adjust expected p-value for number of tests (63)
pnorm <- 0.05 / 63
which(unlist(lapply(
  res,
  function(x) x$p.value
)) < pnorm)
# PC1
# 1

# Correlation between phylogenetic distance and PCs
# Match PCs to phylogenetic distance vector
x <- cophen[match(
  rownames(pca$x),
  names(cophen)
)]

# Correlate
res <- apply(
  pca$x,
  2,
  function(y) {
    cor.test(x, y)
  }
)

# ID which PCs are significantly correlated
# Adjust expected p-value for number of tests (63)
pnorm <- 0.05 / 63
which(unlist(lapply(
  res,
  function(x) x$p.value
)) < pnorm)
# PC1
# 1

###############################
##### Compare Elo ratings#######
##### (Related to figure 7B-C)##
###############################
##### Calculate per-species Elo ratings#####
# Create a matrix of all possible species matchups
matchups <- expand.grid(
  unique(conservation_table$nonref_species),
  unique(conservation_table$nonref_species)
)

# Loop over gene families and collect conservation for species pairs
# Create empty list to save matchup outcomes into
matchup_outcomes <- list()
pb <- txtProgressBar(
  min = 1,
  max = length(conservation_species_best),
  style = 3,
  width = 100,
  char = "."
)
for (i in 1:length(conservation_species_best)) {
  # Update counter
  setTxtProgressBar(pb, i)

  # Get all combinations
  x <- expand.grid(
    unique(conservation_species_best[[i]]$nonref_species),
    unique(conservation_species_best[[i]]$nonref_species)
  )

  # Convert to data frame
  x <- data.frame(
    spp1 = as.character(x$Var1),
    spp2 = as.character(x$Var2),
    spp1_cons = rep(NA, nrow(x)),
    spp2_cons = rep(NA, nrow(x))
  )

  # Remove self matches
  x <- x[!x[, 1] == x[, 2], ]

  # Add outcomes
  for (j in 1:nrow(x)) {
    x$spp1_cons[j] <- min(conservation_species_best[[i]]$rank_trait_dist_norm[
      grep(
        x$spp1[j],
        conservation_species_best[[i]]$nonref_species
      )
    ])

    x$spp2_cons[j] <- min(conservation_species_best[[i]]$rank_trait_dist_norm[
      grep(
        x$spp2[j],
        conservation_species_best[[i]]$nonref_species
      )
    ])
  }

  # Add to list
  matchup_outcomes[[names(conservation_species_best)[i]]] <- x
}

# Get protein pairs for all gene families with at least 10 matchups
protein_sample <- do.call(rbind, lapply(
  matchup_outcomes,
  function(x) {
    if (nrow(x) >= 10) {
      x[sample(1:nrow(x), 10), ]
    }
  }
))

## Calculate Elo ratings over n permutations, using a random match order
## each time
# Set number of permutations
n_perms <- 50

# Create empty lists for saving results
all_elo_iterations <- list()
all_species_elo_distributions <- list()

# Set set seed
set.seed(1234)

# Loop over n permutations
for (g in 1:n_perms) {
  # Update counter
  print(paste(g, "out of", n_perms))

  # Calculate elo scores n times
  mean_elo_scores <- list()

  # Subsample outcomes
  test <- protein_sample[sample(1:nrow(protein_sample), 10000), ]

  for (h in 1:n_perms) {
    # Reorder
    test <- test[sample(1:nrow(test)), ]

    # Remove self
    test <- test[!test[, 1] == test[, 2], ]

    # Create matrix containing elo for each species
    # (to be updated w/ each match)
    unique_species <- unique(test$spp1)
    elo_scores <- as.data.frame(
      rep(
        1500,
        length(unique_species)
      ),
      row.names = unique_species
    )

    # Create species elo lists to keep score records
    species_elo <- split(
      rep(1500, length(unique_species)),
      unique_species
    )

    # Loop over and simulate outcomes for each match
    pb <- txtProgressBar(
      min = 1,
      max = nrow(test),
      style = 3,
      width = 100,
      char = "."
    )

    for (i in 1:nrow(test)) {
      # Update counter
      setTxtProgressBar(pb, i)

      # Get outcome
      outcome <- unlist(test[i, 3:4])
      names(outcome) <- test[i, 1:2]

      # Convert to wins
      if (outcome[1] == outcome[2]) {
        outcome <- c(0.5, 0.5)
      } else {
        x <- which.min(outcome)
        y <- which.max(outcome)
        outcome[x] <- 1
        outcome[y] <- 0
      }

      # Get probabilities
      probs <- elo_scores[match(
        c(test[i, 1:2]),
        rownames(elo_scores)
      ), 1]

      # Calculate elo
      elo_update <- elo.calc(outcome,
        rep(probs[1], 2),
        rep(probs[2], 2),
        k = 4
      )[1, ]
      names(elo_update) <- names(outcome)

      # Update elo matrix
      elo_scores[grep(
        names(elo_update)[1],
        rownames(elo_scores)
      ), 1] <- elo_update[1]

      elo_scores[grep(
        names(elo_update)[2],
        rownames(elo_scores)
      ), 1] <- elo_update[2]

      # Update species scores
      x <- species_elo[[grep(
        names(elo_update)[1],
        names(species_elo)
      )]]
      x <- unlist(c(x, elo_update[1]))
      species_elo[[grep(
        names(elo_update)[1],
        names(species_elo)
      )]] <- x

      x <- species_elo[[grep(
        names(elo_update)[2],
        names(species_elo)
      )]]
      x <- unlist(c(x, elo_update[2]))
      species_elo[[grep(
        names(elo_update)[2],
        names(species_elo)
      )]] <- x

      # Add to list
      mean_elo_scores[[as.character(h)]] <- as.data.frame(elo_scores)
      all_species_elo_distributions[[as.character(h)]] <- species_elo
    }
  }
  # Calculate mean elo score per species
  n <- as.character(unique(matchups$Var1))

  mean_elo_scores <- do.call(
    cbind,
    lapply(
      mean_elo_scores,
      function(x) {
        x[match(
          n,
          rownames(x)
        ), 1]
      }
    )
  )
  rownames(mean_elo_scores) <- n

  # Plot distribution of matches
  plot(
    approx(1:length(species_elo[[1]]),
      species_elo[[1]],
      n = 100
    )$y,
    type = "l",
    ylim = c(1400, 1650),
    xlim = c(0, 100),
    ylab = "Elo score",
    xlab = "Match number",
    cex.axis = 1.5,
    cex.lab = 1.5,
    col = species_colors[match(
      names(species_elo)[1],
      names(species_colors)
    )]
  )
  abline(h = 1500, lty = "dashed", lwd = 1.5)
  for (i in 1:length(species_elo)) {
    lines(
      approx(1:length(species_elo[[i]]),
        species_elo[[i]],
        n = 100
      )$y,
      col = species_colors[match(
        names(species_elo)[i],
        names(species_colors)
      )]
    )
  }
  # Add to list
  all_elo_iterations[[as.character(g)]] <- mean_elo_scores
}

# Combine into matrix
all_elo_iterations <- do.call(
  cbind,
  all_elo_iterations
)

##### Plot mean Elo ratings#####
# Calculate mean elo per species
elo_mean <- rowMeans(all_elo_iterations)

# Sort by value
elo_mean <- sort(elo_mean)

# Compare elo ratings to phylogenetic distance from human
plot(cophen[match(names(elo_mean), names(cophen))],
  elo_mean,
  type = "n",
  xlab = "Cophenetic distance",
  ylab = "Elo score",
  cex.axis = 1.5,
  cex.lab = 1.5
)
abline(
  h = 1500,
  lty = "dashed",
  lwd = 1.5
)
points(cophen[match(names(elo_mean), names(cophen))],
  elo_mean,
  pch = 21,
  bg = species_colors[match(
    names(elo_mean),
    names(species_colors)
  )],
  col = darken_color(species_colors[match(
    names(elo_mean),
    names(species_colors)
  )]),
  cex = 1.5
)
text(
  cophen[match(names(elo_mean), names(cophen))],
  elo_mean,
  names(elo_mean)
)

##### Compare species' relative probabilities#####
# Compare elo ratings of vertebrates to rest
elo_taxa <- metadata$taxogroup2_unieuk[match(
  names(elo_mean),
  metadata$genus.species
)]
median(elo_mean[elo_taxa == "Vertebrata"])
# 1571.394
median(elo_mean[!elo_taxa == "Vertebrata"])
# 1478.744

# Relative probability between Abeoforma and Chimpanzee
elo::elo.prob(
  elo_mean[grep("Pan", names(elo_mean))],
  elo_mean[grep("Abeoforma", names(elo_mean))]
)
# [1] 0.7639637

# Relative probability between Chlorella and Chlamydomonas
elo::elo.prob(
  elo_mean[grep("Chlorella", names(elo_mean))],
  elo_mean[grep("Chlamydomonas", names(elo_mean))]
)
# [1] 0.6745585

##### Compare Elo ratings by phylogenetic age#####
##' young' (cophenetic distance < 1130), 'middle' (>1130, <2000),
## and 'old' (>2000)
# Generate empty vector
elo_age <- rep(
  NA,
  length(elo_mean)
)

# Order phylogenetic distance to match Elo ratings
elo_cophen <- cophen[match(
  names(elo_mean),
  names(cophen)
)]

# Add in age classifications to vector
elo_age[elo_cophen < 1130] <- "young"
elo_age[elo_cophen > 1130 & x < 2000] <- "mid"
elo_age[elo_cophen > 2000] <- "old"

# Split Elo ratings by age classification
elo_age <- split(
  elo_mean,
  elo_age
)

# Compare 'young' species to others
kruskal.test(list(
  elo_age$young,
  c(
    elo_age$mid,
    elo_age$old
  )
))

# Compare 'mid' and 'old'
kruskal.test(list(
  elo_age$mid,
  elo_age$old
))

# Invertebrates compared to 'old' (cophenetic distance >2000)
elo_inv <- rep(
  NA,
  length(elo_mean)
)

# Add labels to vector
elo_inv[elo_cophen < 3000] <- "old"
elo_inv[elo_cophen < 1700 & elo_cophen > 1200] <- "invert"

# Split Elo ratings by classification
elo_inv <- split(
  elo_mean,
  elo_inv
)

# Compare
kruskal.test(list(
  elo_inv$invert,
  elo_inv$old
))

##### Predict Elo ratings from gene family number#####
# Regression predicting Elo ratings using gene family number
mod <- lm(elo_mean[match(
  names(n_ogs_per_species),
  names(elo_mean)
)] ~
  n_ogs_per_species)

# Extract studentized residuals
res <- sort(MASS::studres(mod))

# Extract residuals by taxa
res_taxa <- metadata$taxogroup1_unieuk[match(
  names(res),
  metadata$genus.species
)]

# Split residuals into taxa
res_taxa <- split(
  res,
  res_taxa
)

# Calculate median residual per taxa
sort(unlist(lapply(res_taxa, function(x) median(x))))

############################################
##### Analyze Reactome pathway variation#####
##### (Related to Figure 8)##################
############################################
##### Load reactome#####
# Load reactome pathways
reactome <- read.delim("data/UniProt2Reactome.txt",
  header = FALSE
)

# Extract human pathways
reactome <- reactome[grep(
  "Homo sapiens",
  reactome[, 6]
), ]

# Add column names
colnames(reactome) <- c(
  "uniprot",
  "id",
  "reactome_id",
  "pathway",
  "iea",
  "species"
)

# Split on pathway
reactome <- split(
  reactome,
  reactome$pathway
)

# Create human protein list for each species
conservation_human_proteins <- lapply(
  conservation_species,
  function(x) {
    unique(x$ref_protein)
  }
)
conservation_human_proteins[["Homo-sapiens"]] <-
  unique(unlist(conservation_human_proteins))

##### Calculate %conservation in pathways for#####
##### each species using permutation tests########
# Create empty list to save results into
pathway_conservation <- list()

# Set n permutations
n_perms <- 100

# Set seed
set.seed(666)

# Loop over species and perform permutationt test
for (i in 1:length(conservation_human_proteins)) {
  # Print species name
  print(names(conservation_human_proteins)[i])

  # Extract species proteins
  species <- conservation_human_proteins[[i]]

  # Extract human proteins
  human <- conservation_human_proteins$`Homo-sapiens`

  # Generate empty vectors for results
  p_vals <- c()
  eff_size <- c()
  n_proteins <- c()

  # Set counter
  pb <- txtProgressBar(
    min = 1,
    max = length(reactome),
    style = 3,
    width = 100,
    char = "."
  )

  # Loop over all reactome pathways and compute conservation
  for (j in 1:length(reactome)) {
    # Update counter
    setTxtProgressBar(pb, j)

    # Extract reactome pathway
    pathway <- reactome[[j]]

    # Identify proteins in pathways
    in_pathway <- species[
      species %in% pathway$uniprot
    ]

    # Extract conservation scores
    hits <- conservation_species[[i]][
      conservation_species[[i]]$ref_protein %in% in_pathway,
    ]

    # Take most conserved protein when multiple hits
    hits <- unlist(lapply(
      split(
        hits$trait_dist,
        hits$ref_protein
      ),
      function(x) x[which.min(x)]
    ))

    # Empty vector for permutation results
    perms <- c()

    # Run permutations
    for (k in 1:n_perms) {
      perms <- c(perms, mean(conservation_table$trait_dist
      [sample(
          1:nrow(conservation_table),
          length(hits)
        )]))
    }

    # Identify outliers with cohen's d*
    if (length(hits) > 3) {
      eff_size <- c(eff_size, effsize::cohen.d(hits, perms)$estimate)
    } else {
      eff_size <- c(eff_size, NA)
    }

    # Calculate p-value
    if (is.null(hits)) {
      p_vals <- c(
        p_vals,
        NA
      )
    } else {
      p_vals <- c(
        p_vals,
        sum(mean(hits, na.rm = TRUE) < perms) / n_perms
      )
    }

    # n_proteins
    n_proteins <- c(
      n_proteins,
      length(hits)
    )
  }

  # Calculate percentage conservation for each pathway
  percent_conserved <- unlist(lapply(
    reactome,
    function(x) {
      sum(species %in% x$uniprot) / sum(human %in% x$uniprot)
    }
  ))

  print(mean(eff_size, na.rm = TRUE))

  # Add to results
  pathway_conservation[[names(conservation_human_proteins)[i]]] <- list(
    p_value = p_vals,
    eff_size = eff_size,
    n_proteins = n_proteins,
    percent_conserved = percent_conserved
  )
}

# Save
saveRDS(
  pathway_conservation,
  "out/reactome_pathway_conservation.RDS"
)

# Plot relationship between percent conserved and effect size
plot(
  unlist(lapply(pathway_conservation, function(x) {
    mean(x$percent_conserved, na.rm = TRUE)
  })),
  unlist(lapply(pathway_conservation, function(x) {
    mean(x$eff_size, na.rm = TRUE)
  })),
  ylab = "Mean conservation",
  xlab = "Mean %conserved",
  cex.axis = 1.5,
  cex.lab = 1.5,
  pch = 20,
  col = species_colors[match(
    names(pathway_conservation),
    names(species_colors)
  )]
)
text(
  unlist(lapply(pathway_conservation, function(x) {
    mean(x$percent_conserved, na.rm = TRUE)
  })),
  unlist(lapply(pathway_conservation, function(x) {
    mean(x$eff_size, na.rm = TRUE)
  })),
  names(pathway_conservation)
)
