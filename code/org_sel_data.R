setwd("~/Documents/Research/github/2024-organismal-selection/")
source("code/utils.R")

##### Phylogeny#####
# Load
tree <- read.newick("data/congruified_spprax_species_tree.newick")

# Calculate cophenetic distance (will be used for downstream analyses)
cophen <- cophenetic.phylo(tree)

# Create tree without human (will be used for downstream analyses)
tree_nh <- drop.tip(tree, "Homo-sapiens")

##### Protein conservation#####
# Set working directory (TO DO: CHANGE ONCE DEPOSITED IN ZENODO)
setwd("~/Documents/Research/arcadia-projects/raas-organism-prioritization/conservation_score_v2_06012024/gf-aa-multivar-distances/final_protein_pair_summary_tables/")

# List files (each corresponds to a gene family ('OG...'))
files <- list.files()

# Load all gene families and add to 'conservation'
conservation <- list()
for (i in 1:length(files)) {
  # Load file and add to list
  conservation[[gsub(
    "_final_summary_table.tsv",
    "",
    files[i]
  )]] <- read.delim(files[i])

  # Normalize trait distance rank by gene family size (useful for filtering)
  conservation[[i]]$rank_trait_dist_norm <-
    conservation[[i]]$rank_trait_dist / max(conservation[[i]]$rank_trait_dist)
}

# Combine gene families into a table
conservation_table <- do.call(
  rbind,
  conservation
)

# Split table by species
conservation_species <- split(
  conservation_table,
  conservation_table$nonref_species
)

# Filter keep just the best ortholog per species for each human protein
# Split by human protein
conservation_species_best <- split(
  conservation_table,
  conservation_table$ref_protein
)

# Loop over and select most conserved protein for each species
for (i in 1:length(conservation_species_best)) {
  # Split on species ID
  species_best <- split(
    conservation_species_best[[i]],
    conservation_species_best[[i]]$nonref_species
  )

  # Select most conserved protein per species
  species_best <- do.call(
    rbind,
    lapply(
      species_best,
      function(y) y[which.min(y$trait_dist), ]
    )
  )

  # Add to list
  conservation_species_best[[names(conservation)[i]]] <- species_best
}

##### Save cleaned data#####
saveRDS(
  conservation_table,
  "data/conservation_table.RDS"
)
saveRDS(
  conservation_species,
  "data/conservation_by_spp.RDS"
)
saveRDS(
  conservation_species_best,
  "data/conservation_best_protein_per_spp.RDS"
)
