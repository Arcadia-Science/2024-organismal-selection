message("Loading the protein mv distance calculation source functions...")
# Load some remaining packages
suppressMessages(library(parallel))
suppressMessages(library(RcppParallel))
suppressMessages(library(tidyverse))
# Source the custom functions we use herein
setwd("./src/")
message(" 1. phylo_multivariate_distance_functions.R")
suppressMessages(source("./phylo_multivariate_distance_functions.R"))
message(" 2. protein_dist_permutation_tests.R")
suppressMessages(source("./protein_dist_permutation_tests.R"))
message(" 3. protein_distance_calculation_functions.R")
suppressMessages(source("./protein_distance_calculation_functions.R"))
setwd("../")

# Function to read in and get per-species gene copy number
# for each orthogroup
get_per_spp_og_counts <-
  function(results_dir = NULL, out_dir = NULL) {
    orthogroup_dir <-
      list.files(path = paste0(results_dir, "orthofinder/complete_dataset/"),
                 pattern = "Results_Inflation", full.names = TRUE)

    # read in per-species gene copy number per gene family
    og_counts <-
      read.table(paste0(orthogroup_dir,
                        "/Orthogroups/Orthogroups.GeneCount.tsv"),
                 header = TRUE, check.names = FALSE)
    colnames(og_counts) <- gsub("\\..*", "", colnames(og_counts))

    # And calculate the number of species in each gene family
    og_counts$NumSpecies <-
      rowSums(og_counts[, -c(1, ncol(og_counts))] > 0)

    # Write out to file if an output directory is provided
    if (!is.null(out_dir)) {
      # Create the directory if it doesn't yet exist
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      # And write out to file
      write.table(og_counts, file = paste0(out_dir, "og_copy_num_per_spp.tsv"),
                  sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = TRUE)
    }
    return(og_counts)
  }

pathd8_path <-
  paste0(Sys.getenv("CONDA_PREFIX"), "/bin:", Sys.getenv("CONDA_PREFIX"),
         "/lib64:", Sys.getenv("LD_LIBRARY_PATH"))
Sys.setenv(LD_LIBRARY_PATH = pathd8_path)

message("Done. Time-calibrating our species tree...")
# Set the noveltree results directory to path
noveltree_out_dir <-
  "2024-organismal-selection-zenodo/results-noveltree-model-euks/"

# Read in the asteroid species tree
tree <-
  ladderize(read.tree(paste0(noveltree_out_dir,
                             "speciesrax/species_trees/",
                             "inferred_species_tree.newick")))
tree$tip.label <- gsub("-", "_", tree$tip.label)

# And the timetree.org tree
tt <- read.tree("data/final-species-shortlist-timetree.nwk")

# First make sure the timetree is ultrametric, fully bifurcating,
# with no 0-length branches
tt$edge.length <- tt$edge.length + 0.00001
tt <- multi2di(tt)
log <- capture.output({ # Prevent the function from printing to screen
  tt <- force.ultrametric(tt, method = "extend")
})

# Now, see what congruifying the timetree tree looks like:
cong_tree <-
  geiger::congruify.phylo(reference = tt, target = tree,
                          scale = "treePL")

# Great. Use this one throughout.
spp_tree <- ladderize(cong_tree$phy)
spp_tree$tip.label <- gsub("_", "-", spp_tree$tip.label)

# Save out to file:
write.tree(spp_tree, "data/congruified_spprax_species_tree.newick")

# Read in metadat:
metadat <-
  read.table(paste0(noveltree_out_dir,
                    "pipeline_info/complete_samplesheet.valid.csv"),
             sep = ",", header = TRUE)
species <- metadat$species
taxon <- metadat$taxonomy

###############################################################################
# Get the file paths for all gene family trees, extract the gene
# family ID, and set as variable the base directory of all
# AA statistics
message("Done. Staging all gene family trees and AA statistics...")
gfts <-
  list.files(paste0(noveltree_out_dir, "generax/per_species_rates/"),
             pattern = "_reconciled_gft.newick",
             recursive = TRUE, full.names = TRUE)
gene_families <-
  gsub(".*\\/", "", gfts) |>
  gsub(pattern = "_reconciled_gft.newick", replacement = "")
aa_stat_basedir <-
  "aa-summary-stats/per-family-summaries/aa-physical-properties/"

# Identify which gene families include human proteins:
og_counts <- get_per_spp_og_counts(results_dir = noveltree_out_dir)
human_ogs <- og_counts$Orthogroup[which(og_counts["Homo-sapiens"] > 0)]

# And subset down based on this:
gfts <- gfts[which(gene_families %in% human_ogs)]
gfts <-
  setNames(gfts, gsub("*_reconciled_gft.newick", "", gfts) |>
             gsub(pattern = ".*/", replacement = ""))
gene_families <- gene_families[which(gene_families %in% human_ogs)]
gene_families <- setNames(gene_families, gene_families)

# Combine these into a list:
gene_families <-
  mapply(function(family, gft) c(family = family, gft = gft),
         gene_families, gfts[names(gene_families)],
         SIMPLIFY = FALSE)

message(paste0("Done. Carrying out conservation analysis of ",
               length(gene_families), " gene families using ",
               detectCores(), " cores."))
message("Hang in there, this will take a bit...")

# Now, calculate multivariate distances of each protein, and each species
# from humans and their respective protein gene copies
gf_conservation_res <-
  mclapply(gene_families, genefam_aa_conservation,
           spp_tree = spp_tree, ref_spp = "Homo-sapiens",
           aa_stat_basedir = aa_stat_basedir,
           mc.cores = detectCores(), mc.preschedule = FALSE,
           max_treepl_treesize = 9999)

message("Done!")
