library(elo)
library(ape)
library(phytools)
library(ggtree)
library(scales)
library(dplyr)
library(treeio)
library(ggnewscale)
library(phylobase)
library(phylosignal)
library(here)
library(arcadiathemeR)

# Function get gene symbols for uniprot IDs (relies on the function 'queryMany'
# from the package 'mygene')
get_gene_symbols <- function(uniprot) {
  # Convert to gene symbols
  genes <- mygene::queryMany(uniprot)

  # Add rownames
  rownames(genes) <- genes$query

  # Filter to just unique genes
  genes <- genes[unique(genes$query), ]

  # Match to input list
  genes <- genes[match(uniprot, genes$query), ]

  # Return
  return(genes)
}

# Function to darken color
darken_color <- function(color, factor = 1.4) {
  col <- col2rgb(color)
  col <- col / factor
  col <- rgb(t(col), maxColorValue = 255)
  col
}

#Function for extracting and reformatting tree data for ggtree
#(pulled from https://dmnfarrell.github.io/r/ggplottree)
gettreedata <- function(tree, meta) {
  #get treedata object
  d <- meta[row.names(meta) %in% tree$tip.label, ]
  d$label <- row.names(d)
  y <- full_join(as_tibble(tree), d, by = "label")
  y <- as.treedata(y)
  return(y)
}


# Set up color dictionary for plots
all_colors <- c(
  "#5088C5", "#F28360", "#F7B846", "#97CD78",
  "#7A77AB", "#f898AE", "#3B9886", "#c85152",
  "#73B5E3", "#BAB0A8", "#8A99AD", "#FFB984",
  "#C6E7F4", "#F8C5C1", "#F5E4BE", "#B5BEA4",
  "#DCBFFC", "#B6C8D4", "#DAD3C7", "#DA9085",
  "#2B65A1", "#094468", "#9E3F41", "#FFD0B0",
  "#FFD364", "#D68D22", "#A85E28", "#DCDFEF",
  "#54448C", "#C3E2DB", "#6FBCAD", "#2A6B5E",
  "#09473E", "#FFE3D4", "#E2718F", "#C14C70",
  "#EDE6DA", "#635C5A", "#E6EAED", "#CAD4DB",
  "#ABBAC4", "#687787", "#47784A", "#C1E1AE",
  "#71AC5A", "#1E4812"
)
