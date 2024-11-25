setwd('~/Documents/Research/github/2024-organismal-selection/')
source('code/utils.R')

###################
#####Load data#####
###################
##Phylogeny
# Load phylogeny
tree <- read.newick("data/congruified_spprax_species_tree.newick")

#Calculate cophenetic distance (will be used for downstream analyses)
cophen = cophenetic.phylo(tree)

#Extract human row (used as a gauge of distance from humans per-species)
cophen = cophen[grep('Homo-sapiens',
                     rownames(cophen)),]

# Create tree without human (will be used for downstream analyses)
tree_nh = drop.tip(tree, 'Homo-sapiens')

##Metadata
#Load
metadata = read.csv('data/metadata/metadata.csv')
#Match colors with taxa
taxa_colors = all_colors[1:length(unique(metadata$taxogroup2_unieuk))]
names(taxa_colors) = unique(metadata$taxogroup2_unieuk)

#Match with species
species_colors = taxa_colors[match(metadata$taxogroup2_unieuk,
                                   names(taxa_colors))]
names(species_colors) = metadata$genus.species

##Conservation
#Load
conservation_table = 
  readRDS('data/conservation_table.RDS')
conservation_species = 
  readRDS('data/conservation_by_spp.RDS')
conservation_species_best = 
  readRDS('data/conservation_best_protein_per_spp.RDS')

#Split conservation_table by gene family
conservation = split(conservation_table,
                     conservation_table$gene_family)

###############################################################################
#####Plot phylogeny ###########################################################
#####(Related to figure 2; organism icons added later via illustrator)#########
###############################################################################
#Create metadata for tree (colors with be based on 'taxogroup2...')
meta = data.frame(name = tree$tip.label,
                  taxa = metadata$taxogroup2_unieuk[
                    match(gsub('_', '-', tree$tip.label),
                          metadata$genus.species)],
                  row.names = tree$tip.label)

#Combine with tree via 'gettreedata'
tree_annot = gettreedata(tree,
                         meta)

#Create initial plot
p <- ggtree(tree_annot, 
            layout = 'fan', 
            open.angle = 45)

cols=c('taxa')
col <- cols[1]

#Reorder metadata by taxa
df<-meta[tree$tip.label,][col]

#Match colors to tree order
colors <- species_colors[match(tree$tip.label,
                               names(species_colors))]
#Use taxa as color names
names(colors) = df$taxa

#Plot all
p + new_scale_fill() +
  geom_tippoint(mapping=aes(fill=.data[[col]]),
                size=4,
                shape=21,
                stroke=0) +
  geom_tiplab(offset = 0.5) + 
  scale_fill_manual(values=colors, 
                    na.value="white")

####################################################################
#####Plot phylogeny with n gene families per species################
#####(Related to figure 3)##########################################
####################################################################
#Calculate n ogs per species
n_ogs_per_species = unlist(lapply(conservation_species,
                                  function(x) 
                                    length(unique((x$gene_family)))))

#Combine tree and gene family number via 'phylo4d'
p4d = phylo4d(tree,
              unique_species)

#Match taxa colors to tree
cols = species_colors[match(tree$tip.label,
                            names(species_colors))]

#Plot using 'barplot.phylo4d' from the 'phylosignal' package
barplot.phylo4d(p4d,
                center = FALSE,
                scale = FALSE, 
                bar.col = cols, 
                tree.ladderize = TRUE,
                trait.labels = 'n ogs', 
                use.edge.length = FALSE)

##Compare gene family n of vertebrates to all other species
#Vertebrate gene family n
verts = n_ogs_per_species[names(n_ogs_per_species)%in%
                            metadata$genus.species[
                              metadata$taxogroup2_unieuk == 'Vertebrata']]

#Non-vertebrate gene family n
non_verts = n_ogs_per_species[!names(n_ogs_per_species)%in%
                                metadata$genus.species[
                                  metadata$taxogroup2_unieuk == 'Vertebrata']]

#Compare with Kruskal-Wallis test
verts_ks = kruskal.test(list(verts, 
                             non_verts))
verts_ks

#####################################################################
#####Compare gene family n with phylogenetic distance################ 
#####(Related to figure 4)###########################################
#####################################################################
##Figure 4A: Plot cophenetic distance and gene family number
#Regress gene family number and phylogenetic distance
plot(cophen,
     n_ogs_per_species[match(names(cophen),
                             names(n_ogs_per_species))],
     ylab = 'n gene families',
     xlab = 'Distance from humans',
     cex.lab = 1.5,
     cex.axis = 1.5,
     pch = 20,
     cex = 2,
     col = species_colors[match(names(cophen),
                                names(species_colors))])

#Add species names
text(cophen,
     n_ogs_per_species[match(names(cophen),
                             names(n_ogs_per_species))],
     names(cophen))

#Add regression line
mod = lm(n_ogs_per_species[match(names(cophen),
                                 names(n_ogs_per_species))]~cophen)

abline(mod,
       lty = 'dashed')

##Figure 4B: 
##Compare the observed gene family distribution to expected using permutations
#Calculate n proteins per gene family
n_proteins_per_og = unlist(lapply(conservation,
                                  function(x) 
                                    nrow(x)))

#Calculate n species per OG
n_species_per_og = unlist(lapply(conservation,
                                 function(x) 
                                   length(unique((x$nonref_species)))))

#Extract fitted values from regression comparing gene family n and distance
fitted = mod$fitted.values

#Normalize by sum of all fitted values; converts to a proportion
fitted = fitted/sum(fitted)

#Create n observations per species (i.e. proportion of 10000 observations
#determined by 'fitted'); each permutation will randomly sample from 
#this distribution
n = 100000*fitted

pool = c()
for(i in 1:length(n)){
  pool = c(pool, rep(names(n)[i], n[i]))
}

#Create empty list for saving permutation results
perms = list()

#Counter
pb <- txtProgressBar(
  min = 1,
  max = length(seq(min(n_proteins_per_og),
                   max(n_proteins_per_og), 5)),
  style = 3,
  width = 100,
  char = "."
)

#Set seed
set.seed(666)

#Set up n permutations
n_perms = 100

#Loop over range encompassing the minimum and maximum number of proteins
#per gene family
for(i in seq(min(n_proteins_per_og),
             max(n_proteins_per_og), 5)){
  
  #Update counter
  setTxtProgressBar(pb, i)
  
  #Create empty vector to save individual permutations
  lengths = c()
  
  #Permutate
  for(j in 1:n_perms){
    
    #Calculate the number of species present per sample
    lengths = c(lengths, 
                length(table(pool[sample(1:length(pool), i)])))
    
  }
  
  #Add to results list
  perms[[as.character(i)]] = lengths
}

#Compare saturation point of permuted and real data 
#(i.e. which gene family size contains proteins from all 63 species)
perm_saturation = unlist(lapply(perms, 
                                function(x) length(grep(63, x))))
obs_saturation = split(n_species_per_og, 
                       n_proteins_per_og)
obs_saturation = unlist(lapply(obs_saturation, 
                               function(x) length(grep(63, x))))

#Filter to gene families with at least permutation containing all species
perm_saturation = perm_saturation[perm_saturation>0]
obs_saturation = obs_saturation[obs_saturation>0]

#Find smallest gene family size that contains all species
min(as.numeric(names(perm_saturation)))
#264
min(as.numeric(names(obs_saturation)))
#70

#Compare n proteins per og and n species
smoothScatter(log(n_proteins_per_og), 
              n_species_per_og, 
              nrpoints = 0,
              xlab = 'n proteins per gene family (log)',
              ylab = 'n species per gene family',
              cex.lab = 1.5,
              cex.axis = 1.5,
              bty = 'n',
              nbin = 100,
              colramp = colorRampPalette(c('white', 
                                           rev(
                                             unlist(
                                               arcadia_magma$color_dict)))))
lines(log(as.numeric(names(perms))), 
      unlist(lapply(perms, function(x) mean(x))),
      lwd = 1.5)

##Figure 4C: compare the number of proteins within each gene family
##to evolutionary distance
#Calculate phylogenetic dispersion per OG
phylo_dist_per_og = list()
for(i in 1:length(conservation)){
  
  #Get species in OG
  species_og = unique(conservation[[i]]$nonref_species)
  
  #Get MRCA
  mrca_og = getMRCA(tree, 
                    species_og[species_og%in%tree$tip.label])
  
  #Subset tree
  tree_og = extract.clade(tree,
                          mrca_og)
  
  #Calculate "age" (max branch length)
  phylo_dist_per_og[[names(conservation)[i]]] = max(tree_og$edge.length)
  
}
phylo_dist_per_og = unlist(phylo_dist_per_og)

#Plot
smoothScatter(n_species_per_og, 
              phylo_dist_per_og, 
              nrpoints = 0,
              xlab = 'n species per OG',
              ylab = 'Time span per OG (billion years)',
              cex.lab = 1.5,
              cex.axis = 1.5,
              bty = 'n',
              nbin = 100,
              colramp = colorRampPalette(c('white', 
                                           rev(
                                             unlist(
                                               arcadia_magma$color_dict)))))

#########################################################################
#####Plot molecular conservation distributions for all gene families#####
#####(Related to figure 5)###############################################
#########################################################################
#Calculate distributions of molecular conservation (via histograms)
#across all gene families
conservation_hists = lapply(conservation_species_best,
                            function(x) hist(x$trait_dist, 
                                             plot = FALSE, 
                                             breaks = seq(0, 50, 1))$density)

#Normalize histograms by maximum value (to make comparable across all families)
conservation_hists = lapply(conservation_hists, 
                            function(x) x/max(x))

#Combine
conservation_hists = do.call(rbind,
                             conservation_hists)

#Hierarchically cluster gene family distributions
hcl = hclust(dist(conservation_hists))

#Plot
plot(hcl, labels = FALSE)

##Plot gene family conservation distributions as points
#Calculate median conservation (for initiating plot)
stats = lapply(conservation_species_best, 
               function(x) median(x$trait_dist))
stats = stats[hcl$order]

#Order gene families based on hierarchical clustering
hcl_conservation = conservation_species_best[hcl$order]

#Plot
par(bg=NA)
plot(unlist(lapply(stats, function(x) x$stats[3])),
     ylim = c(0, 50),
     pch = 20,
     ylab = 'Distance from human',
     xlab = '',
     cex.axis = 1.5,
     cex.lab = 1.5)
for(i in 1:length(stats)){
  points(rep(i, length(x[[i]]$trait_dist)),
         hcl_conservation[[i]]$trait_dist,
         pch = 20,
         cex = 2,
         col = scales::alpha('black', 0.1))
}

####################################
#####Plot example gene families#####
#####(Related to figure 6)##########
####################################
#Calculate slopes (will be used to categorize gene families)
pb <- txtProgressBar(
  min = 1,
  max = length(conservation_species_best),
  style = 3,
  width = 100,
  char = "."
)

og_slopes = list()
for(i in 1:length(conservation_species_best)){
  
  #Update counter
  setTxtProgressBar(pb, i)
  
  #Split on species
  s = split(conservation_species_best[[i]],
            conservation_species_best[[i]]$nonref_species)
  
  #Select most conserved protein per species
  s = lapply(s, function(x) x[which.min(x$trait_dist),])
  
  #Recombine
  s = do.call(rbind, s)
  
  #Get species phylo distances
  dists = cophen[match(s$nonref_species,
                       names(cophen))]
  
  #Regression
  mod = lm(s$trait_dist~dists)
  slope = coef(mod)[2]
  r2 = summary(mod)[8]
  
  #Add to list
  og_slopes[[names(conservation_species_best)[i]]] = list(data = s,
                                                          r2 = r2,
                                                          slope = slope)
}

#Extract model fit, slope, and conservation per gene family
r2s = unlist(lapply(og_slopes, function(x) x$r2))[
  unlist(lapply(og_slopes, function(x) nrow(x$data))>40)
]

slopes = unlist(lapply(og_slopes, function(x) x$slope))[
  unlist(lapply(og_slopes, function(x) nrow(x$data))>40)
] 

cons = unlist(lapply(og_slopes, function(x) mean(x$data$trait_dist)))[
  unlist(lapply(og_slopes, function(x) nrow(x$data))>40)
] 

cons_var = unlist(lapply(og_slopes, function(x) 
  sd(x$data$trait_dist)/mean(x$data$trait_dist)))[
  unlist(lapply(og_slopes, function(x) nrow(x$data))>40)
] 

#Combine in a data frame for filtering
res = data.frame(protein = names(og_slopes)
                 [unlist(lapply(og_slopes, function(x) nrow(x$data))>40)],
                 r2 = r2s,
                 cons = cons,
                 cons_var = cons_var,
                 slope = slopes)

#Plot examples
par(mfrow = c(2,2))
#Positive relationship
x = res[res$slope>0,]
x = x[x$cons<10,]
x = x[which.max(x$r2),]
n = x$protein

dists = cophen[
  match(og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
        names(cophen))]

plot(dists,
     og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist,
     xlab = 'Cophenetic distance',
     ylab = 'Conservation with humans',
     cex.axis = 1.5,
     cex.lab = 1.5, 
     pch = 20,
     cex = 2,
     ylim = c(0, 3),
     col = species_colors[
       match(og_slopes[[grep(n, 
                             names(og_slopes))]]$data$nonref_species,
             names(species_colors))])
abline(lm(og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist~dists),
       lty = 'dashed',
       lwd = 1.5)
text(3000,
     3,
     paste('r2 =', signif(x$r2, 2), sep = ''))
title(main = paste(n, '(PTN4)'),
      font.main = 1,
      cex.main = 1.5)

#Positive relationship, human specific
x = res[res$slope>0,]
x = x[which.max(x$r2),]
n = x$protein

dists = cophen[match(og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
                     names(cophen))]

plot(dists,
     og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist,
     xlab = 'Cophenetic distance',
     ylab = 'Conservation with humans',
     ylim = c(0, 20),
     cex.axis = 1.5,
     cex.lab = 1.5, 
     pch = 20,
     cex = 2,
     col = species_colors[
       match(og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
             names(species_colors))])

abline(lm(og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist~dists),
       lty = 'dashed',
       lwd = 1.5)
text(3000,
     20,
     paste('r2 =', signif(x$r2, 2), sep = ''))
title(main = paste(n, '(FOXA1)'),
      font.main = 1,
      cex.main = 1.5)

#Deep conservation
x = res[res$cons<1,]
x = x[which.min(x$cons_var),]
n = x$protein

dists = cophen[match(og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
                     names(cophen))]

plot(dists,
     og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist,
     xlab = 'Cophenetic distance',
     ylab = 'Conservation with humans',
     cex.axis = 1.5,
     cex.lab = 1.5, 
     pch = 20,
     ylim = c(0, 2),
     cex = 2,
     col = species_colors[
       match(og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
             names(species_colors))])
abline(lm(og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist~dists),
       lty = 'dashed',
       lwd = 1.5)
text(3000,
     2,
     paste('r2 =', signif(x$r2, 2), sep = ''))
title(main = paste(n, '(ARF3)'),
      font.main = 1,
      cex.main = 1.5)

#Negative relationship
x = res[res$slope<0,]
x = x[which.max(x$r2),]
n = x$protein

dists = cophen[match(og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
                     names(cophen))]

plot(dists,
     og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist,
     xlab = 'Cophenetic distance',
     ylab = 'Conservation with humans',
     cex.axis = 1.5,
     cex.lab = 1.5, 
     pch = 20,
     cex = 2,
     ylim = c(0, 5),
     col = species_colors[
       match(og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
             names(species_colors))])

abline(lm(og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist~dists),
       lty = 'dashed',
       lwd = 1.5)
text(3000,
     5,
     paste('r2 =', signif(x$r2, 2), sep = ''))
title(main = paste(n, '(3HIDH)'),
      font.main = 1,
      cex.main = 1.5)







