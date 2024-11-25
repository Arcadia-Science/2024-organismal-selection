setwd('~/Documents/Research/github/2024-organismal-selection/')
source('code/utils.R')

#TO DO: For Reactome pathways, estimate age of conservation by identifying where conservation drops off across species; can pair functions with evolutionary scale

###################
#####Load data#####
###################
#####Species complexity#####
#Load
complexity = read.delim('data/metadata/complexity.txt')

#####Tree#####
# Load tree
tree <- read.newick("data/congruified_spprax_species_tree.newick")

#Remove water bear (incorrectly labeled; is actually a duplicate of Chlorella)
tree = drop.tip(tree,
                'Hypsibius-dujardini')

#####Metadata#####
#Load
metadata = read.csv('data/metadata/metadata.csv')[-66,]

#Match colors with taxa
taxa_colors = all_colors[1:length(unique(metadata$taxogroup2_unieuk))]
names(taxa_colors) = unique(metadata$taxogroup2_unieuk)

#Match with species
species_colors = taxa_colors[match(metadata$taxogroup2_unieuk,
                                   names(taxa_colors))]
names(species_colors) = metadata$genus.species

#####Protein conservation scores#####
#Set working directory (will need to replace this once data are deposited)
setwd('~/Documents/Research/arcadia-projects/raas-organism-prioritization/conservation_score_v2_06012024/gf-aa-multivar-distances/final_protein_pair_summary_tables/')

#List files
files = list.files()

#Load all pairs
conservation = list()
for(i in 1:length(files)){
  
  #Load file and add to list
  conservation[[gsub('_final_summary_table.tsv', 
                      '', 
                      files[i])]] = read.delim(files[i])
  
  #Add normalize trait distance column
  conservation[[i]]$rank_trait_dist_norm = 
    conservation[[i]]$rank_trait_dist/max(conservation[[i]]$rank_trait_dist)
}

#Remove Hypsibius-dujardini
#conservation = lapply(conservation, function(x) x[-grep('Hypsibius-dujardini',
#                                                        x$nonref_species),])

#Combine into a table
conservation_table = do.call(rbind, conservation)

#Add complexity to the table
conservation_table$complexity = complexity$complexity[
  match(conservation_table$nonref_species,
        complexity$species)]

#Split by gene family
conservation_split = split(conservation_table,
                           conservation_table$gene_family)

#Split protein pair table by species
conservation_species = split(conservation_table, 
                             conservation_table$nonref_species)

#Create filtered protein pair table with complexity
conservation_species_best = do.call(rbind, conservation_split)
conservation_species_best = split(conservation_species_best, 
                                   conservation_species_best$ref_protein)
for(i in 1:length(conservation_species_best)){
  
  #Select most conserved protein per species
  x = split(conservation_species_best[[i]], 
            conservation_species_best[[i]]$nonref_species)
  x = do.call(rbind, lapply(x, function(y) y[which.min(y$trait_dist),]))
  
  #Add complexity
  x$complexity = complexity$complexity[match(x$nonref_species,
                                             complexity$species)]
  
  #Add to list
  conservation_species_best[[names(conservation)[i]]] = x
}

####################################
#####Plot example gene families#####
####################################
#Calculate slopes
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

r2s = unlist(lapply(og_slopes, function(x) x$r2))[
  unlist(lapply(og_slopes, function(x) nrow(x$data))>40)
]

slopes = unlist(lapply(og_slopes, function(x) x$slope))[
  unlist(lapply(og_slopes, function(x) nrow(x$data))>40)
] 

cons = unlist(lapply(og_slopes, function(x) mean(x$data$trait_dist)))[
  unlist(lapply(og_slopes, function(x) nrow(x$data))>40)
] 

cons_var = unlist(lapply(og_slopes, function(x) sd(x$data$trait_dist)/mean(x$data$trait_dist)))[
  unlist(lapply(og_slopes, function(x) nrow(x$data))>40)
] 

res = data.frame(protein = names(og_slopes)
                 [unlist(lapply(og_slopes, function(x) nrow(x$data))>40)],
                 r2 = r2s,
                 cons = cons,
                 cons_var = cons_var,
                 slope = slopes)

par(mfrow = c(2,2))
#Positive relationship
x = res[res$slope>0,]
x = x[x$cons<10,]
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
     ylim = c(0, 3),
     col = species_colors[match(og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
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
     col = species_colors[match(og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
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
#dists = cophen[match(og_slopes$A0A024RBG1$data$nonref_species,
#                     names(cophen))]

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
     col = species_colors[match(og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
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
     col = species_colors[match(og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
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



plot(dists,
     og_slopes[[grep(n, names(og_slopes))]]$data$trait_dist,
     xlab = 'Cophenetic distance',
     ylab = 'Conservation with humans',
     cex.axis = 1.5,
     cex.lab = 1.5, 
     pch = 20,
     cex = 2,
     ylim = c(0, 3),
     col = species_colors[match(og_slopes[[grep(n, names(og_slopes))]]$data$nonref_species,
                                names(species_colors))])

#####################################################
#####Plot tree and descriptive gene family stats#####
#####################################################
#Create metadata for tree (taxa IDs)
meta = data.frame(name = tree$tip.label,
                  taxa = metadata$taxogroup2_unieuk[match(gsub('_', '-', tree$tip.label),
                                                          metadata$genus.species)],
                  row.names = tree$tip.label)

#Combine with tree via 'gettreedata'
y = gettreedata(tree,
                meta)

#Create initial plot
p <- ggtree(y, 
            layout = 'fan', 
            open.angle = 45)

cols=c('taxa')
col <- cols[1]

#Reorder
df<-meta[tree$tip.label,][col]

#Get colors
colors <- species_colors[match(tree$tip.label,
                               names(species_colors))]
names(colors) = df$taxa

#Plot
p + new_scale_fill() +
  geom_tippoint(mapping=aes(fill=.data[[col]]),
                size=4,
                shape=21,
                stroke=0) +
  geom_tiplab(offset = 0.5) + 
  scale_fill_manual(values=colors, 
                    na.value="white")

###################################
#####Example gene family trees#####
###################################
#Load
setwd('../congruified-gfts/')
files = list.files()
par(mfrow = c(6,12))
for(i in seq(900, 8000, 100)){
  gtf = read.tree(files[i])
  
  #Remove Hypsibius-dujardini
  gtf = drop.tip(gtf, gtf$tip.label[grep('Hypsibius-dujardini', gtf$tip.label)])
  
  species = gtf$tip.label
  species = unlist(lapply(lapply(species, 
                                 function(x) strsplit(x, '_')), 
                          function(y) y[[1]][1]))
  gtf$tip.label = species
  print(max(gtf$edge.length))

  cols = species_colors[match(species,
                              names(species_colors))]
  clades = match(cols, unique(cols))
  clades = split(species, clades)
  
  ecol = phyloch::edge.color(gtf, 
                             groups = clades,
                             col = unique(cols))
  plot.phylo(gtf, 
             edge.color = ecol, 
             show.tip.label = FALSE,
             type = 'unrooted',
             no.margin = TRUE,
             edge.width = 1.5)
}


##################################
#####Distribution of OG sizes#####
##################################
#Calculate n proteins per OG
n_proteins_per_og = unlist(lapply(conservation,
                                  function(x) nrow(x)))

#Calculate n species per OG
n_species_per_og = unlist(lapply(conservation,
                                 function(x) length(unique((x$nonref_species)))))

#Calculate mean number of ogs per species
n_ogs_per_species = unlist(lapply(conservation_species,
                                  function(x) length(unique((x$gene_family)))))

#Calculate n ogs per species
unique_species = lapply(conservation,
                        function(x) unique(x$nonref_species))
unique_species = table(unlist(unique_species))
unique_species = c(unique_species, 
                   length(conservation))
names(unique_species)[length(unique_species)] = 'Homo-sapiens'
unique_species = unique_species[match(tree$tip.label,
                                      names(unique_species))]

#Compare vertebrates and non-verts
v = n_ogs_per_species[names(n_ogs_per_species)%in%
                        metadata$genus.species[metadata$taxogroup2_unieuk == 'Vertebrata']]
n = n_ogs_per_species[!names(n_ogs_per_species)%in%
                        metadata$genus.species[metadata$taxogroup2_unieuk == 'Vertebrata']]


kruskal.test(list(v, n))

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

#Calculate cophenetic distance
cophen = cophenetic(tree)

#Extract human row
cophen = cophen[grep('Homo-sapiens',
                     rownames(cophen)),]

#Calculate null
mod = lm(unique_species[match(names(cophen),
                              names(unique_species))]~cophen)
fitted = mod$fitted.values
fitted = fitted/sum(fitted)

n = 100000*fitted
pool = c()
for(i in 1:length(n)){
  pool = c(pool, rep(names(n)[i], n[i]))
}

perms = list()
pb <- txtProgressBar(
  min = 1,
  max = length(seq(min(n_proteins_per_og),
                     max(n_proteins_per_og), 5)),
  style = 3,
  width = 100,
  char = "."
)
for(i in seq(min(n_proteins_per_og),
    max(n_proteins_per_og), 5)){
  setTxtProgressBar(pb, i)
  lengths = c()
  for(j in 1:100){
    lengths = c(lengths, 
                length(table(pool[sample(1:length(pool), i)])))
    
  }
  perms[[as.character(i)]] = lengths
}

#Compare saturation point of permuted and real data
x = unlist(lapply(perms, function(x) length(grep(64, x))))
y = split(n_species_per_og, n_proteins_per_og)
y = unlist(lapply(y, function(x) length(grep(64, x))))

x = x[x>0]
y = y[y>0]

min(as.numeric(names(x)))
#234
min(as.numeric(names(y)))
#71

#Compare n proteins per og and n species
par(mfrow = c(1,2))
smoothScatter(log(n_proteins_per_og), 
              n_species_per_og, 
              nrpoints = 0,
              xlab = 'n proteins per OG (log)',
              ylab = 'n species per OG',
              cex.lab = 1.5,
              cex.axis = 1.5,
              bty = 'n',
              nbin = 100,
              colramp = colorRampPalette(c('white', 
                                                  rev(unlist(arcadia_magma$color_dict)))))
lines(log(as.numeric(names(perms))), 
      unlist(lapply(perms, function(x) mean(x))),
      lwd = 1.5)

#Compare n proteins per og and evolutionary distance
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
                                                  rev(unlist(arcadia_magma$color_dict)))))


OG0013524 


plotTree.barplot(tree,
                 unique_species)

p4d = phylobase::phylo4d(tree,
                            unique_species)
cols = species_colors[match(tree$tip.label,
                            names(species_colors))]
phylosignal::barplot.phylo4d(p4d,
                             center = FALSE,
                             scale = FALSE, 
                             bar.col = cols, 
                             tree.ladderize = TRUE,
                             trait.labels = 'n ogs', 
                             use.edge.length = FALSE)

#Compare to number of gene families
plot(cophen,
     unique_species[match(names(cophen),
                          names(unique_species))],
     ylab = 'n ogs',
     xlab = 'Distance from humans',
     cex.lab = 1.5,
     cex.axis = 1.5,
     pch = 20,
     cex = 2,
     col = species_colors[match(names(cophen),
                                names(species_colors))])
text(cophen,
     unique_species[match(names(cophen),
                          names(unique_species))],
     names(cophen))

summary(lm(unique_species[match(names(cophen),
                                names(unique_species))]~cophen))

#PGLS
x = unique_species[match(names(cophen),
                         names(unique_species))]

#Create data frame
datos = data.frame(n_ogs = x,
                   Species = names(x))

#Create comparative data object
comp.data<-caper::comparative.data(tree, 
                                   datos, 
                                   names.col="Species", 
                                   vcv.dim=2, 
                                   warn.dropped=TRUE)

#PGLS
mod <- caper::pgls(n_ogs~1, 
                   data=comp.data,
                   lambda = 'ML')

#Plot phylogenetic residuals
res = as.numeric(mod$residuals)
names(res) = rownames(mod$residuals)
plot(mod$res,
     ylim = c(-120, 120))

################################################
#####Plot example gene family aa statistics#####
################################################
#Load
aa_stats_phylo = t(read.delim('data/aa-summary-statistics/OG0000011_phylo_corr_dat.tsv',
                              row.names = 1))
aa_stats = t(read.csv('data/aa-summary-statistics/OG0000011_summary_statistics.csv',
                       row.names = 1))
aa_stats = aa_stats[,match(colnames(aa_stats_phylo),
                           colnames(aa_stats))]
aa_stats = aa_stats[rownames(aa_stats)%in%rownames(aa_stats_phylo),]
aa_dist = t(read.delim('data/aa-summary-statistics/OG0000011_protein_dists.tsv'))
aa_tree = read.tree('data/aa-summary-statistics/OG0000011_congruified.newick')

#Filter tree
aa_tree = keep.tip(aa_tree, colnames(aa_stats)[1:30])

#Reorder
aa_stats = aa_stats[,match(aa_tree$tip.label,
                           colnames(aa_stats))]
aa_stats_phylo = aa_stats_phylo[,match(aa_tree$tip.label,
                                 colnames(aa_stats_phylo))]

#Plot
image(t(cbind(scale(aa_stats[1,1:30]),
            scale(aa_stats_phylo[1,1:30]))),
      col = colorRampPalette(unlist(arcadia_magma$color_dict))(100))

#Plot tree
plot(aa_tree)

#Plot other measures
image(t(as.matrix(scale(aa_stats_phylo[2:10,1:30]))),
      col = colorRampPalette(unlist(arcadia_magma$color_dict))(100))

#NMDS
aa_mds = cmdscale(aa_dist[1:30, 1:30])
plot(aa_mds,
     pch = 20,
     cex = 2,
     xlab = 'MDS 1',
     ylab = 'MDS 2',
     cex.axis = 1.5,
     cex.lab = 1.5,
     col = taxa_colors[order(names(taxa_colors))][4:length(taxa_colors)])

#Distance distribution
x = aa_dist[grep('Aedes.aegypti_tr.A0A1S4FXZ9.A0A1S4FXZ9.AEDAE',
                 rownames(aa_dist)),1:30]

plot(x[order(x)],
     pch = 20,
     cex = 2,
     col = taxa_colors[order(names(taxa_colors))][4:length(taxa_colors)][order(x)],
     ylab = 'Distance',
     xlab = 'Protein',
     cex.axis = 1.5,
     cex.lab = 1.5, 
     las = 2,
     xaxt = 'n')

##################################################################
#####Distribution of conservation scores across gene families#####
##################################################################
#Calculate histograms across all gene families
conservation_hists = lapply(conservation_species_best,
                            function(x) hist(x$trait_dist, plot = FALSE, 
                                             breaks = seq(0, 50, 1))$density)
conservation_hists = lapply(conservation_hists, 
                            function(x) x/max(x))
conservation_hists = do.call(rbind,
                             conservation_hists)

#Plot as heatmap
hcl = hclust(dist(conservation_hists))
plot(hcl, labels = FALSE)

conservation_hists = conservation_hists[hcl$order,]
image(conservation_hists,
      col = colorRampPalette(c('white', rev(unlist(arcadia_magma$color_dict))))(100),
      ylim = c(0, 0.5))

#Plot as points
stats = lapply(conservation_species_best, function(x) boxplot.stats(x$trait_dist))
stats = stats[hcl$order]
x = conservation_species_best[hcl$order]

png('~/Desktop/test.png',
    width = 3600,
    height = 1200)
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
         x[[i]]$trait_dist,
         pch = 20,
         cex = 2,
         col = scales::alpha('black', 0.1))
}
#points(unlist(lapply(stats, function(x) x$stats[3])),
#       pch = 21,
#       cex = 1,
#       bg = 'white',
#       col = 'black')
dev.off()


###################################################
#####Compare conserved species per gene family#####
###################################################
#Get top species per gene family
top_species = lapply(conservation, 
                     function(x) x$nonref_species[which.min(x$rank_trait_dist)])

#Sort
top_species = sort(table(unlist(top_species)))

#Calculate mean conservation ranking per species
rank_trait_dist = list()
for(i in 1:length(conservation_species)){
  
  #Get species name
  n = names(conservation_species)[i]
  
  #Split on gene family
  gene_fams = split(conservation_species[[i]],
                    conservation_species[[i]]$gene_family)
  
  #Extract best conservation ranking for each gene family
  rank_trait_dist[[n]] = unlist(lapply(gene_fams,
                                       function(x) 
                                         min(x$rank_trait_dist_norm)))
}

x = unlist(lapply(rank_trait_dist, function(x) mean(x)))
y = unlist(lapply(rank_trait_dist, function(x) length(x)))

######################################################################
#####Morphospace of gene family conservation values (genome-wide)#####
######################################################################
#Create matrix of conservation values for each gene family
conservation_matrix = as.data.frame(matrix(ncol = length(conservation_species_best),
                                           nrow = nrow(complexity)))

#Column and rownames
colnames(conservation_matrix) = names(conservation_species_best)
rownames(conservation_matrix) = complexity$species

#Add in values
for(i in 1:length(conservation_species_best)){
  x = conservation_species_best[[i]]$trait_dist[match(rownames(conservation_matrix),
                                                      conservation_species_best[[i]]$nonref_species)]
  conservation_matrix[,i] = x
}

#Replace NAs
conservation_matrix[is.na(conservation_matrix)] = max(conservation_table$trait_dist)

#PCA
pca = prcomp(conservation_matrix)

plot(pca$x)
text(pca$x,
     rownames(pca$x))

#UMAP
u = umap::umap(conservation_matrix,
               verbose = TRUE)

plot(u$layout)
text(u$layout,
     rownames(u$layout))

#Phylomorphospace
tree_nh = drop.tip(tree, 'Homo-sapiens')
x = pca$x[match(tree_nh$tip.label,
                rownames(pca$x)),1:2]
#x = u$layout[match(tree_nh$tip.label,
#             rownames(u$layout)),1:2]
cols = species_colors[match(rownames(x),
                            names(species_colors))]
cols = c(cols, rep('grey90', tree_nh$Nnode))
names(cols) = 1:(length(tree_nh$tip)+tree_nh$Nnode)
phylomorphospace(tree_nh,
                 x, 
                 ftype="off",
                 node.by.map=TRUE,
                 node.size = 1.5,
                 control=list(col.node=cols),
                 xlab = paste("PC1", "(45.78%)"),
                 ylab = paste("PC2", ("7.73%")),
                 cex.axis = 1.5,
                 cex.lab = 1.5)

#Correlation between distance and PCs
x = cophen[match(rownames(pca$x),
                 names(cophen))]
res = apply(pca$x, 2, function(y) cor.test(x, y))
plot(unlist(lapply(res, function(x) x$p.value)))

#Correlation between ortholog n and PC1
res = apply(pca$x, 2, function(y) cor.test(n_ogs_per_species, y[match(names(n_ogs_per_species),
                                                                      rownames(pca$x))]))
plot(unlist(lapply(res, function(x) x$p.value)))

#Branch length variation
euc.dist <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))

eucs = c()
for(i in 1:20){
  x = pca$x[match(tree_nh$tip.label,
                  rownames(pca$x)),i:(i+1)]
  out = phylomorphospace(tree_nh,
                         x, 
                         ftype="off",
                         node.by.map=TRUE,
                         node.size = 1.5)
  
  pos = rbind(out$xx, out$yy)
  eucs = c(eucs, 
           sd(apply(pos, 2, function(x) euc.dist(x[1], x[2])))/
             mean(apply(pos, 2, function(x) euc.dist(x[1], x[2]))))
}
plot(eucs)

#Make space using n gene families and average conservation for each species?
n_ogs = lapply(conservation_species,
               function(x) length(unique(x$gene_family)))
n_ogs = unlist(n_ogs)

mean_dist = lapply(conservation_species, 
                   function(x) mean(x$trait_dist))

plot(n_ogs, m[match(names(n_ogs), names(m))])

############################################
#####Per-species conservation Elo score#####
############################################
#Get matchups
matchups = expand.grid(unique(conservation_table$nonref_species),
                       unique(conservation_table$nonref_species))

#Organize outcomes based on most conserved proteins per gene family
matchup_outcomes = list()
pb <- txtProgressBar(
  min = 1,
  max = length(conservation_species_best),
  style = 3,
  width = 100,
  char = "."
)
for(i in 1:length(conservation_species_best)){
  
  # Update counter
  setTxtProgressBar(pb, i)
  
  #Get all combinations
  x = expand.grid(unique(conservation_species_best[[i]]$nonref_species),
                  unique(conservation_species_best[[i]]$nonref_species))
  
  #Convert to data frame
  x = data.frame(spp1 = as.character(x$Var1),
                 spp2 = as.character(x$Var2),
                 spp1_cons = rep(NA, nrow(x)),
                 spp2_cons = rep(NA, nrow(x)))
  
  #Remove self matches
  x = x[!x[,1] == x[,2],]
  
  #Add outcomes
  for(j in 1:nrow(x)){
    x$spp1_cons[j] = min(conservation_species_best[[i]]$rank_trait_dist_norm[
      grep(x$spp1[j],
           conservation_species_best[[i]]$nonref_species)])
    
    x$spp2_cons[j] = min(conservation_species_best[[i]]$rank_trait_dist_norm[
      grep(x$spp2[j],
           conservation_species_best[[i]]$nonref_species)])
  }
  
  #Add to list
  matchup_outcomes[[names(conservation_species_best)[i]]] = x

}

#Save
#saveRDS(matchup_outcomes, 
#        'out/elo_matchups.RDS')

#Get protein pairs
protein_sample = do.call(rbind, lapply(matchup_outcomes, 
                                       function(x) if(nrow(x)>= 10){
                                         x[sample(1:nrow(x),10),]}))

#Loop over n random combinations of proteins
n_perms = 50
all_elo_iterations = list()
all_species_elo_distributions = list()
set.seed(1234)
for(g in 1:n_perms){
  
  #Update counter
  print(paste(g, "out of", n_perms))
  
  #Calculate elo scores n times
  mean_elo_scores = list()
  
  #Subsample outcomes (right now this is from a single gene family, 
  #will need to develop a strategy to sample n matchups from each family)
  ##TO TRY: Choose n proteins from each species to make it an even number per
  test = protein_sample[sample(1:nrow(protein_sample), 10000),]
  
  for(h in 1:n_perms){
    
    #Reorder
    test = test[sample(1:nrow(test)),]
    
    #Remove self
    test = test[!test[,1] == test[,2],]
    
    #Create matrix containing elo for each species 
    #(to be updated w/ each match)
    unique_species = unique(test$spp1)
    elo_scores = as.data.frame(rep(1500,
                                   length(unique_species)),
                               row.names = unique_species)
    
    #Create species elo lists to keep score records
    species_elo = split(rep(1500, length(unique_species)), 
                        unique_species)
    
    #Loop over and simulate outcomes for each match
    pb <- txtProgressBar(
      min = 1,
      max = nrow(test),
      style = 3,
      width = 100,
      char = "."
    )
    
    for(i in 1:nrow(test)){
      
      # Update counter
      setTxtProgressBar(pb, i)
      
      #Get outcome
      outcome = unlist(test[i,3:4])
      names(outcome) = test[i,1:2]
      
      #Convert to wins
      if(outcome[1] == outcome[2]){
        outcome = c(0.5, 0.5)
      }else{
        x = which.min(outcome)
        y = which.max(outcome)
        outcome[x] = 1
        outcome[y] = 0
      }
      
      #Get probabilities
      probs = elo_scores[match(c(test[i,1:2]),
                               rownames(elo_scores)),1]
      
      #Calculate elo
      elo_update = elo.calc(outcome, 
                            rep(probs[1], 2), 
                            rep(probs[2], 2), 
                            k = 4)[1,]
      names(elo_update) = names(outcome)
      
      #Update elo matrix
      elo_scores[grep(names(elo_update)[1],
                      rownames(elo_scores)),1] = elo_update[1]
      
      elo_scores[grep(names(elo_update)[2],
                      rownames(elo_scores)),1] = elo_update[2]
      
      #Update species scores
      x = species_elo[[grep(names(elo_update)[1],
                            names(species_elo))]]
      x = unlist(c(x, elo_update[1]))
      species_elo[[grep(names(elo_update)[1],
                        names(species_elo))]] = x
      
      x = species_elo[[grep(names(elo_update)[2],
                            names(species_elo))]]
      x = unlist(c(x, elo_update[2]))
      species_elo[[grep(names(elo_update)[2],
                        names(species_elo))]] = x
    
    #Add to list
    mean_elo_scores[[as.character(h)]] = elo_scores
    all_species_elo_distributions[[as.character(h)]] = species_elo
  }
  
  #Calculate mean elo score per species
  n = as.character(unique(matchups$Var1))
  
  mean_elo_scores = do.call(cbind, 
                            lapply(mean_elo_scores, 
                                   function(x) x[match(n,
                                                       rownames(x)),1]))
  rownames(mean_elo_scores) = n
  
  #Calculate means and error
  final_elo = rowMeans(mean_elo_scores)
  lower = apply(mean_elo_scores, 1, function(x) quantile(x, prob = 0.05))
  upper = apply(mean_elo_scores, 1, function(x) quantile(x, prob = 0.95))
  
  o = order(final_elo, decreasing = TRUE)
  plot(final_elo[o],
       pch = 20,
       cex = 1.5,
       xaxt = 'n',
       ylab = 'Elo score',
       xlab = '',
       bty = 'n',
       ylim = c(min(lower), 
                max(upper)))
  abline(h = 1500,
         lty = 'dashed',
         lwd = 1.5)
  axis(1,
       1:length(final_elo), 
       labels = names(final_elo)[o],
       las = 2)
  segments(1:length(final_elo),
           lower[o],
           1:length(final_elo),
           upper[o])
  }
  #Add to list
  all_elo_iterations[[as.character(g)]] = mean_elo_scores
  
}

#Combine into matrix
all_elo_iterations = do.call(cbind, 
                             all_elo_iterations)

#Save
#saveRDS(all_elo_iterations, 'out/elo_permutations.RDS')

#Load
all_elo_iterations = readRDS('out/elo_permutations.RDS')

#Plot example outcomes
plot(approx(1:length(species_elo[[1]]),
                     species_elo[[1]], 
                     n = 100)$y,
     type = 'l',
     ylim = c(1400, 1650),
     xlim = c(0, 100),
     ylab = 'Elo score',
     xlab = 'Match number',
     cex.axis = 1.5,
     cex.lab = 1.5,
     col = species_colors[match(names(species_elo)[1],
                                names(species_colors))])
abline(h = 1500, lty = 'dashed', lwd = 1.5)
for(i in 1:length(species_elo)){
  lines(approx(1:length(species_elo[[i]]),
               species_elo[[i]], 
               n = 100)$y,
        col = species_colors[match(names(species_elo)[i],
                                   names(species_colors))])
}

final = unlist(lapply(species_elo, function(x) x[length(x)]))
plot(rep(1, length(final)),
     final,
     xlim = c(0.9, 1.1),
     ylim = c(1400, 1650),
     pch = 20,
     cex = 2,
     col = species_colors[match(names(species_elo),
                                names(species_colors))])

#Calculate cophenetic distance
cophen = cophenetic(tree)

#Extract human row
cophen = cophen[grep('Homo-sapiens',
                     rownames(cophen)),]

#Calculate mean elo per species
m = rowMeans(all_elo_iterations)

#Calculate mean standardized variance
m_sd = apply(all_elo_iterations, 1, function(x) sd(x)/mean(x))

#Plot against complexity
plot(cophen[match(names(m), names(cophen))],
     m,
     type = 'n',
     xlab = 'Cophenetic distance',
     ylab = 'Elo score',
     cex.axis = 1.5, 
     cex.lab = 1.5)
abline(h = 1500,
       lty = 'dashed',
       lwd = 1.5)
points(cophen[match(names(m), names(cophen))],
       m,
       pch = 21,
       bg = species_colors[match(names(m),
                                 names(species_colors))],
       col = darken_color(species_colors[match(names(m),
                                               names(species_colors))]),
       cex = 1.5)
text(cophen[match(names(m), names(cophen))],
     m,
     names(m))

#Violin plot
all_elo_iterations_list = split(all_elo_iterations,
                                rownames(all_elo_iterations))
all_elo_iterations_list = all_elo_iterations_list[order(unlist(lapply(all_elo_iterations_list, 
                                                               function(x)
                                                                 mean(x))))]

#Non vertebrate vs. vetebrate elo
n = metadata$taxogroup2_unieuk[match(names(m),
                                     metadata$genus.species)]
median(m[n == 'Vertebrata'])
#1571.394
median(m[!n == 'Vertebrata'])
#1478.744

#Relative probabilities
m = sort(m)
elo::elo.prob(m[length(m)], m[1])
elo::elo.prob(m[grep('Chlorella', names(m))], 
              m[grep('Chlamydomonas', names(m))])

#Compare to n ogs
plot(n_ogs_per_species,
     m[match(names(n_ogs_per_species),
             names(m))])

mod = lm(m[match(names(n_ogs_per_species),
                 names(m))]~
           n_ogs_per_species)

#Residuals
res = sort(MASS::studres(mod))

#Compare by phylogenetic distance to humans
p = rep(NA, length(m))
x = cophen[match(names(m),
                 names(cophen))]
p[x<1130] = 'young'
p[x>1130 & x<2000] = 'mid'
p[x>2000] = 'old'
p = split(m, p)

kruskal.test(list(p$young, c(p$mid, p$old)))
kruskal.test(list(p$mid, p$old))

p = rep(NA, length(m))
x = cophen[match(names(m),
                 names(cophen))]
p[x<3000] = 'old'
p[x<1700 & x>1200] = 'invert'
p = split(m, p)
kruskal.test(list(p$invert,
                  p$old))

#Compare by taxonomic group
n = metadata$taxogroup1_unieuk[match(names(res),
                                     metadata$genus.species)]
n = split(res, n)
sort(unlist(lapply(n, function(x) median(x))))

#########################################################
#####Network visualizations of gene family distances#####
#########################################################
#Load distance
og_dist = readRDS('~/Documents/Research/arcadia-projects/raas-organism-prioritization/conservation_score_v1/gf-aa-multivar-distances/dist-mats/OG0002302_distance_mat.RDS')

#PCA
pca = prcomp(og_dist$distance_matrix)
plot(pca$x[,1:2])

#Hclust
hcl = hclust(dist(og_dist$distance_matrix))
plot(hcl)

###############################################################
#####Per-species reactome pathway conservation with humans#####
###############################################################
#TO DO: Calculate PGLS on proportion complete/species for each reactome pathway, 
#analyze residuals as a way to rank species enrichments across pathways
#Load reactome
reactome = read.delim('~/Desktop/rationally-selecting-research-organisms/UniProt2Reactome.txt', header = FALSE)

#Filter to just human
reactome = reactome[grep('Homo sapiens', reactome[,6]),]

#Add column names
colnames(reactome) = c('uniprot', 'id', 'reactome_id', 'pathway', 'iea', 'species')

#Split on pathway
reactome = split(reactome, 
                 reactome[,4])

#Create human protein list for each species
conservation_human_proteins = lapply(conservation_species, function(x) unique(x$ref_protein))
conservation_human_proteins[['Homo-sapiens']] = unique(unlist(conservation_human_proteins))

#Calculate %conservation in pathways
pathway_conservation = list()
n = 100
for(i in 1:length(conservation_human_proteins)){
  print(names(conservation_human_proteins)[i])
  y = conservation_human_proteins[[i]]
  z = conservation_human_proteins$`Homo-sapiens`
  
  #Loop over reactome and compute statistic
  p_vals = c()
  eff_size = c()
  n_proteins = c()
  pb <- txtProgressBar(
    min = 1,
    max = length(reactome),
    style = 3,
    width = 100,
    char = "."
  )
  for(j in 1:length(reactome)){
    # Update counter
    setTxtProgressBar(pb, j)
    
    x = reactome[[j]]
    p = y[y%in%x$uniprot]
    hits = conservation_species[[i]][conservation_species[[i]]$ref_protein%in%p,]
    hits = unlist(lapply(split(hits$trait_dist, 
                               hits$ref_protein), 
                         function(x) x[which.min(x)]))
    perms = c()
    for(k in 1:n){
      perms = c(perms, mean(conservation_table$trait_dist
                            [sample(1:nrow(conservation_table),
                                    length(hits))]))
    }
    
   #Calculate cohen's d*
    if(length(hits)>3){
      eff_size = c(eff_size, effsize::cohen.d(hits, perms)$estimate)
    }else{
      eff_size = c(eff_size, NA)
    }
   
   #Calculate p-value
   p_vals = c(p_vals, 
              sum(mean(hits, na.rm = TRUE)<perms)/n)
   
   #n_proteins
   n_proteins = c(n_proteins, 
                  length(hits))
  }
 
  #Calculate percentage conservation for each pathway
  percent_conserved = unlist(lapply(reactome, 
                                    function(x) 
                                      sum(y%in%x$uniprot)/sum(z%in%x$uniprot)))
  
  print(mean(eff_size, na.rm = TRUE))
  
  #Add to results
  pathway_conservation[[names(conservation_human_proteins)[i]]] = list(p_value = p_vals,
                                                                       eff_size = eff_size,
                                                                       n_proteins = n_proteins,
                                                                       percent_conserved = percent_conserved)
}

#Save
saveRDS(pathway_conservation,
        'out/reactome_pathway_conservation.RDS')

#Plot relationship between percent conserved and effect size
plot(unlist(lapply(pathway_conservation, function(x) 
  mean(x$percent_conserved, na.rm = TRUE))),
     unlist(lapply(pathway_conservation, function(x) 
       mean(x$eff_size, na.rm = TRUE))),
  ylab = 'Mean conservation',
  xlab = 'Mean %conserved',
  cex.axis = 1.5,
  cex.lab = 1.5,
  pch = 20,
  col = species_colors[match(names(pathway_conservation),
                             names(species_colors))])
text(unlist(lapply(pathway_conservation, function(x) 
  mean(x$percent_conserved, na.rm = TRUE))),
  unlist(lapply(pathway_conservation, function(x) 
    mean(x$eff_size, na.rm = TRUE))),
  names(pathway_conservation))

#Convert to histogram frequencies
pathway_conservation_hist = lapply(pathway_conservation, 
                                   function(x) hist(x$eff_size, 
                                                    breaks = seq(-10, 10, 0.05))$counts)

#Normalize
pathway_conservation_hist = 
  lapply(pathway_conservation_hist, function(x) x/max(x))

#Combine
pathway_conservation_hist = do.call(rbind, pathway_conservation_hist)

#Reorder
pathway_conservation_hist = pathway_conservation_hist[match(tree_nh$tip.label,
                                                            rownames(pathway_conservation_hist)),]

#PCA
pca = prcomp(pathway_conservation_hist)

x = pca$x[match(tree_nh$tip.label,
          rownames(pca$x)),]
phylomorphospace(tree_nh, x[,1:2])

#Calculate cophenetic distance
cophen = ape::cophenetic.phylo(tree)[1,]
cophen = cophen[match(tree$tip.label,
                      names(cophen))]

plot(cophen,
     pathway_conservation_hist[,ncol(pathway_conservation_hist)])
text(cophen,
     pathway_conservation_hist[,ncol(pathway_conservation_hist)],
     tree$tip.label)

##################################################
#####Individual reactome pathway conservation#####
##################################################
#Get all species associated with pathway
reactome_conservation = list()
for(i in 1:length(reactome)){
  print(paste(i, 'out of', length(reactome)))
  
  #Percent
  percent = lapply(pathway_conservation, function(x) 
    x$percent_conserved[grep(paste('^', names(reactome)[i], '$', sep = ''),
                             names(x$percent_conserved))])
  
  #Effect size
  eff = lapply(pathway_conservation, function(x) 
    x$eff_size[grep(paste('^', names(reactome)[i], '$', sep = ''), 
                             names(x$percent_conserved))])
  
  #P-value
  p = lapply(pathway_conservation, function(x) 
    x$p_value[grep(paste('^', names(reactome)[i], '$', sep = ''), 
                    names(x$percent_conserved))])
  
  #Phylo distance
  p_dist = cophen[match(names(percent), 
                        names(cophen))]
  
  reactome_conservation[[names(reactome)[i]]] = list(percent_conserved = unlist(percent),
                                                     eff_size = unlist(eff),
                                                     p_value = unlist(p),
                                                     phylogenetic_distance = p_dist)
}

#Fit linear models
mods = list()
plot(reactome_conservation[[1]]$phylogenetic_distance, 
     reactome_conservation[[1]]$percent_conserved,
     type = 'n',
     xlab = 'Cophenetic distance',
     ylab = 'Percent conserved',
     cex.lab = 1.5,
     cex.axis = 1.5,
     ylim = c(0, 1))
for(i in 1:length(reactome_conservation)){
  if(length(reactome_conservation[[i]]$percent_conserved)>0){
    if(sum(is.na(reactome_conservation[[i]]$percent_conserved))<length(reactome_conservation[[i]]$percent_conserved)){
      
      y = reactome_conservation[[i]]$percent_conserved
      x = reactome_conservation[[i]]$phylogenetic_distance
      
      x = x[!is.na(y)]
      y = y[!is.na(y)]
      
      y = y[!is.na(x)]
      x = x[!is.na(x)]
      

      mod = lm(y~poly(x, 2))
      
      newdat = data.frame(x = seq(0, 
                                  max(cophen), 
                                  length.out = 100))
      pred = predict(mod, 
                     newdata = newdat)
      pred[pred>1] = 1
      pred[pred<0] = 0

      lines(newdat$x, 
            pred, 
            col=alpha('gray50', 0.15)) 
      
      mods[[names(reactome_conservation)[i]]] = pred
    }
  }
}

#NMDS
d = dist(do.call(rbind, mods))
u = umap::umap(do.call(rbind, mods), 
               verbose = TRUE)

#Plot all
plot(reactome_conservation[[1]]$phylogenetic_distance, 
     reactome_conservation[[1]]$percent_conserved,
     type = 'n',
     xlab = 'Cophenetic distance',
     ylab = 'Percent conserved',
     cex.lab = 1.5,
     cex.axis = 1.5)
lines(mods[[i]]$model)
for(i in 1:length(mods)){
  x = predict()
  abline(mods[[i]])
}



#PGLS
i = 20
x = reactome_conservation[[i]]$percent_conserved
names(x) = names(reactome_conservation[[i]]$eff_size)
x[is.na(x)] = 0
tree_pgls = keep.tip(tree, tree$tip.label[tree$tip.label%in%names(x)])

#PGLS
x = x[match(tree_pgls$tip.label, names(x))]

#Create data frame
datos = data.frame(n_ogs = x,
                   Species = names(x))

#Create comparative data object
comp.data<-caper::comparative.data(tree_pgls, 
                                   datos, 
                                   names.col="Species", 
                                   vcv.dim=2, 
                                   warn.dropped=TRUE)

#PGLS
mod <- caper::pgls(n_ogs~1, 
                   data=comp.data,
                   lambda = 'ML')

#Plot phylogenetic residuals
plot(x, mod$residuals)


hist(sort(unlist(lapply(reactome_conservation, function(x) 
     max(x$phylogenetic_distance[x$percent_conserved>=0.95],
         na.rm = TRUE)))))

completeness = unlist(lapply(reactome_conservation, function(x) 
  max(x$phylogenetic_distance[x$percent_conserved>=0.95],
      na.rm = TRUE)))
completeness = data.frame(x = completeness)

x = pathway_conservation[[1]]$p_value
names(x) = names(pathway_conservation[[1]]$percent_conserved)

