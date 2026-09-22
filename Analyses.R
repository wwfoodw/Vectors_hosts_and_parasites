
library(readxl)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(phytools)
library(ape)
library(dplyr)
library(microbiomeutilities)
library(RColorBrewer)

#.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#.
#.  Data ----
#. 
#.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

load("COI/dada2_coi.RData")

# rename rownames
rownames(cat_seqtabAll) <- sapply(strsplit(rownames(cat_seqtabAll), "_"), `[`, 1)

metadata <- read.csv("rangitahua_metadata.csv")
rownames(metadata) <- metadata$Sequence_ID

# taxonomy
s <- read.table("COI/taxonomy/s.txt", sep=",", header=FALSE)
colnames(s) <- c("ASV", "s")
g <- read.table("COI/taxonomy/g.txt", sep=",", header=FALSE)
colnames(g) <- c("ASV", "g")
f <- read.table("COI/taxonomy/f.txt", sep=",", header=FALSE)
colnames(f) <- c("ASV", "f")
o <- read.table("COI/taxonomy/o.txt", sep=",", header=FALSE)
colnames(o) <- c("ASV", "o")
c <- read.table("COI/taxonomy/c.txt", sep=",", header=FALSE)
colnames(c) <- c("ASV", "c")
p <- read.table("COI/taxonomy/p.txt", sep=",", header=FALSE)
colnames(p) <- c("ASV", "p")
k <- read.table("COI/taxonomy/k.txt", sep=",", header=FALSE)
colnames(k) <- c("ASV", "k")

#arrange data frames into list
df_list <- list(p,c,o,f,g,s) # note, if memory is exhausted this needs to be broken up

#merge
taxa <- k[k$ASV %in% colnames(cat_seqtabAll), ]

for ( .df in df_list ) {
  taxa <-merge(taxa, .df, by.x="ASV", by.y="ASV", all=T)
}

taxa <- taxa[!duplicated(taxa$ASV), ]
rownames(taxa) <- taxa$ASV
taxa <- taxa[,-1]

gc() # clear unused memory
rm(k,p,c,o,f,g,s)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

arthro_taxa <- taxa[(taxa$Phylum == "Arthropoda"), ]
arthro_seqtabAll <- cat_seqtabAll[ , (colnames(cat_seqtabAll) %in% rownames(arthro_taxa)) ]

#

metadata <- metadata[rownames(metadata) %in% rownames(arthro_seqtabAll), ]
arthro_seqtabAll <- arthro_seqtabAll[rownames(arthro_seqtabAll) %in% rownames(metadata), ]

arthro_seqtabAll <- arthro_seqtabAll[ , colSums(arthro_seqtabAll) > 10 ]
arthro_seqtabAll <- arthro_seqtabAll[rowSums(arthro_seqtabAll) > 10 , ]

spider_metadata <- metadata[metadata$Classification == "spider", ]
spider_seqtabAll <- arthro_seqtabAll[rownames(arthro_seqtabAll) %in% rownames(spider_metadata), ]
spider_seqtabAll <- spider_seqtabAll[ , colSums(spider_seqtabAll) > 0]
spider_taxa <- taxa[rownames(taxa) %in% colnames(spider_seqtabAll), ]

##################
#
# remove obvious artefacts ----
#
###################

# Define a function to get the highest and second-highest values
get_top_two <- function(x) {
  sorted_values <- sort(x, decreasing = TRUE)
  return(sorted_values[1:2])
}
# Apply the function to each row of the dataframe
top2 <- apply(arthro_seqtabAll[, -1], 1, get_top_two)
# Create a new dataframe with the results
top2 <- data.frame(t(top2))
colnames(top2) <- c("Highest", "Second_Highest")
#

write.csv(arthro_seqtabAll, "arthro_seqtabAll.csv")
arthro_seqtabAll[arthro_seqtabAll/(apply(arthro_seqtabAll, 1, max, na.rm=TRUE)) != 1] <- 0 # remove all but highest value in row apply(seqtabAll, 1, max, na.rm=TRUE)
arthro_seqtabAll <- arthro_seqtabAll[ , colSums(arthro_seqtabAll) > 0 ]
write.csv(arthro_seqtabAll, "top_arthro_seqtabAll.csv")

#.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#.
#.  Make phyloseq object ----
#. 
#.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Create phyloseq object:

taxa$Domain <- NULL # reduce dimensions

ps <- phyloseq(otu_table(arthro_seqtabAll, taxa_are_rows=FALSE),
               tax_table(as.matrix(arthro_taxa)),
               sample_data(metadata),
               refseq(DNAStringSet(coi_cat_sequences))
)

ps@otu_table[ps@otu_table <= 10] <- 0
ps@otu_table <- ps@otu_table[rowSums(ps@otu_table) >0, ]

# Get row sums samples (might need to use this later to offset proportional data)
ps@sam_data$read_number <- rowSums(ps@otu_table)

# Check the updated sample data
head(ps@sam_data)

# Remove empty taxa
rownames(ps@tax_table)[rownames(ps@tax_table) != colnames(ps@otu_table)] # check!!!! --> SHOULD BE EMPTY
names(ps@refseq)[names(ps@refseq) != colnames(ps@otu_table)] # check!!!! --> SHOULD BE EMPTY
ps@tax_table <- ps@tax_table[colSums(ps@otu_table) >0,  ]
ps@refseq <- ps@refseq[colSums(ps@otu_table) >0,  ]
ps@otu_table <- ps@otu_table[, colSums(ps@otu_table) >0 ]

sort(colSums(ps@otu_table))

# Resolve the taxonomy of the Culex sp. assignments using Smith and Fonseca 2004 method
BrowseSeqs(ps@refseq)
ps@tax_table[,"Species"][rownames(ps@tax_table) == "Hap_2"] <- "Culex pervigilans"
ps@tax_table[,"Species"][rownames(ps@tax_table) == "Hap_1"] <- "Culex quinquefasciatus"

plot_taxa_heatmap(ps,
                  subset.top = 20,
                  VariableA = "Location_Description",
                  heatcolors = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(100),
                  transformation = "compositional")

library(reshape2)
# annotate sam data with mosquito ID
hap <- melt(ps@otu_table)
hap <- hap[hap$value != 0, ]
colnames(hap) <- c("Sample", "ASV", "value")

sample_names(ps) == hap$Sample
sample_names(ps) %in% hap$Sample
hap$Sample %in% sample_names(ps)
rownames(hap) <- hap$Sample
hap$Species <- ifelse(hap$ASV == "Hap_1",
                               "C quinquefasciatus",
                               ifelse(hap$ASV == "Hap_2",
                                      "C pervigilans",
                                      NA))

# Reorder to match sample names in the phyloseq object
mos_hap_2 <- hap[sample_names(ps),]
sample_names(ps) == mos_hap_2$Sample
ps@sam_data$Hap <- mos_hap_2$ASV
ps@sam_data$Species <- mos_hap_2$Species
View(ps@sam_data)

coi_alignment <- AlignSeqs(ps@refseq, anchor=NA, verbose=FALSE)
BrowseSeqs(coi_alignment)
writeXStringSet(coi_alignment, filepath = "coi_alignment.fasta")
write.csv(ps@otu_table, "coi_alignment.csv")

mosquito_ps <- subset_taxa(ps, Family == "Culicidae")
mosquito_ps@otu_table <- mosquito_ps@otu_table[rowSums(mosquito_ps@otu_table) > 0]
mosquito_ps@sam_data <- mosquito_ps@sam_data[rownames(mosquito_ps@sam_data) %in% rownames(mosquito_ps@otu_table) , ]

save.image(file = "Mosquito.1.RData")
load(file = "Mosquito.1.RData")

#.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#.
#.  PLANT Data ----
#. 
#.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

load("trnL/Plant.1.RData")

metadata <- mosquito_ps@sam_data

# taxonomy
s <- read.table("trnL/taxonomy/s.txt", sep=",", header=FALSE)
colnames(s) <- c("ASV", "s")
g <- read.table("trnL/taxonomy/g.txt", sep=",", header=FALSE)
colnames(g) <- c("ASV", "g")
f <- read.table("trnL/taxonomy/f.txt", sep=",", header=FALSE)
colnames(f) <- c("ASV", "f")
o <- read.table("trnL/taxonomy/o.txt", sep=",", header=FALSE)
colnames(o) <- c("ASV", "o")
c <- read.table("trnL/taxonomy/c.txt", sep=",", header=FALSE)
colnames(c) <- c("ASV", "c")
p <- read.table("trnL/taxonomy/p.txt", sep=",", header=FALSE)
colnames(p) <- c("ASV", "p")
k <- read.table("trnL/taxonomy/k.txt", sep=",", header=FALSE)
colnames(k) <- c("ASV", "k")

#arrange data frames into list
df_list <- list(p,c,o,f,g,s) # note, if memory is exhausted this needs to be broken up

#merge
taxa <- k[k$ASV %in% colnames(seqtabAll), ]

for ( .df in df_list ) {
  taxa <-merge(taxa, .df, by.x="ASV", by.y="ASV", all=T)
}
rownames(taxa) <- taxa$ASV
taxa <- taxa[,-1]

gc() # clear unused memory
rm(k,p,c,o,f,g,s)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

plant_taxa <- taxa[(taxa$Kingdom == "Viridiplantae"), ]
plant_seqtabAll <- seqtabAll[ , (colnames(seqtabAll) %in% rownames(plant_taxa)) ]

write.csv(as.data.frame(sort(rowSums(plant_seqtabAll))), "plant_read_counts.csv")

#

metadata <- metadata[rownames(metadata) %in% rownames(plant_seqtabAll), ]
plant_seqtabAll <- plant_seqtabAll[rownames(plant_seqtabAll) %in% rownames(metadata), ]

plant_seqtabAll[plant_seqtabAll <= 10] <- 0
plant_seqtabAll <- plant_seqtabAll[ , colSums(plant_seqtabAll) > 10 ]
plant_seqtabAll <- plant_seqtabAll[rowSums(plant_seqtabAll) > 10 , ]

#.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#.
#.  Make phyloseq object ----
#. 
#.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Create phyloseq object:

taxa$Domain <- NULL # reduce dimensions

ps <- phyloseq(otu_table(plant_seqtabAll, taxa_are_rows=FALSE),
               tax_table(as.matrix(plant_taxa)),
               sample_data(metadata),
               refseq(DNAStringSet(trnl_sequences))
)

# Get row sums samples (might need to use this later to offset proportional data)
ps@sam_data$read_number <- rowSums(ps@otu_table)

# Check the updated sample data
head(ps@sam_data)
View(ps@sam_data)

# Remove empty taxa
rownames(ps@tax_table)[rownames(ps@tax_table) != colnames(ps@otu_table)] # check!!!! --> SHOULD BE EMPTY
names(ps@refseq)[names(ps@refseq) != colnames(ps@otu_table)] # check!!!! --> SHOULD BE EMPTY
ps@tax_table <- ps@tax_table[colSums(ps@otu_table) >0,  ]
ps@refseq <- ps@refseq[colSums(ps@otu_table) >0,  ]
ps@otu_table <- ps@otu_table[, colSums(ps@otu_table) >0 ]

sort(colSums(ps@otu_table))

#

plant_alignment <- AlignSeqs(ps@refseq, anchor=NA, verbose=FALSE)
BrowseSeqs(plant_alignment)
writeXStringSet(plant_alignment, filepath = "plant_alignment.fasta")
write.csv(ps@otu_table, "plant_alignment.csv")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# make phylogenetic tree ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

phangAlign <- phyDat(as(plant_alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
# best molecular evolutionary model
(mT <- modelTest(phangAlign, treeNJ)) # find model with lowest BIC for the alignment object
# choose best model from the table according to BIC
bestmodel <- mT$Model[which.min(mT$BIC)]
bestmodel
# search for a better tree using NNI rearrangements
fit = pml(treeNJ, data=phangAlign, model=bestmodel)
fitGTR <- update(fit, k=4, inv=0.2) # use update to change parameters
# search for a better tree using NNI rearrangements
fitter <- optim.pml(fitGTR, optNni=TRUE)
fitGTR <- optim.pml(fitter, model="TVMe", multicore=TRUE, optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
# Bootstrap
bs = bootstrap.pml(fitGTR, bs=100, optNni=TRUE, control = pml.control(trace = 0))
BA <- plotBS(midpoint(fitGTR$tree), bs, p = 50, type="n") # option type=="n" just assigns the bootstrap values and return the tree without plotting it.

#Merge phylogenetic tree to phyloseq object:
ps<-merge_phyloseq(ps,BA)
ps<-merge_phyloseq(ps,plant_alignment)

plant_ps <- tip_glom(ps, h = 0.1, hcfun = cluster::agnes)
plot(plant_ps@phy_tree)

# summary stats ----
# How many of each mosquito species contained plant material?
table(plant_ps@sam_data$Species)
sort(rowSums(otu_table(plant_ps)))

#

save.image(file="plant_1.RData")
load(file="plant_1.RData")

####################################################################
#
#
# Remove sequences ----
#
#
####################################################################

# Ngāti Kurī must be notified of the intention of users for use DNA sequences 
# Consequently, sequence data is remove from the files prior to publishing

ls()
rm(arthro_seqtabAll,arthro_taxa,BA,bestmodel,bs,cat_seqtabAll,coi_alignment,coi_cat_sequences,df_list,dm,fit,fitGTR,fitter,get_top_two,hap,metadata,mos_hap_2,mosquito_ps,mT,phangAlign,plant_alignment,plant_seqtabAll,plant_taxa,ps,seqtabAll,spider_metadata,spider_seqtabAll,spider_taxa,taxa,top2,treeNJ,trnl_sequences)

plant_ps@refseq <- NULL
save.image(file="rangitahua_export_to_github.RData")

####################################################################
#
#
# Ordination ----
#
#
####################################################################

eso.dpcoa <- DPCoA(plant_ps)
eso.dpcoa
plot_ordination(plant_ps, eso.dpcoa, "samples", color = "Location_Description", shape = "Species")
plot_ordination(plant_ps, eso.dpcoa, "species", color = "Order")

# replot as vectors:

eso.dpcoa <- ordinate(plant_ps, method = "DPCoA")

library(grid)

sp <- plot_ordination(
  plant_ps,
  eso.dpcoa,
  type = "species",
  color = "Order",
  justDF = TRUE
)

# Identify the first two ordination-coordinate columns
axes <- names(sp)[1:2]

ggplot(sp) +
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = .data[[axes[1]]],
      yend = .data[[axes[2]]],
      colour = Order
    ),
    arrow = arrow(
      type = "closed",
      length = unit(2, "mm")
    ),
    linewidth = 0.5
  ) +
  coord_equal() +
  labs(x = axes[1], y = axes[2]) +
  theme_bw()

plot_ordination(plant_ps, eso.dpcoa, "samples", color = "Location_Description", shape = "Species") +  theme_bw()

####################################################################
#
#
# Heatmap ----
#
#
####################################################################

x <- plant_ps
x <- tip_glom(plant_ps, h = 0.2, hcfun = cluster::agnes)
x@otu_table <- 100 * x@otu_table/rowSums(x@otu_table) # convert to percent

library(phytools)

loc <- x@sam_data$Location_Description
names(loc) <- rownames(x@sam_data)

mat <- t(as.data.frame(x@otu_table))
loc <- loc[colnames(mat)]

sample_order <- order(loc)
mat2 <- mat[, sample_order]

phylo.heatmap(
  x@phy_tree,
  mat2
)

# The above, but split by location

split_mats <- lapply(
  split(colnames(mat), loc),
  function(s) mat[, s, drop = FALSE]
)

blank <- matrix(
  NA,
  nrow = nrow(mat),
  ncol = 1,
  dimnames = list(rownames(mat), " ")
)

mat_sep <- do.call(
  cbind,
  c(rbind(split_mats, list(blank)))
)

phylo.heatmap(
  x@phy_tree,
  mat_sep
)

plot_tree(x, label.tips = "Order", color="Order", plot.margin=0.1, title = "By Height")
plot_tree(x, label.tips = "taxa_names", color="Location_Description", plot.margin=0.1, title = "By Height")
