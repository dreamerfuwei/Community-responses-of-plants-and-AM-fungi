# Community-responses-of-plants-and-AM-fungi
#####################################################################
### Data processing code
## DADA2 ITS workflow
rm(list=ls())
library(scales)
packageVersion("scales")
library(Rcpp)
packageVersion("Rcpp")
library("stringi")
packageVersion("stringi")
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

## Change me to the directory containing the fataq files
path <- "D:/EDH2017_soil"
list.files(path)

## generate matched lists of the forward and reverse read files, as well as parsing out the sample name
fnFs <- sort(list.files(path, pattern = ".R1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = ".R2.fq", full.names = TRUE))

## cHANGE ME to your forward primer sequence
FWD <- "GTGARTCATCGAATCTTTG"
REV <- "TCCTCCGCTTATTGATATGC"

## To ensure we have the right primers, and the correct orienttion of the primers on the reads, we will verify the presence and orientation of these primers in the data
allOrients <- function(primer){
  # Create all orientation of the input sequence
  require(Biostrings)
  dna <- DNAString(primer) # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, complement = complement(dna), Reverse = reverse(dna), RevComp = reverseComplement(dna))
  return(sapply(orients, toString)) # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

## The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult, so we'll filter them in this step (only do this in this step)
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, compress=FALSE, multithread = TRUE)

## We are now ready to count the number of times the primer appear in the forward and reverse read, while considering all possible primer orientations.
primerHits <- function(primer, fn){
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

## After install cutadapt, we need to tell R the path to the cutadapt command
cutadapt <- "C:/Users/Wei Fu/Miniconda2/Scripts/cutadapt"
system2(cutadapt, args = "--version")

## We now create output filenames for the cutadapt-ed files, and define the parameters we are going to give the cutadapt commond.
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

## As a sanity check, we will count the presence of primers in the first cutadapt-ed sample:
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

## Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = ".R1.fq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = ".R2.fq", full.names = TRUE))

## Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "[.]")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

## Inspect read quality profiles
# We sart by visualiziling the quality profiles of the forward reads:
plotQualityProfile(cutFs[1:2])
# Now we visualize the quality profile of the reverse reads:
plotQualityProfile(cutRs[1:2])

## Filter and trim
# Assigning the filenames for the output of the filtered reads to be stored as fastq.gz files.
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = FALSE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

## Learn the error rates
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

# Visualize the estimated error rates as a sanity check.
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

## Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

## Sample inference
dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)

## Merge paired reads 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

## Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

## Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread=FALSE, verbose = TRUE)
table(nchar(getSequences(seqtab.nochim)))

## Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denosiedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

## Assign taxonomy with UNITE database
unite.ref <- "D:/EDH2017/sh_general_release_dynamic_s_all_02.02.2019.fasta"  # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

# Inspecting the taxonomic assignments:
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Write out
setwd("D:/EDH2017_soil")
write.csv(seqtab.nochim, file = "EDH2017.seqtab.nochim.csv")
write.csv(taxa, file = "EDH2017.taxonomy.csv")

###########################################################
## Compute alpha diversity indices
library(vegan)                                       
setwd("D:/EDH2017_soil/results")                                             # setting your working directory
asv.table <- read.csv(file = "asv_fungi_table_52000.csv", header = TRUE, row.names = 1)
head(asv.table)
t.asv.table <- t(asv.table)

richness <- specnumber(t.asv.table)                  # species richness
shannon_diversity <- diversity(t.asv.table)           # shannon entropy (base e)
simpson_diversity <- diversity(t.asv.table, "simpson")    # simpson diversity
inverse_simpson <- diversity(t.asv.table, "invsimpson")
pielou_evenness <- shannon_diversity/log(richness)      # pielou evenness
shannon_evenness <- shannon_diversity/richness        # shannon evenness
simpson_evenness <- simpson_diversity/richness        # simpson evenness

###########################################################
## Model fit using random forest
library(randomForest)
library(randomForestExplainer)
setwd("D:/Writing/Extreme drought/HST.AMF.2017/Data analyses/results/AMF/Model.fit")
data <- read.csv(file = "env.csv", row.names = 1)
head(data)

forest <- randomForest(AMF ~ Plant.richness + Moisture + ANPP + Bbiomass + C + N + AP + CN + pH, 
                       data = data, localImp = TRUE, importance = TRUE)
forest

min_depth_frame <- min_depth_distribution(forest)
head(min_depth_frame)

plot_min_depth_distribution(min_depth_frame)

importance_frame <- measure_importance(forest)
importance_frame

plot_multi_way_importance(importance_frame, size_measure = "no_of_nodes")

plot_multi_way_importance(importance_frame, x_measure = "mse_increase", 
                          y_measure = "node_purity_increase", 
                          size_measure = "p_value", no_of_labels = 10)

importance(forest, type = 1)

##################################################################
# Redundancy Analysis (RDA)
library(adespatial);packageVersion("adespatial")
library(vegan)
library(ggplot2)
setwd("D:/HST.AMF.2017/results/AMF/RDA")

# Read data
spe.amf <- read.csv(file = "AMF.Genus.csv", header = TRUE, row.names = 1)
spe.amf <- t(spe.amf)
head(spe.amf)

env <- read.csv(file = "env.csv", header = TRUE, row.names = 1)
head(env)
envplant <- read.csv(file = "envplant.csv", header = TRUE, row.names = 1)
head(envplant)
envchem <- read.csv(file = "envchem.csv", header = TRUE, row.names = 1)
head(envchem)
envtreament <- read.csv(file = "envtreatment.csv", header = TRUE, row.names = 1)
head(envtreament)

# db-RDA
(spe.amf.dbRDA <- capscale(spe.amf ~ Treatment + ANPP + C + CN + N + AP + pH + Richness, 
                           env, distance = "bray", comm = spe.amf))

anova(spe.amf.dbRDA, permutations = how(nperm = 9999))
anova(spe.amf.dbRDA, permutations = how(nperm = 9999), by = "axis")
summary(spe.amf.dbRDA)

(R2a.amf.dbRDA <- RsquareAdj(spe.amf.dbRDA)$adj.r.squared)
(R2 <- RsquareAdj(spe.amf.dbRDA)$r.squared)
vif.cca(spe.amf.dbRDA)

# using ggplot2 make figures
spe.amf.dbRDA.scaling1 <- summary(spe.amf.dbRDA, scaling = 1) # extract scaling 1
spe.amf.dbRDA.site <- data.frame(spe.amf.dbRDA.scaling1$sites)[, 1:2]
# spe.amf.dbRDA.site.constraints <- data.frame(spe.amf.dbRDA.scaling1$constraints)[, 1:2]
spe.amf.dbRDA.env <- data.frame(spe.amf.dbRDA.scaling1$biplot)[, 1:2]

# read group data
group <- read.delim('group.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
group$Treatment <- factor(group$Treatment, 
                          levels =c("Control", "CHR", "CHRR", "INT", "INTR"))
# merge group information
spe.amf.dbRDA.site$sample <- row.names(spe.amf.dbRDA.site)
spe.amf.dbRDA.site <- merge(spe.amf.dbRDA.site, group, by = 'sample')
head(spe.amf.dbRDA.site)

spe.amf.dbRDA.env$sample <- NA
spe.amf.dbRDA.env$group <- row.names(spe.amf.dbRDA.env)
spe.amf.dbRDA.env <- spe.amf.dbRDA.env[5:12, ]
head(spe.amf.dbRDA.env)

spe.amf.dbRDA.env.tre <- data.frame(spe.amf.dbRDA.scaling1$centroids)[, 1:2]
colnames(spe.amf.dbRDA.env.tre) <- c("T1", "T2")
head(spe.amf.dbRDA.env.tre)

p <- ggplot(spe.amf.dbRDA.site, aes(CAP1, CAP2)) +
            geom_point(aes(color = Treatment), size = 5) +
            scale_y_continuous(limits = c(-0.5, 0.5)) +
            scale_color_manual(values = c("#8DC63F", "#FBB040", "#F06B22", "#39B54A", "#009444")) +
            geom_vline(xintercept = 0, color = "gray", size = 0.5) +
            geom_hline(yintercept = 0, color = "gray", size = 0.5) +
            geom_point(data = spe.amf.dbRDA.env.tre, aes(x = T1, y = T2), color ="blue", size = 5) +
            geom_text(data = spe.amf.dbRDA.env.tre, aes(x = T1, y = T2, label = rownames(spe.amf.dbRDA.env.tre)), color ="blue", size = 3) +
            geom_segment(data = spe.amf.dbRDA.env, 
                         aes(x = 0,y = 0, xend = CAP1, yend = CAP2),
                         arrow = arrow(length = unit(0.1, "cm")), 
                         size = 0.3, color = "blue") +
            geom_text(data = spe.amf.dbRDA.env, aes(CAP1 * 1.1, CAP2 * 1.1, label = group), color ="blue", size = 3) +
            labs(x = 'db-RDA1 (18.81%)', y = 'db-RDA2 (3.50%)') + 
            theme(panel.grid = element_blank(), 
                  panel.background = element_rect(color = "black", fill ="transparent"), 
                  legend.title = element_text(), 
                  legend.key = element_rect(fill = 'transparent'), 
                  legend.position = c(0.85, 0.85))
p

p1 <- p + stat_ellipse(aes(color = Treatment), level = 0.65, linetype = 2)
p1

##########################################################################
# Network construction
library(igraph)
library(psych)
library(vegan)
library(Matrix)
library(RColorBrewer)

setwd("D:/Study/PhD/Writing/Extreme drought/HST.AMF.2017/Data analyses/results/AMF/networks/species.net")

otu.table <- read.csv(file = "AMF_plant_ASVtab.csv", header = T, row.names = 1)
head(t(otu.table))

# Filtering low abundacne OTUs, maintaining OTUs appear more than 6 times of all the samples
otu.table.filter <- otu.table[, specnumber(t(otu.table)) >= 3]
print(c(ncol(otu.table), "versus", ncol(otu.table.filter))) # OTUs left

# caculate pairwise 
occor <- corr.test(otu.table.filter, use = "pairwise", method = "spearman", adjust = "BH", alpha = .05)
occor.r <- occor$r
occor.p <- occor$p

#  Filter the association based on p-values and level of correlations
occor.r[occor.p > 0.01 | abs(occor.r) < 0.3] = 0

# Create igraph object
net.graph <- graph_from_adjacency_matrix(occor.r, mode = "undirected", weighted = TRUE, diag = FALSE)
net.graph

# Creating a vector to remove the isolated nodes (nodes with no interactions)
bad.vs <- V(net.graph)[degree(net.graph) == 0] 

# Removing the isolated nodes from the graph object using the function delete.vertices()
net.graph <- delete.vertices(net.graph, bad.vs)
net.graph

# Assign the igraph weight attribute to net.graph.weight
net.graph.weight <- E(net.graph)$weight

# remove the weight before plotting (may influence the layouts)
E(net.graph)$weight <-NA

edge.list <- as_edgelist(net.graph)
write.csv(edge.list, file = "species.net.edgelist.filter3.csv")
# Plot the graph object
set.seed(123)
plot(net.graph, main = "AMF-PLANT NETWORK FILTER5", edge.lty = 1, margin=c(0,0,0,0),
     vertex.size = 5, vertex.frame.color = NA, edge.curved = T, edge.width = 1, 
     layout = layout_with_fr(net.graph), vertex.label = NA)

# compute node degrees (links) and use that to set node size:
deg <- degree(net.graph, mode = "all")
V(net.graph)$size <- deg*2
plot(net.graph, layout = layout_with_fr(net.graph), main = "AMF-PLANT NETWORK FILTER5", vertex.label = NA)

mycolors <- brewer.pal(12, "Paired") #library(RColorBrewer)
mycolors
# "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"

node.attr <- read.csv(file = "groups.csv", header = T, row.names = 1)
# set vertices color
net.graph.col <- node.attr[V(net.graph)$name, ]
net.graph.col$Genus <- factor(net.graph.col$Genus, levels = c("Ambispora", "Archaeospora", "Claroideoglomus", "Diversispora", "Dominikia", "Glomus", 
                                                              "Kamienskia", "Paraglomus", "Funneliformis", "Rhizophagus", "Septoglomus", "Unclassified", "plant"))
levels(net.graph.col$Genus) <- c("#A6CEE3", "#FFFF99", "gray80", "gray30", "#1F78B4", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FB9A99", "#B15928", "#33A02C")
V(net.graph)$color <- as.character(net.graph.col$Genus)

plot(net.graph, edge.curved = T, vertex.frame.color = NA, layout = layout_with_fr(net.graph), 
     main = "Plant-AM fungal co-occurrence network", vertex.label = NA, edge.width	= net.graph.weight*3, 
     edge.color = "gray80", vertex.size = deg*1.5)

legend(x = -1.5, y = -0.5, c("Ambispora", "Archaeospora", "Claroideoglomus", "Diversispora", "Dominikia", "Glomus", "Kamienskia", "Paraglomus", "Funneliformis", "Rhizophagus", "Septoglomus", "Unclassified", "Plant"), pch=21, col="#777777", pt.bg = levels(net.graph.col$Genus), pt.cex=2, cex=.8, bty="n", ncol=1)

# network attributes
# 1. number of edges
num.edges <- length(E(net.graph))
num.edges

positive.edges <- sum(net.graph.weight > 0) # positive edges
positive.edges

negative.edges <- sum(net.graph.weight < 0) # negative edges
negative.edges

# 2. number of vertices
num.vertices <- length(V(net.graph))
num.vertices

# 3. number of connectance
connectance <- edge_density(net.graph, loops = FALSE) # loops = TRUE means self-connectance i.e A-A, B-B
connectance

# 4. average degree
average.degree <- mean(igraph::degree(net.graph))
average.degree

# 5. average path length
average.path.length <- average.path.length(net.graph) # or mean_distance(net.graph)
average.path.length

# 6.diameter
diameter <- diameter(net.graph, directed = FALSE, unconnected = TRUE, weights = NULL)
diameter

# 7. adge connectivity / group adhesion
edge.connectivity <- edge_connectivity(net.graph)
edge.connectivity

# 8. clustering coefficient
clustering.coefficient <- transitivity(net.graph)
clustering.coefficient

no.clusters <- no.clusters(net.graph)
no.clusters

# 9. betweenness centralization
centralization.betweenness <- centralization.betweenness(net.graph)$centralization
centralization.betweenness

# 10. degree centralization
centralization.degree <- centralization.degree(net.graph)$centralization
centralization.degree
