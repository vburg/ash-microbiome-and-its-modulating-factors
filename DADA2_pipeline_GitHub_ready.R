# angepasst ausgehend von DADA2 Pipeline Tutorial (1.12)
# https://benjjneb.github.io/dada2/tutorial.html

# Nutzung der primer_clipped Daten, aber vorher L?ufe nach R1 und R2 sortiert/getrennt
# conda activate mg37
# ?ber Script: bash commands-split.sh


# srun -N 1 -n 1 -p highmem --time=2-00:00:00 --pty bash -i
# conda activate dada2-v1.14

#library(Biostrings)
#library(ShortRead)
#library(reshape2)
#library(gridExtra)


library(dada2)
library(ggplot2)
library(phyloseq)

path <- "seq/MiSeq_IBF/PrimerClipped_split_subset"
list.files(path)

r1.fwd <- list.files(path, pattern="_R1_fwd.fastq", full.names = TRUE)
r1.rvs <- list.files(path, pattern="_R1_rvs.fastq", full.names = TRUE)
r2.fwd <- list.files(path, pattern="_R2_fwd.fastq", full.names = TRUE)
r2.rvs <- list.files(path, pattern="_R2_rvs.fastq", full.names = TRUE)

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(r1.fwd), "_"), `[`, 1)
sample.names
sample.names1 <- sapply(strsplit(sample.names, "-"), `[`, 3)
sample.names2 <- sapply(strsplit(sample.names, "-"), `[`, 4)
sample.names <- paste(sep="-", sample.names1,sample.names2)
sample.names

pdf("QualityProfile_R1_fwd.pdf")
plotQualityProfile(r1.fwd[1:2])
dev.off()

pdf("QualityProfile_R1_rev.pdf")
plotQualityProfile(r1.rvs[1:2])
dev.off()

pdf("QualityProfile_R2_fwd.pdf")
plotQualityProfile(r2.fwd[1:2])
dev.off()

pdf("QualityProfile_R2_rev.pdf")
plotQualityProfile(r2.rvs[1:2])
dev.off()


# trim and filter
filt.r1.fwd <- file.path(path, "filtered", paste0(sample.names, "_R1_fwd_filt.fastq"))
filt.r1.rvs <- file.path(path, "filtered", paste0(sample.names, "_R1_rvs_filt.fastq"))
filt.r2.fwd <- file.path(path, "filtered", paste0(sample.names, "_R2_fwd_filt.fastq"))
filt.r2.rvs <- file.path(path, "filtered", paste0(sample.names, "_R2_rvs_filt.fastq"))

names(filt.r1.fwd) <- sample.names
names(filt.r1.rvs) <- sample.names
names(filt.r2.fwd) <- sample.names
names(filt.r2.rvs) <- sample.names

orient1.out <- filterAndTrim(r1.fwd, filt.r1.fwd, r2.rvs, filt.r2.rvs, truncLen=c(225,155),
                             maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                             compress=TRUE, multithread=TRUE)
head(orient1.out)


orient2.out <- filterAndTrim(r2.fwd, filt.r2.fwd, r1.rvs, filt.r1.rvs, truncLen=c(150,225),
                             maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                             compress=TRUE, multithread=TRUE)
head(orient2.out)

pdf("QualityProfile_filt.pdf")
plotQualityProfile(filt.r1.fwd[1:2])
dev.off()

pdf("QualityProfile_filt_Rs.pdf")
plotQualityProfile(filt.r2.rvs[1:2])
dev.off()

# Learn errors
err.orient1.fwd <- learnErrors(filt.r1.fwd, multithread=TRUE)
err.orient1.rvs <- learnErrors(filt.r2.rvs, multithread=TRUE)

err.orient2.fwd <- learnErrors(filt.r2.fwd, multithread=TRUE)
err.orient2.rvs <- learnErrors(filt.r1.rvs, multithread=TRUE)

pdf("errors_orient1_forward.pdf")
plotErrors(err.orient1.fwd, nominalQ=TRUE)
dev.off()

pdf("errors_orient1_reverse.pdf")
plotErrors(err.orient1.rvs, nominalQ=TRUE)
dev.off()

pdf("errors_orient2_forward.pdf")
plotErrors(err.orient2.fwd, nominalQ=TRUE)
dev.off()

pdf("errors_orient2_reverse.pdf")
plotErrors(err.orient2.rvs, nominalQ=TRUE)
dev.off()

dada.o1.fwd <- dada(filt.r1.fwd, err=err.orient1.fwd, multithread=TRUE)
dada.o1.rvs <- dada(filt.r2.rvs, err=err.orient1.rvs, multithread=TRUE)
mergers.o1 <- mergePairs(dada.o1.fwd, filt.r1.fwd, dada.o1.rvs, filt.r2.rvs, verbose=TRUE)
print("Mergers Orient1")
head(mergers.o1[[1]])

dada.o2.fwd <- dada(filt.r2.fwd, err=err.orient2.fwd, multithread=TRUE)
dada.o2.rvs <- dada(filt.r1.rvs, err=err.orient2.rvs, multithread=TRUE)
mergers.o2 <- mergePairs(dada.o2.fwd, filt.r2.fwd, dada.o2.rvs, filt.r1.rvs, verbose=TRUE)
print("Mergers Orient2")
head(mergers.o2[[1]])

saveRDS(dada.o1.fwd, file.path(path, "Esche_MiSeq_dada.o1.fwd.rds"))
saveRDS(dada.o1.rvs, file.path(path, "Esche_MiSeq_dada.o1.rvs.rds"))
saveRDS(dada.o2.fwd, file.path(path, "Esche_MiSeq_dada.o2.fwd.rds"))
saveRDS(dada.o2.rvs, file.path(path, "Esche_MiSeq_dada.o2.rvs.rds"))

dada.o1.fwd
dada.o1.rvs
dada.o2.fwd
dada.o2.rvs

# construct sequence tables
seqtab.o1 <- makeSequenceTable(mergers.o1)
dim(seqtab.o1)
print("Sequence Length distribution, Orient1")
table(nchar(getSequences(seqtab.o1)))

seqtab.o2 <- makeSequenceTable(mergers.o2)
dim(seqtab.o2)
print("Sequence Length distribution, Orient2")
table(nchar(getSequences(seqtab.o2)))

# Track reads through the pipeline mitz neuer Funktion
getN <- function(x) sum(getUniques(x))
getTrack <- function(out, dadaFs, dadaRs, mergers, sample.names){
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
  # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
  rownames(track) <- sample.names
  return(track)
}
track1 <- getTrack(orient1.out, dada.o1.fwd, dada.o1.rvs, mergers.o1, sample.names)
track2 <- getTrack(orient2.out, dada.o2.fwd, dada.o2.rvs, mergers.o2, sample.names)
track1
track2

# merge
seqtab<- mergeSequenceTables(seqtab.o1, seqtab.o2, repeats="sum")
print("Merged sequence table")
dim(seqtab)
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
print("merged sequenceTable, with chimera removal")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim, file.path(path, "Esche_MiSeq_IBF_no_chimeras.rds"))

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "reference_databases/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa.RDP <- assignTaxonomy(seqtab.nochim, "reference_databases/rdp_train_set_18.fa.gz", multithread=TRUE)

#  species level assignments based on exact matching between ASVs and sequenced reference strains.
taxa <- addSpecies(taxa, "reference_databases/silva_species_assignment_v138.1.fa.gz")
taxa.RDP <- addSpecies(taxa.RDP, "reference_databases/rdp_species_assignment_18.fa.gz")

# Let's inspect the taxonomic assignments:
taxa.print <- taxa.RDP # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# save data
saveRDS(taxa, file.path(path, "Esche_MiSeq_IBF_taxa.rds"))
saveRDS(taxa.RDP, file.path(path, "Esche_MiSeq_IBF_taxa.RDP.rds"))
write.csv(seqtab.nochim, "Esche_MiSeq_IBF_no_chimeras.csv")
write.csv(taxa, "Esche_MiSeq_IBF_taxa.csv")
write.csv(taxa.RDP, "Esche_MiSeq_IBF_taxa.RDP.csv")

# ----

## R-studio

library(vegan)
library(ggplot2)
library(phyloseq)
theme_set(theme_bw())

taxa <- readRDS("Esche_MiSeq_IBF_taxa.RDP.rds") 
seqtab.nochim <- readRDS("Esche_MiSeq_IBF_no_chimeras.rds")

rownames(seqtab.nochim)
# rownames(seqtab.nochim) <- substring(rownames(seqtab.nochim),1,3)
# rownames(seqtab.nochim)
samples.out <- rownames(seqtab.nochim)
site <- substring(samples.out,1,3)
replicate <- sapply(strsplit(samples.out, "-"), `[`, 2)
replicate
samdf <- data.frame(Site=site, Replicate=replicate, Sample=samples.out)
rownames(samdf) <- samples.out
samdf
#write.csv(samdf, "Esche_MiSeq_meta.csv")

#samdf1 <- read.csv("Esche_MiSeq_meta.csv")
# rownames(samdf1) <- samples.out
# samdf1
str(samdf)

# We now construct a phyloseq object directly from the dada2 outputs.

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps		
sample_data(ps)

# Removing mt and cp
subset_taxa(ps, Family=="Poaceae")
subset_taxa(ps, Order=="Chloroplast")
subset_taxa(ps, Phylum=="Cyanobacteria/Chloroplast")
ps1 <-  subset_taxa(ps, (Order!="Chloroplast") | is.na(Order))
ps1 <-  subset_taxa(ps1, (Phylum!="Cyanobacteria/Chloroplast") | is.na(Phylum))
ps1 <-  subset_taxa(ps1, (Family!="Poaceae") | is.na(Family))
ps1
subset_taxa(ps1, Kingdom=="Eukaryota") # zur Kontrolle
sort(sample_sums(ps1))

# short names for our ASVs and keep the full DNA sequences for other purposes
dna <- Biostrings::DNAStringSet(taxa_names(ps1))
names(dna) <- taxa_names(ps1)
ps2 <- merge_phyloseq(ps1, dna)
taxa_names(ps2) <- paste0("ASV", seq(ntaxa(ps2)))
ps2
head(sample_data(ps2))
(asv_tab <- data.frame(otu_table(ps2)[1:10, 1:10]))
head(tax_table(ps2))
head(refseq(ps2))

(ps3 <- prune_taxa(taxa_sums(ps2) > 1, ps2))
sort(sample_sums(ps3))
saveRDS(ps3, "phyloseq3.rds")

# Execute prevalence filter, using `prune_taxa()` function
# alle raus, die nur in einer Probe vorkommen

prevdf = apply(X = otu_table(ps3),
               MARGIN = ifelse(taxa_are_rows(ps3), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps3),
                    tax_table(ps3))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# nur ?bersicht phyla

keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= 2)]
ps01 = prune_taxa(keepTaxa, ps3)
ps01
# Anzahl ASVs weniger als mit combined Auswertung; MiSeq_Bact: 1284
sort(sample_sums(ps01))

prevdf = apply(X = otu_table(ps01),
               MARGIN = ifelse(taxa_are_rows(ps01), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps01),
                    tax_table(ps01))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
saveRDS(ps01, "phyloseq01.rds")


## am cluster:
path <- "seq/MiSeq_IBF/PrimerClipped_split_subset"
ps01 <- readRDS(file.path(path, "phyloseq01.rds"))

# install.packages("remotes")
# remotes::install_github("vmikk/metagMisc")
#
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("metagenomeSeq")
# BiocManager::install("DECIPHER")

library(phyloseq)
sort(sample_sums(ps01))

# ---
# phylogenetic tree ?ber R - Nachteil tree ist nicht gerootet, daher werden
# immer andere ASVs als root gew?hlt und die Werte und die Ordination schwankt

#install.packages("phangorn")
library (phangorn)
library (DECIPHER)
# Here we first construct a neighbor-joining tree, and then fit a GTR+G+I 
# maximum likelihood tree using the neighbor-joining tree as a starting point.
alignment <- AlignSeqs(refseq(ps01), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align, "F81")
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
fitGTR
tree.rooted <- midpoint(fitGTR$tree)
detach("package:phangorn", unload=TRUE)

ps4 <- merge_phyloseq(ps01,phy_tree(fitGTR$tree))
ps4_rooted <- merge_phyloseq(ps01,phy_tree(tree.rooted))
ps4
ps4_rooted
phy_tree(ps4_rooted)
saveRDS(ps4_rooted, file.path(path, "phyloseq4-rooted.rds"))
saveRDS(ps4, file.path(path, "phyloseq4.rds"))
saveRDS(ps4, file.path(path, "phyloseq5.rds"))
