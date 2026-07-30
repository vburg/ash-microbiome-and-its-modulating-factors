# =============================================================================
# Ash phyllosphere microbiome analyses
#
# Purpose:
#   Explore taxonomic composition, alpha diversity, beta diversity,
#   ordination, pairwise community comparisons, and selected taxon models
#   using phyloseq and microViz.
#
# Required objects before running this script:
#   ps4         - phyloseq object containing ASV counts, taxonomy, sample data,
#                 and a phylogenetic tree
#   blatt_umw   - sample-level leaf/environmental metadata
#   phylo_blatt - phyloseq object containing leaf/environmental metadata
#   var.ino     - data frame used to define groups for convex hulls
#   neuenv      - metadata object exported near the end of the script
#
# Expected metadata variables include:
#   Site, Sample, plot, and Vitalitaet_Sommer
#
# Repository-relative folders used below:
#   data/       - input metadata
#   data/processed/ - processed R objects
#   results/    - exported statistical results
#
# Create these folders before running the script if they do not already exist.
# Run the script from the repository root or from an RStudio project.
#
# Several analyses use permutations or stochastic procedures. For exact
# reproducibility, set a fixed random seed before those analyses, for example:
# set.seed(1234)
#
# Package versions should be recorded with renv::snapshot() or sessionInfo().
# =============================================================================

# Load packages used throughout the analysis ---------------------------------
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(vegan)
library(microbiome)
library(microViz)

# Remove samples that were excluded from the final analysis ------------------
ps4<-subset_samples(ps4,!(sample_names(ps4) %in% c("TH2-17","ST1-06","BB1-12")))
# The following lines show the corresponding exclusions from blatt_umw.
# They remain commented out because the metadata may already have been filtered.
#blatt_umw<-blatt_umw[!grepl("TH2-17",blatt_umw$Sample), ]
#blatt_umw<-blatt_umw[!grepl("ST1-06",blatt_umw$Sample), ]
#blatt_umw<-blatt_umw[!grepl("BB1-12",blatt_umw$Sample), ]

# Extract taxonomy and abundance tables for inspection -----------------------
tax<- as.data.frame(tax_table(ps4))
otu<-as.data.frame(otu_table(ps4))

# Summarize and agglomerate taxa at phylum and genus levels ------------------
table(tax_table(ps4)[, "Phylum"], exclude = NULL)
abund.phyla = sort(tapply(taxa_sums(ps4), tax_table(ps4)[, "Phylum"], sum), TRUE)
abund.phyla
(ps4_phylum <- tax_glom(ps4, taxrank="Phylum") )
(ps4_genus <- tax_glom(ps4, taxrank="Genus") )
plot_bar(ps4_phylum, x="Sample", fill="Phylum")

p = plot_bar(ps4_phylum, x="Site", fill="Phylum")
p+ geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

#### Visualization of relative bacterial phylum and genus abundances #########

# Plot the 40 most abundant phyla per sample, grouped by site.
ps4 %>%
  #ps_filter(site == 1) %>%
comp_barplot(
  tax_level = "Phylum", n_taxa = 40,
  bar_outline_colour = NA,
  sample_order = "bray",
  facet_by="Site",
  bar_width = 0.5,
  taxon_renamer = toupper
) + coord_flip()

####

# Plot the 40 most abundant genera per sample, grouped by site.
ps4 %>%
  #ps_filter(site == 1) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 40,
    bar_outline_colour = NA,
    sample_order = "bray",
    facet_by="Site",
    merge_other = FALSE,
    bar_width = 0.5,
    taxon_renamer = toupper
  ) + coord_flip() +
  theme(legend.position = "bottom")

#######################

# Duplicate genus-level plot retained from the original script.
ps4 %>%
  #ps_filter(site == 1) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 40,
    bar_outline_colour = NA,
    sample_order = "bray",
    facet_by="Site",
    merge_other = FALSE,
    bar_width = 0.5,
    taxon_renamer = toupper
  ) + coord_flip() +
  theme(legend.position = "bottom")

# Merge samples by site and visualize site-level genus composition.
phyloseq::merge_samples(ps4, group = "Site") %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 40,
    sample_order = c("TH2","TH1","ST1","SN1","NI1","MV1","HE1","BY3","BY2","BY1","BW1","BB1"),
    bar_width = 0.8
  ) +
  coord_flip() + labs(x = NULL, y = NULL)


###########################################

####

# Visualize genus composition for one selected site.
ps4 %>%
  ps_filter(Site == "BW1") %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 40,
    bar_outline_colour = NA,
    sample_order = "bray",
    #facet_by="Site",
    bar_width = 0.5,
    taxon_renamer = toupper
  ) + coord_flip() +
  theme(legend.position = "bottom")

#######################

# Filter the relative abundance of a selected genus in each sample.
mass-rel<-ps4 %>%
      transform_sample_counts(function(x) x / sum(x)) %>%
      psmelt() %>%
      filter(Genus == "Massilia") %>%
      group_by(Sample) %>%
      summarise(
             Massilia_percent = sum(Abundance) * 100,
             .groups = "drop"
        )

#################

# Apply a log(1 + x) transformation for selected diversity analyses.
pslog <- transform_sample_counts(ps4, function(x) log(1 + x))

# Alpha diversity -----------------------------------------------------------

# Calculate and visualize selected richness and diversity indices.
richness_data <-plot_richness(ps4, x="Site", measures=c("Chao1","Shannon"), color="Site")

plot_richness(ps4, x="Site", measures=c("Observed","Chao1","Shannon","Simpson","Fisher"), color="Site") + 
  theme_bw() + 
  geom_boxplot(lwd=0.9, alpha=0.7, aes(fill="Site")) 

richness_data + geom_boxplot()

alpharich<-estimate_richness(ps4,split = TRUE, measures=c("Observed","Chao1","Shannon","Simpson","Fisher"))
#

# Merge samples by site and calculate site-level diversity indices.
ps4_site <- merge_samples(ps4, "Site")

alpharich_site <- estimate_richness(ps4_site, measures=c("Observed", "Chao1", "Shannon", "Simpson", "Fisher"))

print(alpharich_site)

# Calculate descriptive statistics for each site.
summary_alpha <- alpharich %>%
  group_by(Site) %>%
  summarise(across(c(Observed, Chao1, se.chao1, Shannon, Simpson, Fisher),
                   list(
                     median = ~ median(.),
                     mean = ~ mean(.),
                     min = ~ min(.),
                     max = ~ max(.)
                   ),
                   .names = "{col}_{fn}"))

# Print and export the summary statistics.
print(summary_alpha)
write.csv(summary_alpha, "results/alpharich_summary.csv")

# Optional export/import steps retained as comments.
#write.csv(alpharich, "results/alpharich.csv")
#alpharich<-read.csv("results/alpharich.csv")
#alpharich<-tibble::column_to_rownames(alpharich,"X")

alpharich<-as.data.frame(alpharich)
plot_richness(pslog, x="Site", measures=c("Shannon", "Simpson"), color="Site")

# Add row names and selected metadata columns to the alpha-diversity table.
# The row order of alpharich and blatt_umw must match before using cbind().
rownames(alpharich) <- factor(rownames(alpharich))
alpharich <- mutate(alpharich, rownames(alpharich))
alpharich<-cbind(alpharich, blatt_umw$plot,blatt_umw$Sample)

# Inspect distributions and test assumptions before inferential analyses.
hist(alpharich$Fisher)
shapiro.test(alpharich$Simpson)

# One-way ANOVA and Tukey post hoc test for Shannon diversity.
anosha<-aov(Shannon~Site,alpharich)
summary(anosha)
tukey_anosha<-TukeyHSD(anosha)
print(tukey_anosha)

library(multcompView)
multcompLetters4(anosha,tukey_anosha)

# Kruskal-Wallis and Dunn tests for Simpson diversity.
anosimp<-kruskal.test(Simpson~Site,alpharich)
dunsimp<-FSA::dunnTest(alpharich$Simpson, alpharich$Site, method = "bonferroni") # Alternative: dunn.test::dunn.test

# Kruskal-Wallis and Dunn tests for Chao1 richness.
hist(alpharich$Chao1)
anochao<-kruskal.test(Chao1~Site,alpharich)
dunchao1<-FSA::dunnTest(alpharich$Chao1, alpharich$Site, method = "bonferroni")

####

# Generate compact-letter displays for Dunn test results.
PT = dunsimp$res

PT

library(rcompanion)

cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)

###################################################


# Kruskal-Wallis and Dunn tests for observed richness.
shapiro.test(alpharich$Observed)
aovobs<-kruskal.test(Observed~Site,alpharich)
dunobs<-dunn.test::dunn.test(alpharich$Observed, alpharich$Site, method = "bonferroni")
PT = dunobs$res

PT

library(rcompanion)

cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)

# Additional alpha-diversity calculations.
divo::li(ps4)

head(estimate_richness(ps4, measures=c("Chao1","Shannon", "Simpson","Fisher")))
div1 <- estimate_richness(ps4, measures=c("Observed","Chao1","Shannon", "Simpson","Fisher"))
write.table(div1, "Div_ps4.txt", sep = ";")

# Plot alpha diversity by site.
ggplot(alpharich, aes(x = Site, y = Observed, fill = Site)) +
       geom_boxplot(alpha = 0.7) +
       geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) +
       labs(title = "Observed-Werte nach Site",
                       y = "Observed-Wert",
                       x = "Site") +
       theme_minimal() +
       theme(axis.text.x = element_text(angle = 45, hjust = 1))


#################################################################

# Overall ANOSIM and PERMANOVA -----------------------------------------------

# Extract the site grouping variable.
presabs = get_variable(pslog, "Site")

# Test overall site differences using weighted UniFrac distances.
anosim(phyloseq::distance(pslog,"wuniFrac"),presabs) # Preference for the UniFrac method
#anosim(phyloseq::distance(pslog,"bray"),presabs)

presabs = get_variable(pslog, "Site")
#adonis(phyloseq::distance(pslog,"bray") ~ presabs)
vegan::adonis2(phyloseq::distance(pslog,"wuniFrac") ~ presabs)

# This function call was retained from the original script.
# The original note states that it does not work with matrix, phyloseq, or dist objects.
adonisplus::adonispost(data = ps4, which = presabs, alpha = 0.05)

# Network exploration -------------------------------------------------------

plot_net(pslog,  maxdist = 0.4, point_label = "Sample", color = "Site")
ig<-make_network(pslog, max.dist = 0.6)
plot_network(ig,pslog,  line_weight = 0.4, color = "Site")

# Principal coordinates analysis -------------------------------------------

# The object name is retained from the original script.
# The active distance is Bray-Curtis, although alternative comments mention weighted UniFrac.
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "Bray") # Alternatives noted originally: wuniFrac, pslog
evals <- out.wuf.log$values$Eigenvalues

plot_ordination(pslog, out.wuf.log, color = "Site") +
  labs(col = "Site") +
  coord_fixed(sqrt(evals[2] / evals[1]))

plot_ordination(pslog, out.wuf.log, color = "Site") +
  labs(col = "Site") +
  coord_fixed(sqrt(evals[2] / evals[1])) + geom_text(mapping = aes(label = Sample), size = 2, vjust = 1.2)

# Extract ordination coordinates and attach site labels.
pCoA.scores <- as.data.frame(out.wuf.log$vectors)
data.scores = subset(pCoA.scores, select = c(Axis.1,Axis.2))
data.scores$name <- rownames(data.scores)
data.scores$sample <- sample_data(pslog)$Site
data.scores

# Optional import of group definitions used for convex hulls.
#var.ino <- read.csv("data/Varianten.csv")

# var.ino must contain a column named Bakt with the 12 site/group labels.
str(var.ino)
var.ino$Bakt[2]
i <- 1

# Calculate the convex hull for the first group.
hull.data <- data.scores[data.scores$sample == var.ino$Bakt[i], ][chull(data.scores[data.scores$sample == 
                                                                                       var.ino$Bakt[i], c("Axis.1", "Axis.2")]), ]  # Hull values for group 1

# Calculate and combine convex hulls for all 12 groups.
for(i in 1:12) {

  grp.i <- data.scores[data.scores$sample == var.ino$Bakt[i], ][chull(data.scores[data.scores$sample == 
                                                                                     var.ino$Bakt[i], c("Axis.1", "Axis.2")]), ]
  hull.data <- rbind(hull.data, grp.i)

}
hull.data

library(RColorBrewer)
palette <- brewer.pal(12, "Sel3")


ggplot() +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_polygon(data=hull.data,aes(x=Axis.1,y=Axis.2,fill=sample,group=sample),alpha=0.30) + # Add convex hulls
  geom_point(data=data.scores,aes(x=Axis.1,y=Axis.2,colour=sample),size=2) +
  #scale_fill_manual(values = palette) + # Change the fill colour palette
  #scale_color_manual(values = palette) +
  theme_bw() 


#####

# Pairwise ANOSIM between selected site pairs --------------------------------
# Multiple alternatives are kept as comments; the last active assignment is used.
# Multiple-testing adjustment is not required for a single pairwise comparison.

filterVariant = c("BB1", "BW1")
filterVariant = c("BB1", "BY1")
#filterVariant = c("BB1", "BY2")
#filterVariant = c("BB1", "BY3")
#filterVariant = c("BB1", "HE1")
#filterVariant = c("BB1", "MV1")
#filterVariant = c("BB1", "NI1")
#filterVariant = c("BB1", "SN1")
#filterVariant = c("BB1", "ST1")
#filterVariant = c("BB1", "TH1")
#filterVariant = c("BB1", "TH2")

ps5 <- subset_samples(pslog, Site %in% filterVariant)
sample_data(ps5)

presabs = get_variable(ps5, "Site")
anosim(phyloseq::distance(ps5,"wuniFrac"),presabs)
#adonis2(phyloseq::distance(ps5,"wuniFrac") ~ presabs)


########


# Pairwise ANOSIM loop ------------------------------------------------------

anosim_list <- list()
sites <- c("BB1","BW1","BY1","BY2","BY3","HE1","MV1","NI1","SN1","ST1","TH1","TH2")

# Loop over every unique pair of sites.
for (i in 1:(length(sites) - 1)) {
  site1 <- sites[i]
  print(site1)
  for (j in (i + 1):length(sites)) {
    site2 <- sites[j]
    print(site2)

    # Subset the phyloseq object for the current site pair.
    ps5 <- subset_samples(pslog, Site %in% c(site1, site2))

    # Extract the "Site" variable.
    # This line is preserved from the original script.
    presabs <- get_variable(pslog, "Site")

    # Calculate the distance matrix.
    dist_matrix <- phyloseq::distance(ps5, "wuniFrac")

    # Run ANOSIM.
    anosim_result <- anosim(dist_matrix, presabs)

    # Store the result in a named list.
    result_name <- paste(site1, "vs", site2, sep = "_")
    anosim_list[[result_name]] <- anosim_result
  }
}

anres<-print(anosim_list)
plot(anosim_list$BB1_vs_TH1)
plot(anosim_list)

# Export pairwise ANOSIM and metadata results.
writeLines(anres, "results/pairanosim_loop.txt")

write.csv(neuenv, "results/Meta_dat.csv")


