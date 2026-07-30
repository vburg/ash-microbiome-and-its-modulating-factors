# =============================================================================
# Ash phyllosphere microbiome: environmental ordination and diversity analyses
#
# Purpose:
#   Prepare abundance and environmental data, remove excluded samples, perform
#   indicator-species analysis, PCoA, RDA, dbRDA, diversity comparisons, and
#   export tables for Cytoscape and downstream analyses.
#
# Required objects before running this script:
#   ps2, pslog, blatt_umw, impute_data, num_env, num_env_filt, plot_env,
#   rda_env, col.vector, pCoA.scores, otu3, y_utm32, and MAP
#
# Some of these objects are created in earlier analysis scripts. Their origin,
# preprocessing, and sample order should be documented in the repository README.
#
# Repository-relative input files:
#   data/processed/ps4_umwelt_blatt.rds
#   data/Blatt_umwelt_data.csv
#   data/ps4_ASV.csv
#   data/ps4_taxa.csv
#   data/Div_ps4.csv
#   data/alpharich.csv
#
# Repository-relative outputs:
#   data/processed/final_ps_IBF.rds
#   results/cytoscape/otu_table.csv
#   results/cytoscape/tax_table.csv
#   results/cytoscape/sample_table.csv
#   results/cytoscape/CoreMicro/abundance_genus.csv
#   results/cytoscape/CoreMicro/taxa_genus.csv
#
# Create data/processed/, results/cytoscape/, and
# results/cytoscape/CoreMicro/ before running the corresponding export commands.
#
# Reproducibility notes:
#   - Run this script from the repository root or an RStudio project.
#   - Record package versions with renv::snapshot() or sessionInfo().
#   - Set a fixed random seed before permutation or stochastic procedures when
#     exact repeatability is required, for example set.seed(1234).
#   - Verify that abundance tables and metadata have identical sample names and
#     ordering after every filtering step.
#
# Analysis commands, parameters, object names, and active code were preserved.
# Only comments and computer-specific file paths were changed.
# =============================================================================

# Load packages --------------------------------------------------------------
library(phyloseq)
library(vegan)
library(tidyverse)
library(microViz)
library(microbiome)
library(ggplot2)

library(ecodist)
library(ape)

# Load the processed phyloseq object created by the preceding workflow.
ps1<-readRDS("data/processed/ps4_umwelt_blatt.rds")

sum(sample_sums(ps1))
sort(sample_sums(ps1))
ps1<-subset_samples(ps1,!(sample_names(ps1) %in% c("TH2-17","ST1-06","BB1-12"))) # Remove the specified samples from the phyloseq object.
sample_data(ps2)

 # Import environmental metadata. The Sample column must match sample names.
 env_dat<- read.csv("data/Blatt_umwelt_data.csv")


# Import the ASV abundance and taxonomy tables.
# Confirm their orientation before calculating distances or constructing objects.
otu<-read.csv("data/ps4_ASV.csv", header = TRUE, sep=",")
tax<-read.csv("data/ps4_taxa.csv", header=TRUE, sep = ",")



otu<-tibble::column_to_rownames(otu,"X")
otu_new<-otu[order(row.names(otu)),]
otu<-otu_new
otu<-otu[!grepl("TH2-17",rownames(otu)), ]
otu<-otu[!grepl("ST1-06",rownames(otu)), ]
otu<-otu[!grepl("BB1-12",rownames(otu)), ]

impute_data<-impute_data[!grepl("TH2-17",impute_data$Sample), ]
impute_data<-impute_data[!grepl("ST1-06",impute_data$Sample), ]
impute_data<-impute_data[!grepl("BB1-12",impute_data$Sample), ]

num_env<-num_env[!grepl("TH2-17",num_env$Sample), ]
num_env<-num_env[!grepl("ST1-06",num_env$Sample), ]
num_env<-num_env[!grepl("BB1-12",num_env$Sample), ]

num_env_filt<-num_env_filt[!grepl("ST1-06",num_env_filt$Sample), ]
num_env_filt<-num_env_filt[!grepl("BB1-12",num_env_filt$Sample), ]

otu_filt <- otu[!grepl("TH2", rownames(otu)), ]
otu_filt <- otu_filt[!grepl("HE1", rownames(otu_filt)), ]


plot_env <- plot_env[!grepl("TH2", plot_env$Sample), ]
plot_env <- plot_env[!grepl("HE1", plot_env$Sample), ]

# Export abundance, taxonomy, and sample tables for Cytoscape.
write.csv(otu, "results/cytoscape/otu_table.csv")
write.csv(tax, "results/cytoscape/tax_table.csv")
write.csv(impute_data, "results/cytoscape/sample_table.csv")

library(rcompanion)
library(FSA)

# Import previously calculated alpha-diversity values.
div.dat<-read.csv("data/Div_ps4.csv")
head(div.dat)

kruskal.test(div.dat, Chao1~Sample)
dunnTest(div.dat, Chao1~plot)
t.test(div.dat, Chao1~plot)

keep_cols<- c("Sample","plot")

plot_var_cols<- c("Sample","plot","Probennr","Baumnr","bodentyp","x_utm32","y_utm32","ph","KAK",
                  "ca_vorrat","k_vorrat","mg_vorrat","NA_vorrat","corg_vorrat","n_vorrat","c_carboNAt_vorrat",
                  "hangrichtung","MAT","MAP","HHN","dBHD_ashes","Tjul","RHJul")

tree_var_cols<-c("Sample","plot","Vitalitaet_Sommer","TM_BLATT","BLATTFLAECHE","VERHAELTNIS_BL_BB","SPEZFLAECHE", "SPEZMASSE",
                 "OSMOLALITAET", "GChl","Carotin","GChlA","GChl_Caro","KH","STAERKE","GASCORBAT","GPH","VAN","PC")
                 
AS_var_cols<-c("plot","AS","Prolin_AS_PRZ","PROTEIN","Asp","Thr","Ser","Glu","Gln","Gly","Ala","Cys","Val", 
                 "Met","Ile","Leu","Tyr","Phe","Gaba","His","Trp","Orn","Lys","NH4","Arg","Pro")

dbrda_var_cols<-c("plot","MAP","KAK","y_utm32","Tjul", "dBHD_ashes" , "NA_vorrat", "ca_vorrat" , "GChl" , "RHJul" , "x_utm32" , "PC")
  
num_env<-impute_data %>% select(where (~ is.numeric(.x)))#,
                            #all_of(keep_cols))
num_env<-impute_data %>% select(where (~ is.numeric(.x)), all_of(keep_cols))

plot_env<-impute_data %>% select(all_of(plot_var_cols))
  
tree_env<-impute_data %>% select(all_of(tree_var_cols))
  
AS_env<-  impute_data %>% select(all_of(AS_var_cols))

dbrda_env<- num_env_filt %>% select(all_of(dbrda_var_cols))

kruskal_env<- impute_data %>% select(all_of(dbrda_var_cols))

str(num_env)
  
fac_env <- impute_data %>% select_if(is_character)
num_env<-cbind(num_env, env_dat$plot,env_dat$Sample)
colnames(num_env)[colnames(num_env)=="plot"] <- "env_dat$plot"
num_env<-tibble::column_to_rownames(num_env,"env_dat$Sample")


num_env<-num_env %>%
  rename("plot"="env_dat$plot")

num_env_filt<-num_env[!grepl("HE1", num_env$Sample), ]
num_env_filt<-num_env_filt[!grepl("TH2", num_env_filt$Sample), ]

plot_env[!grepl("HE1", plot_env$Sample), ]

test<-sample_data(ps1)

numenv$RHJul<-NULL # Delete columns.

rda()


######################################################################################
##### dbRDA as an alternative:

# Fit a distance-based redundancy analysis using plot-level predictors.
db_rda_plot<-capscale(otu_dist ~.,data=plot_env,distance="bray", na.action =na.omit)


#dbrda(formula, data, distance = "euclidean", sqrt.dist = FALSE,
 #     add = FALSE, dfun = vegdist, metaMDSdist = FALSE,
  #    na.action = na.fail, subset = NULL, ...)


# dbRDA with amino-acid variables
db_rda_full<-capscale(otu_dist ~ AS + Prolin_AS_PRZ + PROTEIN + Asp + Thr + Ser + Glu + Gln + Gly + Ala + Cys + Val + Met + Ile + Leu + Tyr + Phe +
                        Gaba + His + Trp + Orn + Lys + NH4 + Arg + Pro,data=AS_env, na.action =na.omit)
db_rda_null<-capscale(otu_dist ~ 1, data=AS_env,distance="bray", na.action =na.omit)
db_rda_res<-ordiR2step(db_rda_null,scope=formula(db_rda_full), direction = "forward",na.action =na.omit)
aov.res<-anova(db_rda_res, by = "term")
plot(db_rda_res)
db_rda_res$call

ggordiplots::gg_envfit(
  db_rda_res,
  AS_env,
  groups = AS_env$plot,
  scaling = 2,
  choices = c(1, 2),
  perm = 999,
  alpha = 0.05,
  angle = 20,
  len = 0.5,
  unit = "cm",
  arrow.col = "black",
  pt.size = 2,
  plot = TRUE
)
#
# dbRDA with tree variables
db_rda_full2<-capscale(otu_dist ~Vitalitaet_Sommer + TM_BLATT + BLATTFLAECHE + VERHAELTNIS_BL_BB + SPEZFLAECHE + SPEZMASSE + OSMOLALITAET + GChl +
                         Carotin + GChlA + GChl_Caro + KH + STAERKE + GASCORBAT + GPH + VAN + PC,data=tree_env, na.action =na.omit)
db_rda_null2<-capscale(otu_dist ~ 1, data=tree_env,distance="bray", na.action =na.omit)
db_rda_res2<-ordiR2step(db_rda_null2,scope=formula(db_rda_full2), direction = "forward",na.action =na.omit)
aov.res2<-anova(db_rda_res2, by = "term")
db_rda_res2$call
plot(db_rda_res2)
#
# dbRDA with plot-level variables
db_rda_full3<-capscale(otu_dist_filt ~ x_utm32 + y_utm32 + ph + KAK + ca_vorrat + k_vorrat + mg_vorrat + NA_vorrat +
                         corg_vorrat + n_vorrat + c_carboNAt_vorrat + MAT + MAP + HHN + dBHD_ashes + Tjul + RHJul,data=plot_env, na.action =na.exclude)
db_rda_null3<-capscale(otu_dist_filt ~ 1, data=plot_env,distance="bray", na.action =na.omit)
db_rda_res3<-ordiR2step(db_rda_null3,scope=formula(db_rda_full3), direction = "forward",na.action =na.exclude)
aov.res3<-anova(db_rda_res3, by = "term")
plot(db_rda_res3)
db_rda_res3$call

#capscale(formula = otu_dist_filt ~ MAP + KAK + y_utm32 + Tjul + 
 #          dBHD_ashes + NA_vorrat + ca_vorrat + RHJul, data = plot_env, 
  #       distance = "bray", na.action = na.omit)


ggordiplots::gg_envfit(
  db_rda_res3,
  plot_env,
  groups = plot_env$plot,
  scaling = 2,
  choices = c(1, 2),
  perm = 999,
  alpha = 0.05,
  angle = 20,
  len = 0.5,
  unit = "cm",
  arrow.col = "black",
  pt.size = 2,
  plot = TRUE
)
#########
# dbRDA with the complete table
# The filtered object is used because TH2 and HE1 are removed during site-level imputation.

db_rda_full4<-capscale(otu_dist_filt ~ AS + Prolin_AS_PRZ + PROTEIN + Asp + Thr + Ser + Glu + Gln + Gly + Ala + Cys + Val + Met + Ile + Leu + Tyr + Phe +
                        Gaba + His + Trp + Orn + Lys + NH4 + Arg + Pro + Vitalitaet_Sommer + TM_BLATT + BLATTFLAECHE + VERHAELTNIS_BL_BB + SPEZFLAECHE + SPEZMASSE + OSMOLALITAET + GChl +
                        Carotin + GChlA + GChl_Caro + KH + STAERKE + GASCORBAT + GPH + VAN + PC + x_utm32 + y_utm32 + ph + KAK + ca_vorrat + k_vorrat + mg_vorrat + NA_vorrat +
                        corg_vorrat + n_vorrat + c_carboNAt_vorrat + MAT + MAP + HHN + dBHD_ashes + Tjul + RHJul,
                      data=num_env_filt, na.action =na.omit)
db_rda_null4<-capscale(otu_dist_filt ~ 1, data=num_env_filt,distance="bray", na.action =na.omit)
db_rda_res4<-ordiR2step(db_rda_null4,scope=formula(db_rda_full4), direction = "forward",na.action =na.omit)
aov.res4<-anova(db_rda_res4, by = "term")
plot(db_rda_res4)
db_rda_res4$call

ggordiplots::gg_envfit(
  db_rda_res4,
  dbrda_env,
  groups = num_env_filt$plot,
  scaling = 2,
  choices = c(1, 2),
  perm = 999,
  alpha = 0.9,
  angle = 25,
  len = 0.5,
  unit = "cm",
  arrow.col = "black",
  #pch.fill = "blue",
  pt.size = 1.5,
  plot = TRUE
)

ggord::ggord(db_rda_res4, grp_in = NULL, axes = c("1", "2"), xlims = NULL,ylims = NULL)


#######

anova.cca(rdares,step=1000, by = "term")
ordiplot(rdares, scaling = 2, type = "text")

ps_genus <- tax_glom(ps1, taxrank="Genus") 
ps_phylum <- tax_glom(ps1, taxrank="Phylum") 

ord_explore(ps_genus)

tps1 %>% tax_fix()



#####


# Save the final processed phyloseq object using a repository-relative path.
saveRDS(ps2, "data/processed/final_ps_IBF.rds")

(ps2_genus <- tax_glom(ps2, taxrank="Genus") )
write.csv(otu_table(ps2_genus), "results/cytoscape/CoreMicro/abundance_genus.csv")
write.csv(tax_table(ps2_genus), "results/cytoscape/CoreMicro/taxa_genus.csv")



library(fastANCOM)
otu<-as.data.frame(otu)
impute_data<-as.data.frame(impute_data)
impute_data$plot <- as.factor(impute_data$plot)

head(otu3)
# Run fastANCOM using the abundance matrix and site factor.
# The order of plot_factor must exactly match the rows or columns used in otu3.
fastANCOM(Y=otu3,x=plot_factor)
otu3<-as.matrix(otu)
otu
otu2<-t(otu)
tax
impute_data
plot_factor<-c(rep(1,26),rep(2,28),rep(3,28),rep(4,25),rep(5,18),rep(6,28),rep(7,29),rep(8,27),rep(9,29),rep(10,27),rep(11,29),rep(12,27))
plot_factor<-as.factor(plot_factor)
str(impute_data$plot)



cor(y_utm32,MAP, data=dbrda_env)
cor(dbrda_env$y_utm32, dbrda_env$MAP)
cor.test(dbrda_env$y_utm32, dbrda_env$MAP)