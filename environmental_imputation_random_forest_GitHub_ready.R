# =============================================================================
# Environmental-data imputation and random-forest analyses
#
# Purpose:
#   Import and clean leaf/environmental measurements, remove variables with
#   excessive missingness, impute missing numeric values, inspect collinearity,
#   prepare environmental and microbiome predictors, and fit random-forest
#   models for ash vitality.
#
# Required input files:
#   data/Datenbank_ZALF.xlsx
#   data/disease_state.tsv
#   data/processed/phyloseq4.rds
#
# Required objects before running all sections:
#   env_dat, num_env, ps, table, df, and any metadata columns referenced below.
#   Several of these originate from earlier scripts and should be documented
#   in the repository README.
#
# Repository-relative output:
#   results/Umwelt_data_interpoliert.csv
#
# Reproducibility notes:
#   - Run this script from the repository root or an RStudio project.
#   - Record package versions with renv::snapshot() or sessionInfo().
#   - The random seed is fixed at 42 before data splitting and model fitting.
#   - Document the number of samples and predictors retained after filtering,
#     imputation, merging, and train/test partitioning.
#   - Verify that microbiome tables, environmental data, and response variables
#     use identical sample identifiers and ordering before merging.
#   - Record the exact versions of caret, ranger, randomForest, fastDummies,
#     microbiome, and phyloseq used for the analysis.
#
# Important:
#   Executable commands, formulas, thresholds, object names, and model settings
#   were preserved. Only comments and computer-specific paths were changed.
# =============================================================================

# Load packages --------------------------------------------------------------
library(readxl)
library(dplyr)
library(randomForest) 
library(caret)
library(fastDummies)
library(phyloseq)
library(ggplot2)
library(vegan)
library(reshape2)

# Read the pre-cleaned Excel file.
# Import the cleaned database and retain only samples relevant to the microbiome study.
data <- read_excel("data/Datenbank_ZALF.xlsx") %>% 
  filter(MATERIAL == "Blatt") %>% # Keep only leaf material.
  filter(Site != "SN2") %>% # Remove site SN2 because it is not represented in the microbiome samples.
  tibble::column_to_rownames("SampleName")
head(env_dat)

# Column names of the initial data frame:
# [1] "Probennr"                             "Baumnr"                               "Site"                                 "DATUM"                                "MATERIAL"                             "BHD [cm]"                             "Tendenz_Vitalität_Winter_davor"      
# [8] "Vitalität_Winter_davor"               "Wasserreiser_Krone_1"                 "Wasserreiser_Stamm_1"                 "Insektenbefall_1"                     "Stammfußnekrose_1"                    "Pilzfruchtkörper_1"                   "Rhizomorphen_1"                      
# [15] "Tendenz_Vitalität_Sommer"             "Vitalität_Sommer"                     "LFE_Tendenz_Vitalität_Sommer"         "LFE_Vitalität_Sommer"                 "Tendenz_Vitalität_Winter_danach"      "Vitalität_Winter_danach"              "Wasserreiser_Krone_2"                
# [22] "Wasserreiser_Stamm_2"                 "Insektenbefall_2"                     "Stammfußnekrose_2"                    "Pilzfruchtkörper_2"                   "Rhizomorphen_2"                       "FM_BLATT[mg]"                         "TM_BLATT[mg]"                        
# [29] "TM[%]"                                "BLATTFLÄCHE[mm²]"                     "BLATTLÄNGE[mm]"                       "BLATTBREITE[mm]"                      "VERHÄLTNIS_BL_BB"                     "SPEZ.FLÄCHE[cm²/g]"                   "SPEZ.MASSE[mg/cm²]"                  
# [36] "OSMOLALITÄT[osmol/kg]"                "GChl[mg/g TM]"                        "Chla[mg/g TM]"                        "Chlb[mg/g TM]"                        "Carotin[mg/g TM]"                     "GChlA[mg/mm²]"                        "ChlaA"                               
# [43] "ChlbA"                                "CarotinA"                             "GChl_Caro"                            "Chla/b"                               "GChl/BLATT[µg]"                       "KH[mg/g TM]"                          "KHE[µmol Glycosyleinh./g TM]"        
# [50] "KH/GChl"                              "KH/BLATT[µg/Blatt]"                   "AS[µmol/g TM]"                        "KHE/AS"                               "PROLIN[µmol/g TM]"                    "Prolin_AS[% AS]"                      "STÄRKE[mg/g TM]"                     
# [57] "STÄRKEE[µmol Glycosyleinh./g TM]"     "SUM_KHESTÄE[µmol Glycosyleinh./g TM]" "KH/STÄRKE"                            "KHE/STÄRKEE"                          "STÄRKE/BLATT[µg]"                     "PROTEIN[mg/g TM]"                     "PROTEIN/AS"                          
# [64] "GASCORBAT[mg/g TM]"                   "GASC/CHLges"                          "GPH[µmol/g TM]"                       "E280[E280/mg TM/ml]"                  "VAN[µmol/g TM]"                       "VAN_608"                              "PC[µmol/g TM]"                       
# [71] "ODHP[µmol/g TM]"                      "Asp[µmol/g TM]"                       "Thr[µmol/g TM]"                       "Ser[µmol/g TM]"                       "Glu[µmol/g TM]"                       "Gln[µmol/g TM]"                       "Gly[µmol/g TM]"                      
# [78] "Ala[µmol/g TM]"                       "Cys[µmol/g TM]"                       "Val[µmol/g TM]"                       "Met[µmol/g TM]"                       "Ile[µmol/g TM]"                       "Leu[µmol/g TM]"                       "Tyr[µmol/g TM]"                      
# [85] "Phe[µmol/g TM]"                       "Gaba[µmol/g TM]"                      "His[µmol/g TM]"                       "Trp[µmol/g TM]"                       "Orn[µmol/g TM]"                       "Lys[µmol/g TM]"                       "NH4[µmol/g TM]"                      
# [92] "Arg[µmol/g TM]"                       "Pro[µmol/g TM]"                       "SUMASA[µmol/g TM]"                    "Asp[%AS]"                             "Thr[%AS]"                             "Ser[%AS]"                             "Glu[%AS]"                            
# [99] "Gln[%AS]"                             "Gly[%AS]"                             "Ala[%AS]"                             "Cys[%AS]"                             "Val[%AS]"                             "Met[%AS]"                             "Ile[%AS]"                            
# [106] "Leu[%AS]"                             "Tyr[%AS]"                             "Phe[%AS]"                             "Gaba[%AS]"                            "His[%AS]"                             "Trp[%AS]"                             "Orn[%AS]"                            
# [113] "Lys[%AS]"                             "NH4[%AS]"                             "Arg[%AS]"                             "Pro[%AS]"                             "Temp_Probennahme°C"                   "Luftf_Probennahme%"                  


# Identify columns with more than 20% missing values.
# Calculate the proportion of missing values in each column.
na_percentage <- colMeans(is.na(env_dat), na.rm = TRUE)
columns_to_remove <- names(na_percentage[na_percentage > 0.2])
columns_to_remove

# Remove the 18 columns exceeding the missing-value threshold. 
# [1] "BHD [cm]"                       "Tendenz_Vitalität_Winter_davor" "Vitalität_Winter_davor"         "Wasserreiser_Krone_1"           "Wasserreiser_Stamm_1"           "Insektenbefall_1"               "Stammfußnekrose_1"              "Pilzfruchtkörper_1"             "Rhizomorphen_1"                
# [10] "Wasserreiser_Krone_2"           "Wasserreiser_Stamm_2"           "Insektenbefall_2"               "Stammfußnekrose_2"              "Pilzfruchtkörper_2"             "Rhizomorphen_2"                 "VAN_608"                        "Temp_Probennahme°C"             "Luftf_Probennahme%"    

data2 <- env_dat %>% dplyr::select(-all_of(columns_to_remove))
# Count missing values by site and column.
na_by_site <- data2 %>%
  group_by(plot) %>%
  summarize(across(everything(), ~ sum(is.na(.))))
# TH2 has no records for previous vitality, and NI1 has no records for subsequent vitality.
# Identify sites with the highest numbers of missing values.
tibble::column_to_rownames(na_by_site, "plot") %>% rowSums()

colnames(data2)
# [1] "Probennr"                             "Baumnr"                               "Site"                                 "DATUM"                                "MATERIAL"                             "Tendenz_Vitalität_Sommer"             "Vitalität_Sommer"                    
# [8] "LFE_Tendenz_Vitalität_Sommer"         "LFE_Vitalität_Sommer"                 "Tendenz_Vitalität_Winter_danach"      "Vitalität_Winter_danach"              "FM_BLATT[mg]"                         "TM_BLATT[mg]"                         "TM[%]"                               
# [15] "BLATTFLÄCHE[mm²]"                     "BLATTLÄNGE[mm]"                       "BLATTBREITE[mm]"                      "VERHÄLTNIS_BL_BB"                     "SPEZ.FLÄCHE[cm²/g]"                   "SPEZ.MASSE[mg/cm²]"                   "OSMOLALITÄT[osmol/kg]"               
# [22] "GChl[mg/g TM]"                        "Chla[mg/g TM]"                        "Chlb[mg/g TM]"                        "Carotin[mg/g TM]"                     "GChlA[mg/mm²]"                        "ChlaA"                                "ChlbA"                               
# [29] "CarotinA"                             "GChl_Caro"                            "Chla/b"                               "GChl/BLATT[µg]"                       "KH[mg/g TM]"                          "KHE[µmol Glycosyleinh./g TM]"         "KH/GChl"                             
# [36] "KH/BLATT[µg/Blatt]"                   "AS[µmol/g TM]"                        "KHE/AS"                               "PROLIN[µmol/g TM]"                    "Prolin_AS[% AS]"                      "STÄRKE[mg/g TM]"                      "STÄRKEE[µmol Glycosyleinh./g TM]"    
# [43] "SUM_KHESTÄE[µmol Glycosyleinh./g TM]" "KH/STÄRKE"                            "KHE/STÄRKEE"                          "STÄRKE/BLATT[µg]"                     "PROTEIN[mg/g TM]"                     "PROTEIN/AS"                           "GASCORBAT[mg/g TM]"                  
# [50] "GASC/CHLges"                          "GPH[µmol/g TM]"                       "E280[E280/mg TM/ml]"                  "VAN[µmol/g TM]"                       "PC[µmol/g TM]"                        "ODHP[µmol/g TM]"                      "Asp[µmol/g TM]"                      
# [57] "Thr[µmol/g TM]"                       "Ser[µmol/g TM]"                       "Glu[µmol/g TM]"                       "Gln[µmol/g TM]"                       "Gly[µmol/g TM]"                       "Ala[µmol/g TM]"                       "Cys[µmol/g TM]"                      
# [64] "Val[µmol/g TM]"                       "Met[µmol/g TM]"                       "Ile[µmol/g TM]"                       "Leu[µmol/g TM]"                       "Tyr[µmol/g TM]"                       "Phe[µmol/g TM]"                       "Gaba[µmol/g TM]"                     
# [71] "His[µmol/g TM]"                       "Trp[µmol/g TM]"                       "Orn[µmol/g TM]"                       "Lys[µmol/g TM]"                       "NH4[µmol/g TM]"                       "Arg[µmol/g TM]"                       "Pro[µmol/g TM]"                      
# [78] "SUMASA[µmol/g TM]"                    "Asp[%AS]"                             "Thr[%AS]"                             "Ser[%AS]"                             "Glu[%AS]"                             "Gln[%AS]"                             "Gly[%AS]"                            
# [85] "Ala[%AS]"                             "Cys[%AS]"                             "Val[%AS]"                             "Met[%AS]"                             "Ile[%AS]"                             "Leu[%AS]"                             "Tyr[%AS]"                            
# [92] "Phe[%AS]"                             "Gaba[%AS]"                            "His[%AS]"                             "Trp[%AS]"                             "Orn[%AS]"                             "Lys[%AS]"                             "NH4[%AS]"                            
# [99] "Arg[%AS]"                             "Pro[%AS]" 
# Orn represents ornithine, a non-proteinogenic amino acid.
# Amino-acid columns used in data2.

cols.aa.gTM <- c("PROTEIN","PROTEIN_AS","GPH","E280","VAN","PC","ODHP","Asp","Thr","Ser","Glu","Gln","Gly","Ala",              
                 "Cys","Val","Met","Ile","Leu","Tyr","Phe","Gaba","His","Trp","Orn","Lys","NH4","Arg","Pro","SUMASA")

cols.sitevar<- c("Probennr", "bodentyp","x_utm32","y_utm32","ph","KAK","ca_vorrat","k_vorrat","mg_vorrat","NA_vorrat","corg_vorrat",
                 "n_vorrat","c_carboNAt_vorrat","hangrichtung","MAT","MAP","HHN","dBHD_ashes","Tjul","RHJul","ecoG")


cols.vitality <- c(  "Vitalität_Sommer" 
                     # "Tendenz_Vitalität_Sommer","LFE_Tendenz_Vitalität_Sommer" ,"LFE_Vitalität_Sommer"  , "Tendenz_Vitalität_Winter_danach" ,
)
cols.leaf.measures <- c("TM_BLATT" , "BLATTFLAECHE","VERHAELTNIS_BL_BB" , 
                        "SPEZFLAECHE" 
                        ,"TM", "OSMOLALITAET" #"BLATTLAENGE","BLATTBREITE" ,"FM_BLATT",
)
cols.molecules <- c( "Carotin" ,"GChlA" ,"CarotinA" ,
                     "GChl_Caro" ,"KH", 
                     "AS", "PROLIN"  ,"STAERKE",
                     "KH" ,"PROTEIN" , 
                     "GASCORBAT","E280" ,"VAN", "PC" ,
                     "NH4",  "Gaba"
                     #"NH4[%AS]","Gaba[%AS]","PROTEIN/AS" , "GASC/CHLges" , "Prolin_AS[% AS]" ,"Chla/b","KH/GChl" ,"KH/BLATT[µg/Blatt]", "GChl_Caro" , 
                     #"STÄRKE/BLATT[µg]" ,"SUM_KHESTÄE[µmol Glycosyleinh./g TM]" ,"KHE/AS" ,"KHE[µmol Glycosyleinh./g TM]" ,"ODHP[µmol/g TM]" ,
                     #"STÄRKEE[µmol Glycosyleinh./g TM]"  ,"KHE/STÄRKEE"  ,"GChl[mg/g TM]" , "Chla[mg/g TM]",
)
# Impute missing values using the median of each column.
numerical_cols <- data2 %>%  select_if(is.numeric) %>% colnames()
# Create a globally median-imputed version of the numeric variables.
data2.imputed.all <- data2 %>% mutate_at(vars(one_of(numerical_cols)), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))
# Alternatively, impute missing values using the median within each site.
# Create a site-wise median-imputed data set and use Sample as row names.
data2.imputed.site <- data2 %>% group_by(plot) %>%  
  mutate_at(vars(one_of(numerical_cols)), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>% 
  ungroup() %>% tibble::column_to_rownames("Sample")

# Check collinearity within different predictor subsets.

data2.imputed.site %>%  #.site
  dplyr::select(cols.molecules) %>% 
  #select_if(is.numeric) %>% 
  cor(use="na.or.complete") %>% abs() %>%
  melt() %>%
  ggplot(aes(x=Var1, y = Var2, fill=value, label=round(value, 2))) + geom_tile() + geom_text(color="black") +
  scale_fill_gradient2(low="lightblue", mid="white", high="red", midpoint=0.5, limits=c(0, 1), name="Correlation") +
  theme(axis.text.x = element_text(angle = 90))

data2.imputed.site %>%  #.site
  dplyr::select(cols.leaf.measures) %>% 
  #select_if(is.numeric) %>% 
  cor(use="na.or.complete") %>% abs() %>%
  melt() %>%
  ggplot(aes(x=Var1, y = Var2, fill=value, label=round(value, 2))) + geom_tile() + geom_text(color="black") +
  scale_fill_gradient2(low="lightblue", mid="white", high="red", midpoint=0.5, limits=c(0, 1), name="Correlation") +
  theme(axis.text.x = element_text(angle = 90))

data2.imputed.site %>%  #.site
  dplyr::select(cols.sitevar) %>% 
  #select_if(is.numeric) %>% 
  cor(use="na.or.complete") %>% abs() %>%
  melt() %>%
  ggplot(aes(x=Var1, y = Var2, fill=value, label=round(value, 2))) + geom_tile() + geom_text(color="black") +
  scale_fill_gradient2(low="lightblue", mid="white", high="red", midpoint=0.5, limits=c(0, 1), name="Correlation") +
  theme(axis.text.x = element_text(angle = 90))

data2.imputed.site %>%  #.site
  dplyr::select(cols.aa.gTM) %>% 
  #select_if(is.numeric) %>% 
  cor(use="na.or.complete") %>% abs() %>%
  melt() %>%
  ggplot(aes(x=Var1, y = Var2, fill=value, label=round(value, 2))) + geom_tile() + geom_text(color="black") +
  scale_fill_gradient2(low="lightblue", mid="white", high="red", midpoint=0.5, limits=c(0, 1), name="Correlation") +
  theme(axis.text.x = element_text(angle = 90))


num_env %>%  #.site
  #dplyr::select(cols.aa.gTM) %>% 
  #select_if(is.numeric) %>% 
  cor(use="na.or.complete") %>% abs() %>%
  melt() %>%
  ggplot(aes(x=Var1, y = Var2, fill=value, label=round(value, 2))) + geom_tile() + geom_text(color="black") +
  scale_fill_gradient2(low="lightblue", mid="white", high="red", midpoint=0.5, limits=c(0, 1), name="Correlation") +
  theme(axis.text.x = element_text(angle = 90))


# Export the site-wise imputed environmental data.
write.csv(data2.imputed.site, "results/Umwelt_data_interpoliert.csv")
impute_data<-read.csv("results/Umwelt_data_interpoliert.csv")
numenv$RHJul<-NULL # Delete the RHJul column.

# Vitality correlations are moderate, approximately 0.5 between summer and subsequent winter vitality.
# For amino acids, use either percentage values or µmol/g dry mass values.
# TM_BLATT and FM_BLATT, specific mass and specific area, and leaf-area variables are correlated. 
# ODHP is strongly correlated with other variables, and STÄRKEE correlates with STÄRKE.
# Derive a vitality-trend variable. 
data2.imputed.site <- data2.imputed.site %>% 
  mutate(vitalityTrend = case_when(
    Vitalität_Sommer > Vitalität_Winter_danach ~ 1,
    Vitalität_Sommer < Vitalität_Winter_danach ~ -1,
    TRUE ~ 0))

data2.imputed.site <- data2.imputed.site %>% 
  mutate(diseaseState = ifelse(Vitalität_Sommer < 1.5,  
                               "one", 
                               ifelse(Vitalität_Sommer < 2.5, 
                                      "two",
                                      ifelse(Vitalität_Sommer < 3.5,
                                             "three", 
                                             "four-five"))))

# Standardize selected predictor groups and prepare the environmental predictor table.
rf.env.in <- data2.imputed.site %>% 
  mutate_at(vars(one_of(c(cols.leaf.measures, cols.molecules))), ~decostand(., "standardize")) %>%
  dplyr::select(c("plot",#"vitalityTrend", 
                  "Vitalitaet_Sommer", #"Vitalität_Winter_danach", 
                  cols.aa.gTM , cols.leaf.measures, cols.molecules,)) 

# Apply dummy encoding to the site variable.
rf.env.in <- rf.env.in %>% tibble::rownames_to_column("Sample") %>%
  dummy_cols(select_columns = "Site") %>% select(-c("Site")) %>% tibble::column_to_rownames("Sample")


## Microbiome data
disease.state <- read.table(file = "data/disease_state.tsv", sep="\t", header=TRUE)
# Add disease-category columns.
disease.state <- disease.state %>% mutate(disease_cat = ifelse(disease_state_short < 3, 
                                                               "good", 
                                                               ifelse(disease_state_short == 3, "middle", "bad")))
# Merge metadata; the original workflow notes that 251 samples remain. 
metadata.table <- merge(sample_data(ps) %>% data.frame(), disease.state, by.x="Sample", by.y= "sample_id", all.y=T, all.x=T) %>%
  filter(Sample != "") %>% filter(!is.na(disease_state_short))  %>% filter(!is.na(Replicate)) %>% 
  select(-c("Site", "baum_id")) 
rownames(metadata.table) <- metadata.table$Sample
# Load the processed phyloseq object used for microbiome feature generation.
ps <- readRDS("data/processed/phyloseq4.rds")
# Retain samples with at least 19,000 reads.
ps <- prune_samples(sample_sums(ps)>=19000, ps)
# Retain taxa meeting the prevalence and abundance thresholds defined below.
ps.filt <- ps %>%   filter_taxa(function(x) (sum(x >= 1) >= nsamples(.)* 0.1 ) && sum(x >= 3) > 10, TRUE)
ps.filt
sample_data(ps.filt) <- sample_data(metadata.table)
# Tree-based models do not require scaling for prediction, although scaling may aid feature-importance interpretation.

#rownames(rf.df.scale) <- rownames(rf.df)
#summary(rf.df.scale)

# Add microbiome data.
ps.genus <- ps %>% tax_glom("Genus")  
table.genus.hel <- ps.genus %>% filter_taxa(function(x) sum(x > 30) > (0.1*length(x)), TRUE) %>% 
  microbiome::transform("hellinger") %>% otu_table() %>% data.frame() 
# Rename ASVs using concatenated taxonomy strings.
tax2ASV <- tax_table(ps.genus) %>% data.frame() %>% 
  tidyr::unite(taxonomy, all_of(c("Phylum","Class", "Order", "Family", "Genus")), sep = ";") %>% 
  dplyr::select(taxonomy)
colnames(table.genus.hel) <- tax2ASV[colnames(table.genus.hel), "taxonomy"]


# Prepare filtered ASV-level data.
table.hel <- ps.filt %>%   microbiome::transform("hellinger") %>% otu_table() %>% data.frame() 
# Rename ASVs using concatenated taxonomy strings.
tax2ASV <- tax_table(ps.filt) %>% data.frame() 
tax2ASV$ASVID <-rownames(tax2ASV)
tax2ASV <- tax2ASV %>%
  tidyr::unite(taxonomy, all_of(c("Phylum","Class", "Order", "Family", "Genus", "ASVID")), sep = ";") %>% 
  dplyr::select( "taxonomy")
colnames(table.hel) <- tax2ASV[colnames(table.hel), "taxonomy"]


# Add vitality-assessment metadata to the phyloseq object.
meta <- ps %>% sample_data() %>% data.frame()
meta <- merge(meta, df[, c("Vitalit?t_Sommer", "LFE_Tendenz_Vitalit?t_Sommer", "LFE_Vitalit?t_Sommer", "Vitalit?t_Winter_danach" )] , all.x = T, by="row.names") %>%
  tibble::column_to_rownames("Row.names")
sample_data(ps) <- sample_data(meta)
