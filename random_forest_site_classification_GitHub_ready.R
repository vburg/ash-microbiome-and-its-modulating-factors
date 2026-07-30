# =============================================================================
# Random-forest classification of sampling sites using genus-level abundances
#
# Purpose:
#   Aggregate ASVs to genus level, convert counts to relative abundances,
#   filter genera by mean abundance and prevalence, classify sampling sites
#   with a random forest, and extract genus-level variable importance.
#
# Required object before running this script:
#   ps1 - phyloseq object containing:
#         - an abundance table
#         - taxonomy including the rank "Genus"
#         - sample metadata including the variable "Site"
#
# Reproducibility notes:
#   - Record package versions with renv::snapshot() or sessionInfo().
#   - The random seed is fixed at 123 before fitting the model.
#   - Document the number of samples and genera before and after filtering.
#   - Confirm that sample metadata and abundance-table rows remain aligned.
#   - The filtering thresholds retain genera with:
#       mean relative abundance > 0.001
#       prevalence in at least 10% of samples
#
# Important:
#   The executable code, parameters, object names, and workflow below are
#   preserved exactly. Only explanatory comments were added.
# =============================================================================

# Load phyloseq --------------------------------------------------------------
library(phyloseq)

# Aggregate ASVs at genus level ----------------------------------------------
# NArm = FALSE retains taxa with missing genus assignments during agglomeration.
ps_genus <- tax_glom(
  ps1,
  taxrank = "Genus",
  NArm = FALSE
)

# Inspect available taxonomic ranks and the number of genus-level taxa.
rank_names(ps_genus)

ntaxa(ps_genus)


# Convert counts to relative abundance within each sample --------------------
ps_genus <- transform_sample_counts(
  ps_genus,
  function(x) x / sum(x)
)

# Extract the genus-level abundance table as a matrix ------------------------
otu_genus <- as(otu_table(ps_genus), "matrix")

# Ensure that samples are rows and genera are columns.
if (taxa_are_rows(ps_genus)) {
  otu_genus <- t(otu_genus)
}

# Convert the matrix to a data frame for downstream filtering and modeling.
otu_genus <- as.data.frame(otu_genus)

# Extract the taxonomy table.
tax <- tax_table(ps_genus)

# Replace ASV/taxon identifiers with genus names as column names.
colnames(otu_genus) <- as.character(tax[, "Genus"])

# Inspect the resulting genus names.
head(colnames(otu_genus))


# Count duplicated genus names after renaming.
# Duplicates should be checked because duplicated predictor names can make
# interpretation of random-forest importance values ambiguous.
sum(duplicated(colnames(otu_genus)))


# Filter genera by abundance and prevalence ---------------------------------
# Retain genera with mean relative abundance greater than 0.1% and presence
# in at least 10% of samples.
keep <- colMeans(otu_genus) > 0.001 &
  colSums(otu_genus > 0) / nrow(otu_genus) >= 0.10

otu_rf <- otu_genus[, keep]

# Report the number of genera retained for random-forest classification.
ncol(otu_rf)


# Extract sample metadata and define the response variable -------------------
meta <- as(sample_data(ps_genus), "data.frame")

# Convert sampling site to a factor for classification.
site <- factor(meta$Site)

# Load the randomForest package ----------------------------------------------
library(randomForest)

# Fix the random seed for reproducibility.
set.seed(123)

# Fit the random-forest classifier ------------------------------------------
# ntree = 1000 specifies the number of trees.
# importance = TRUE calculates predictor-importance measures.
rf_site <- randomForest(
  x = otu_rf,
  y = site,
  ntree = 1000,
  importance = TRUE
)

# Print the fitted model, including the confusion matrix and OOB error.
rf_site


# Extract variable-importance measures --------------------------------------
imp <- importance(rf_site)

# Convert importance values into a data frame.
imp_df <- data.frame(
  Genus = rownames(imp),
  MeanDecreaseAccuracy = imp[, "MeanDecreaseAccuracy"],
  MeanDecreaseGini = imp[, "MeanDecreaseGini"]
)

# Rank genera by Mean Decrease Accuracy.
imp_df <- imp_df |>
  dplyr::arrange(desc(MeanDecreaseAccuracy))

# Display the 20 most important genera.
head(imp_df, 20)


# Inspect original taxon identifiers and genus assignments -------------------
taxa_names(ps_genus)[1:10]

head(tax_table(ps_genus)[, "Genus"])
