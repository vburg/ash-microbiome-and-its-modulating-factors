# =============================================================================
# Variance partitioning and partial dbRDA
#
# Purpose:
#   Partition variation in the microbial community among geographical,
#   climatic, soil, leaf, and tree-related explanatory variables and test
#   unique contributions with partial distance-based redundancy analyses.
#
# Required objects before running this script:
#   ps1         - phyloseq object containing the ASV table and sample metadata
#   ps          - phyloseq object referenced by taxa_are_rows()
#   impute_data - sample-level environmental and host metadata
#
# Required metadata columns in impute_data:
#   x_utm32, y_utm32, MAP, RHJul, Tjul, ca_vorrat, NA_vorrat, KAK,
#   GChl, PC, and dBHD_ashes
#
# Reproducibility requirements:
#   - Ensure that rows in impute_data are in exactly the same order as samples
#     in the community matrix before filtering or modeling.
#   - Record package versions with renv::snapshot() or sessionInfo().
#   - Use a fixed random seed before permutation tests when exact repeatability
#     is required, for example set.seed(1234).
#   - Report the number of samples retained after complete-case filtering.
#   - Document whether the input community table contains raw, filtered,
#     rarefied, or otherwise transformed counts.
#
# Important:
#   The executable analysis code, formulas, object names, and parameters below
#   are preserved exactly. Comments were translated and expanded only.
# =============================================================================

# Load required packages -----------------------------------------------------
library(phyloseq)
library(vegan)
library(dplyr)

# 1. Extract the ASV table ---------------------------------------------------
otu <- as(otu_table(ps1), "matrix")

# If taxa are stored as rows, transpose the matrix so samples become rows.
if (taxa_are_rows(ps)) {
  otu <- t(otu)
}

# 2. Apply a Hellinger transformation to the community data -----------------
otu_hel <- decostand(otu, method = "hellinger")

# 3. Extract sample metadata -------------------------------------------------
# Confirm that the resulting metadata rows correspond exactly to otu_hel rows.
meta <- as(sample_data(otu), "data.frame")

# 4. Define explanatory-variable groups -------------------------------------

# Geographical coordinates.
geo <- impute_data %>%
  select(x_utm32, y_utm32)

# Climatic variables.
climate <- impute_data %>%
  select(MAP, RHJul, Tjul)

# Soil variables.
soil <- impute_data %>%
  select(ca_vorrat, NA_vorrat, KAK)

# Leaf variables.
leaf <- impute_data %>%
  select(GChl, PC)

# Tree-related variable.
tree <- impute_data %>%
  select(dBHD_ashes)

# 5. Remove samples with missing values -------------------------------------
# complete_samples is TRUE only for samples complete across all variable groups.
complete_samples <- complete.cases(geo, climate, soil, leaf, tree)

otu_hel2 <- otu_hel[complete_samples, ]
geo2     <- geo[complete_samples, ]
climate2 <- climate[complete_samples, ]
soil2    <- soil[complete_samples, ]
leaf2    <- leaf[complete_samples, ]
tree2    <- tree[complete_samples, ]

# 6. Standardize explanatory variables -------------------------------------
# scale() centers each variable and divides it by its standard deviation.
geo2     <- scale(geo2)
climate2 <- scale(climate2)
soil2    <- scale(soil2)
leaf2    <- scale(leaf2)
tree2    <- scale(tree2)


# Combine geographical and climatic variables as site-level predictors.
site_env <- cbind(geo2, climate2)

# Combine leaf and tree variables as host-level predictors.
host_env <- cbind(leaf2, tree2)

# First variance-partitioning option ----------------------------------------
# Partition variation among site-level, soil, and host-level predictor groups.
vp <- varpart(
  otu_hel2,
  site_env,
  soil2,
  host_env
)

plot(vp)
vp

# Alternative variance-partitioning option ----------------------------------
# Partition variation among geography, climate, and combined local predictors.
vp <- varpart(
  otu_hel2,
  geo2,
  climate2,
  cbind(soil2, leaf2, tree2)
)

plot(vp)
vp


##

# Test geographical effects independently of climate and soil/host effects.
anova.cca(
  capscale(otu_hel2 ~ .,
           data = as.data.frame(geo2),
           add = TRUE),
  permutations = 999
)

###

# Preferred combined workflow -----------------------------------------------

# Combine all standardized explanatory variables into one data frame.
env_all <- data.frame(
  geo2,
  climate2,
  soil2,
  leaf2,
  tree2
)

# Fit the full constrained ordination model.
mod_all <- capscale(
  otu_hel2 ~ .,
  data = env_all
)

# Test the full model and each term using 999 permutations.
anova(mod_all, permutations = 999)
anova(mod_all, by = "terms", permutations = 999)

# Test geography after accounting for all remaining variables.
mod_geo <- capscale(
  otu_hel2 ~ x_utm32 + y_utm32 +
    Condition(MAP + RHJul + Tjul +
                ca_vorrat + NA_vorrat + KAK +
                GChl + PC + tree2),
  data = env_all,
  distance = "bray"
)

anova(mod_geo, permutations = 999)

# Repeat the partial dbRDA for each explanatory-variable group ---------------

# Test climate after accounting for geography, soil, leaf, and tree variables.
mod_climate <- capscale(
  otu_hel2 ~ MAP + RHJul + Tjul +
    Condition(x_utm32 + y_utm32 +
                ca_vorrat + NA_vorrat + KAK +
                GChl + PC + tree2),
  data = env_all,
  distance = "bray"
)

anova(mod_climate, permutations = 999)


# Test combined local soil, leaf, and tree effects after accounting for
# geography and climate.
mod_local_unique <- capscale(
  otu_hel2 ~ ca_vorrat + NA_vorrat + KAK +
    GChl + PC + tree2 +
    Condition(x_utm32 + y_utm32 +
                MAP + RHJul + Tjul),
  data = env_all
)

anova(mod_local_unique, permutations = 999)
