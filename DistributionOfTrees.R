##############################################################################
## Code for "Stingray Spine Diversity Reflects Performance Trade-Offs 
## Linked to Puncture and Breakability"
## Author: Emily Poulin
##############################################################################

##############################################################################
## OPEN PACKAGES & SET WORKING DIRECTORY
##############################################################################
# Load required R packages for phylogenetic and morphological data analysis.
##############################################################################
library(phytools)
library(dplyr)
library(nlme)
library(ape)

##############################################################################
## READ IN DATA
##############################################################################
# Read in morphological trait dataset and assign species names as row names.
# Log-transform ratio-based traits (columns 9 to 12) to normalize data.
# Subset relevant trait columns for downstream analyses (columns 6 to 12).
##############################################################################
# Trait dataset
spine.all= read.csv(file="AllTaxa.csv")
row.names(spine.all)= spine.all$PhyloSpp

spine.all.log= spine.all
spine.all.log[,9:12]= log(spine.all[,9:12])

spine.data= spine.all.log[,c(6:12)]

##############################################################################
## READ IN TREES
##############################################################################
# Read in the full distribution of pruned phylogenetic trees to account for phylogenetic uncertainty.
# Randomly sample 500 trees (without replacement) for use in downstream analyses.
##############################################################################
#pruned.trees <- read.nexus("DistributionOfTrees.nex")

# Sample subset of trees
#subset_trees <- sample(pruned.trees, 500, replace = FALSE)

#write.nexus(subset_trees, file = "subset_trees.nex")

#read in subset
subset_trees <- read.nexus("Subset_trees.nex")

##############################################################################
## CREATE OBJECTS FOR TRAITS & CATEGORICAL DATA
##############################################################################
# Extract and name trait vectors (performance and morphology) and the categorical 
# variable (habitat) for use in comparative analyses. 
# Continuous traits are scaled for multivariate tests.
##############################################################################
# Univariate traits
Lateral= spine.all.log$Lateral 
names(Lateral)= spine.all.log$PhyloSpp
Dorsoventral= spine.all.log$Dorsoventral 
names(Dorsoventral)= spine.all.log$PhyloSpp
SM50= spine.data$SM50
names(SM50)= spine.all$PhyloSpp
IncludedAngle= spine.data$IncludedAngle
names(IncludedAngle)= spine.all$PhyloSpp
SerrationAngle= spine.data$SerrationAngle
names(SerrationAngle)= spine.all$PhyloSpp
PropSerrated= spine.data$PropSerrated
names(PropSerrated)= spine.all$PhyloSpp
SerrationLength= spine.data$SerrationLength
names(SerrationLength)= spine.all$PhyloSpp
SA.Vol= spine.data$SA.Vol
names(SA.Vol)= spine.all$PhyloSpp
SpineAR= spine.data$SpineAR
names(SpineAR)= spine.all$PhyloSpp

# Multivariate traits
spine.data.scaled <- as.data.frame(scale(spine.data))

# Categorical traits
Habs= spine.all$Habitat
names(Habs)= spine.all$PhyloSpp

##############################################################################
## PCA OF MULTIVARIATE DATASET
##############################################################################
# Perform principal component analysis (PCA) on the full trait dataset
##############################################################################

#### PCA
pca.funct= gm.prcomp(spine.data, scale=TRUE)

##############################################################################
## MULTIVARIATE PHYLOGENETIC MANOVA ACROSS TREE DISTRIBUTION
##############################################################################
# Runs multivariate phylogenetic MANOVAs (procD.pgls) across 500 sampled trees
# Extracts R², F, Z, and P-values per tree, and summarizes significance outcomes.
##############################################################################

#Set seed
set.seed(123)

# Define analysis function for one tree
run_procD_pgls <- function(tree) {
  res <- procD.pgls(as.matrix(spine.data.scaled) ~ Habs, phy = tree, iter = 9999)
  
  # Extract statistics for the Habitat effect
  R2 <- res$aov.table$Rsq[1]
  F_value <- res$aov.table$F[1]
  Z_score <- res$aov.table$Z[1]
  P_value <- res$aov.table$`Pr(>F)`[1]
  
  return(c(R2 = R2, F = F_value, Z = Z_score, P = P_value))
}

# Run analysis across all sampled trees
results_pgls <- sapply(subset_trees, run_procD_pgls)

# Organize results into a data frame
results_df <- data.frame(
  tree_index = seq_along(subset_trees),
  t(results_pgls)
)

# Count P-values above and below threshold
p_breakdown <- data.frame(
  P_greater_0.05 = sum(results_df$P > 0.05, na.rm = TRUE),
  P_less_0.05    = sum(results_df$P <= 0.05, na.rm = TRUE)
)

# Save outputs
write.csv(results_df, file = "MultivariateFunctMorph_detailed.csv", row.names = FALSE)
write.csv(p_breakdown, file = "MultivariateFunctMorph_P_Breakdown.csv", row.names = FALSE)

##############################################################################
## UNIVARIATE PHYLOGENETIC ANOVAs ACROSS TREE DISTRIBUTION
##############################################################################
## # Runs phylogenetic ANOVAs (phylANOVA) across 500 sampled trees for univariate traits. 
# Extracts F and P-values per tree, and summarizes significance outcomes by trait for each tree
##############################################################################

# Define list of traits
trait_list <- list(
  SM50 = SM50,
  IncludedAngle = IncludedAngle,
  SerrationAngle = SerrationAngle,
  PropSerrated = PropSerrated,
  SerrationLength = SerrationLength,
  SA.Vol = SA.Vol,
  SpineAR = SpineAR
)

# Define analysis function for one trait on one tree
run_phylANOVA <- function(tree, trait, index) {
  current_trait <- trait[tree$tip.label]
  res <- phylANOVA(tree, Habs, current_trait, nsim = 9999, posthoc = FALSE)
  return(c(F = res$F, P = res$Pf))
}

# Run analyses across all traits and trees
results_trait <- lapply(names(trait_list), function(trait_name) {
  trait <- trait_list[[trait_name]]
  
  results <- sapply(seq_along(subset_trees), function(i) {
    run_phylANOVA(subset_trees[[i]], trait, i)
  })
  
  df <- data.frame(
    tree_index = seq_along(subset_trees),
    F = results["F", ],
    P = results["P", ],
    Trait = trait_name
  )
  return(df)
})

# Combine results into a single data frame
all_results_df <- bind_rows(results_trait)

# Save detailed results
write.csv(all_results_df,
          file = "UnivariateFunctMorph_detailed.csv",
          row.names = FALSE)

# Compute P-value breakdown per trait
p_breakdown <- all_results_df %>%
  group_by(Trait) %>%
  summarise(
    P_greater_0.05 = sum(P > 0.05, na.rm = TRUE),
    P_less_0.05    = sum(P <= 0.05, na.rm = TRUE)
  )

# Save P-value breakdown
write.csv(p_breakdown,
          file = "UnivariateFunctMorph_P_Breakdown.csv",
          row.names = FALSE)

##############################################################################
## PGLS REGRESSIONS WITH PERFORMANCE ACROSS TREE DISTRIBUTION
##############################################################################
# Runs phylogenetic generalized least squares regressions (PGLS) across 500 sampled 
# trees using a Brownian motion correlation structure. 
# Test relationships between performance metrics from two bending scenarios (Dorsoventral, Lateral)
# and morphological predictors (univariate traits + PC axes). 
# Extracts and saves p-values and summarizes significance outcomes by trait for each tree
##############################################################################

# Extract PC scores from PCA
PC1.scr = pca.funct$x[,1]
names(PC1.scr) = rownames(spine.data)

PC2.scr= pca.funct$x[,2]
names(PC2.scr)= rownames(spine.data)

PC3.scr= pca.funct$x[,3]
names(PC3.scr)= rownames(spine.data)

PC4.scr= pca.funct$x[,4]
names(PC4.scr)= rownames(spine.data)

PC5.scr= pca.funct$x[,5]
names(PC5.scr)= rownames(spine.data)

PC6.scr= pca.funct$x[,6]
names(PC6.scr)= rownames(spine.data)

PC7.scr= pca.funct$x[,7]
names(PC7.scr)= rownames(spine.data)

# Define predictors
predictors <- list(
  SM50 = SM50,
  IncludedAngle = IncludedAngle,
  SerrationAngle = SerrationAngle,
  PropSerrated = PropSerrated,
  SerrationLength = SerrationLength,
  SA.Vol = SA.Vol,
  SpineAR = SpineAR,
  PC1 = PC1.scr,
  PC2 = PC2.scr,
  PC3 = PC3.scr,
  PC4 = PC4.scr,
  PC5 = PC5.scr,
  PC6 = PC6.scr,
  PC7 = PC7.scr
)

# Responses
responses <- list(
  Dorsoventral = Dorsoventral,
  Lateral = Lateral
)

# Collect results
pvalue_results <- list()

for (i in seq_along(subset_trees)) {
  current_tree <- subset_trees[[i]]
  spp <- rownames(spine.data)
  corBM <- corBrownian(phy = current_tree, form = ~spp)
  
  for (resp_name in names(responses)) {
    resp <- responses[[resp_name]]
    
    for (pred_name in names(predictors)) {
      pred <- predictors[[pred_name]]
      
      # Build temporary data frame 
      dat <- data.frame(resp = resp, pred = pred, spp = spp)
      
      model <- gls(resp ~ pred, correlation = corBM, data = dat)
      pval <- summary(model)$tTable[2, "p-value"]
      
      pvalue_results[[length(pvalue_results) + 1]] <- data.frame(
        Iteration = i,
        Response = resp_name,
        Predictor = pred_name,
        P = pval
      )
    }
  }
}

# Combine into one data frame
pvalues_df <- bind_rows(pvalue_results)


write.csv(pvalues_df,
          file = "PGLS_pvalues_detailed.csv",
          row.names = FALSE)

# P-value breakdown
p_breakdown <- pvalues_df %>%
  group_by(Response, Predictor) %>%
  summarise(
    P_greater_0.05 = sum(P > 0.05, na.rm = TRUE),
    P_less_0.05    = sum(P <= 0.05, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(p_breakdown,
          file = "PGLS_pvalues_breakdown.csv",
          row.names = FALSE)

sessionInfo()
