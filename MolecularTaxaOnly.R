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
library(geomorph)
library(geiger)
library(phytools)
library(nlme)
library(ape)

##############################################################################
## READ IN DATA
##############################################################################
# Read in morphological trait dataset and assign species names as row names.
# Log-transform ratio-based traits (columns 9 to 12) to normalize data.
# Subset relevant trait columns for downstream analyses (columns 6 to 12).
##############################################################################

# read in phylogeny trimmed to dataset (Stein et al.)
spine.all= read.csv(file="MolecularTaxaOnly.csv")
row.names(spine.all)= spine.all$PhyloSpp

spine.all.log= spine.all
spine.all.log[,9:12]= log(spine.all[,9:12])

spine.data= spine.all.log[,c(6:12)]

##############################################################################
## READ IN TREE
##############################################################################
# Read in a single molecular phylogenetic tree from the VertLife database 
# (Stein et al. 2018; https://vertlife.org/phylosubsets/).
# The tree has been pruned to match species present in the trait dataset.
# Source: Stein et al. 2018, Nature Ecology & Evolution. https://doi.org/10.1038/s41559-017-0448-4
##############################################################################

ray.tree= read.tree("MolecularTaxaOnly.tre") # single tree

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
## PHYLOGENETIC (M)ANOVAs
##############################################################################
# Tests for differences in spine morphology across habitats (marine vs freshwater) 
# while accounting for phylogenetic relationships. 
# A multivariate phylogenetic  MANOVA is performed on the scaled trait dataset (9999 iterations).
# Phylogenetic ANOVAS are performed for univatiate traits (9999 iterations).
##############################################################################

set.seed(123)

## Multivariate functional morphology
habsAOV.funct=procD.pgls(as.matrix(spine.data.scaled) ~ Habs, phy=ray.tree, iter=9999)

## Univariate functional morphology traits
phylANOVA(ray.tree, Habs, SM50, nsim=9999, posthoc=FALSE)
phylANOVA(ray.tree, Habs, IncludedAngle, nsim=9999, posthoc=FALSE) 
phylANOVA(ray.tree, Habs, SerrationAngle, nsim=9999, posthoc=FALSE) 
phylANOVA(ray.tree, Habs, PropSerrated, nsim=9999, posthoc=FALSE) 
phylANOVA(ray.tree, Habs, SerrationLength, nsim=9999, posthoc=FALSE)
phylANOVA(ray.tree, Habs, SA.Vol, nsim=9999, posthoc=FALSE) 
phylANOVA(ray.tree, Habs, SpineAR, nsim=9999, posthoc=FALSE) 

##############################################################################
## VARIANCE TESTS
##############################################################################
# Tests for differences in variance of stingray spine morphology across habitats 
# (marine vs freshwater). 
# A multivariate disparity test is run on the scaled trait dataset (9,999 iterations). 
# Univariate disparity tests are run for each individual trait (9,999 iterations).
##############################################################################

## Multivariate functional morphology
habsvariance= morphol.disparity(as.matrix(spine.data.scaled) ~ Habs, groups=Habs, iter=9999)

## Univariate functional morphology traits
morphol.disparity(spine.data$SM50 ~ Habs, groups=Habs, iter=9999)
morphol.disparity(spine.data$IncludedAngle ~ Habs, groups=Habs, iter=9999)
morphol.disparity(spine.data$SerrationAngle ~ Habs, groups=Habs, iter=9999)
morphol.disparity(spine.data$PropSerrated~ Habs, groups=Habs, iter=9999)
morphol.disparity(spine.data$SerrationLength ~ Habs, groups=Habs, iter=9999)
morphol.disparity(spine.data$SA.Vol ~ Habs, groups=Habs, iter=9999)
morphol.disparity(spine.data$SpineAR ~ Habs, groups=Habs, iter=9999)


##############################################################################
## PGLS REGRESSIONS WITH PERFORMANCE 
##############################################################################
# Runs phylogenetic generalized least squares (PGLS) regressions using a Brownian 
# motion correlation structure. 
# Tests associations between spine morphology (univariate traits and PC axes) and 
# performance metrics from two bending scenarios (Dorsoventral and Lateral).
##############################################################################

# Set up Brownian Motion correlation matrix
spp= rownames(spine.data)
corBM= corBrownian(phy=ray.tree, form=~spp)

PC1.scr= pca.funct$x[,1]
names(PC1.scr)= rownames(spine.data)

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

## Dorsoventral
SM50.Dorsoventral.pgls<-gls(Dorsoventral ~ SM50, correlation=corBM)
AR.Dorsoventral.pgls<-gls(Dorsoventral ~ SpineAR, correlation=corBM)
PC1.Dorsoventral.pgls<-gls(Dorsoventral ~ PC1.scr, correlation=corBM)
IA.Dorsoventral.pgls<-gls(Dorsoventral ~ IncludedAngle, correlation=corBM)
SA.Dorsoventral.pgls<-gls(Dorsoventral ~ SerrationAngle, correlation=corBM)
PS.Dorsoventral.pgls<-gls(Dorsoventral ~ PropSerrated, correlation=corBM)
SL.Dorsoventral.pgls<-gls(Dorsoventral ~ SerrationLength, correlation=corBM)
SA.Vol.Dorsoventral.pgls<-gls(Dorsoventral ~ SA.Vol, correlation=corBM)
PC2.Dorsoventral.pgls<-gls(Dorsoventral ~ PC2.scr, correlation=corBM)
PC3.Dorsoventral.pgls<-gls(Dorsoventral ~ PC3.scr, correlation=corBM)
PC4.Dorsoventral.pgls<-gls(Dorsoventral ~ PC4.scr, correlation=corBM)
PC5.Dorsoventral.pgls<-gls(Dorsoventral ~ PC5.scr, correlation=corBM)
PC6.Dorsoventral.pgls<-gls(Dorsoventral ~ PC6.scr, correlation=corBM)
PC7.Dorsoventral.pgls<-gls(Dorsoventral ~ PC7.scr, correlation=corBM)

## Lateral
SM50.Lateral.pgls<-gls(Lateral ~ SM50, correlation=corBM)
AR.Lateral.pgls<-gls(Lateral ~ SpineAR, correlation=corBM)
PC1.Lateral.pgls<-gls(Lateral ~ PC1.scr, correlation=corBM)
IA.Lateral.pgls<-gls(Lateral ~ IncludedAngle, correlation=corBM)
SA.Lateral.pgls<-gls(Lateral ~ SerrationAngle, correlation=corBM)
PS.Lateral.pgls<-gls(Lateral ~ PropSerrated, correlation=corBM)
SL.Lateral.pgls<-gls(Lateral ~ SerrationLength, correlation=corBM)
SA.Vol.Lateral.pgls<-gls(Lateral ~ SA.Vol, correlation=corBM)
PC2.Lateral.pgls<-gls(Lateral ~ PC2.scr, correlation=corBM)
PC3.Lateral.pgls<-gls(Lateral ~ PC3.scr, correlation=corBM)
PC4.Lateral.pgls<-gls(Lateral ~ PC4.scr, correlation=corBM)
PC5.Lateral.pgls<-gls(Lateral ~ PC5.scr, correlation=corBM)
PC6.Lateral.pgls<-gls(Lateral ~ PC6.scr, correlation=corBM)
PC7.Lateral.pgls<-gls(Lateral ~ PC7.scr, correlation=corBM)

sessionInfo()
