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
library(ggplot2)
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

# read in phylogeny trimmed to dataset (Stein et al.)
spine.all= read.csv(file="AllTaxa.csv")
row.names(spine.all)= spine.all$PhyloSpp

spine.all.log= spine.all
spine.all.log[,9:12]= log(spine.all[,9:12])

spine.data= spine.all.log[,c(6:12)]

# test traits for deviation from normality
shapiro.test(spine.data$SM50) 
shapiro.test(spine.data$IncludedAngle) 
shapiro.test(spine.data$SerrationAngle)
shapiro.test(spine.data$PropSerrated)
shapiro.test(spine.data$SerrationLength)
shapiro.test(spine.data$SA.Vol) 
shapiro.test(spine.data$SpineAR) 


##############################################################################
## READ IN TREE
##############################################################################
# Read in a single summary phylogenetic tree from the VertLife database 
# (Stein et al. 2018; https://vertlife.org/phylosubsets/).
# The summary tree was generated using TreeAnnotator (BEAST v1.10.5; Drummond & Rambaut, 2007)
# from the full distribution of 10,000 time-calibrated trees with 10 fossil calibrations.
# Ladderize for consistent plotting and visualization.
# Sources:
#   - Stein et al. 2018, *Nature Ecology & Evolution*, https://doi.org/10.1038/s41559-017-0448-4
#   - Drummond & Rambaut 2007, *BMC Evolutionary Biology*, https://doi.org/10.1186/1471-2148-7-214
##############################################################################

ray.tree= read.nexus("SummaryTree") # single tree
ray.tree= ladderize(ray.tree)
plot(ray.tree)

##############################################################################
## CREATE OBJECTS FOR TRAITS & CATEGORICAL DATA
##############################################################################
# Extract and name trait vectors (performance and morphology) and categorical 
# variables (habitat and family) for use in comparative analyses. 
# Continuous traits are scaled for multivariate tests.
# Custom color schemes are defined for habitat and family groups.
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
Fams= spine.all$Family
names(Fams)= spine.all$PhyloSpp

##Color schemes
cols.hab= setNames(c("steelblue1","springgreen4"), unique(Habs))
cols.fam= setNames(c("#820000","#4A4A4A","#6A8C7E","#ff497d","#380099"), unique(Fams))

##############################################################################
## STOCHASTIC CHARACTER MAPPING (HABITAT)
##############################################################################
# Performs stochastic character mapping of habitat (marine vs freshwater) 
# on the summary tree. 
# Compares symmetric (SYM) and all‑rates‑different (ARD) models of discrete character 
# evolution, then generates 20 stochastic maps under the preferred model. 
# One representative map is selected at random for plotting.
##############################################################################
set.seed(123)

fitSYM<-fitDiscrete(ray.tree,Habs,model="SYM") # preferred over ARD
print(fitSYM,digits=3)

fitARD<-fitDiscrete(ray.tree,Habs,model="ARD")
print(fitARD,digits=3)

simmap.hab= make.simmap(tree=ray.tree, x=Habs, model="SYM", nsim=20, pi="estimated")
simmap.hab.plot <- simmap.hab[[ sample(1:length(simmap.hab), 1) ]]

##############################################################################
## PCA OF MULTIVARIATE DATASET
##############################################################################
# Perform principal component analysis (PCA) on the full trait dataset.
# Extract PC scores and percent variance explained for plotting. 
# Combine PC scores with categorical variables (habitat, family) for 
# visualization and group comparisons.
##############################################################################

pca.funct= gm.prcomp(spine.data, scale=TRUE)

# extract % variance on PCs 1 and 2
pc1.var= round(((pca.funct$sdev[1]^2) / sum(pca.funct$sdev^2) * 100), 1) 
pc2.var= round(((pca.funct$sdev[2]^2) / sum(pca.funct$sdev^2) * 100), 1)

######## COMBINE CATEGORICAL TRAITS WITH PC TABLE ########

# table with first 4 PC axes
pc.data= as.data.frame(pca.funct$x[,1:2])
names(pc.data)= c("PC.1", "PC.2")
row.names(pc.data)= spine.all$PhyloSpp

# table with categorical data
grps= spine.all[,2:3]
names(grps)= c("Family","Habitat")

##############################################################################
## PHYLOMORPHOSPACE BY HABITAT
##############################################################################
# Plot phylomorphospace using PC1 and PC2 scores with a habitat-mapped tree.
# Add convex hulls for marine and freshwater groups to highlight habitat-specific regions.
##############################################################################

# Plot phylomorphospace
par(mar = c(5, 5, 2, 2))
plot(NA, xlim=range(pc.data$PC.1), ylim=range(pc.data$PC.2), asp=1, axes=FALSE, xlab=paste("PC 1 (",pc1.var,"%)", sep=""), ylab=paste("PC 2 (",pc2.var,"%)", sep="")) +
axis(1, at=seq(-2, 2, by=1), cex.axis=0.8) +
axis(2, at=seq(-2, 1, by=1), cex.axis=0.8)
abline(v=seq(-2, 3, by=1),col="grey", lty="dotted")
abline(h=seq(-2, 3, by=1),col="grey", lty="dotted")

phylomorphospace(tree=simmap.hab.plot, pc.data, bty="n", colors=cols.hab, lwd=4, ftype="off",node.by.map=TRUE, node.size=c(0,1.2), add=TRUE)

legend(x=-2.8, y=2.9, c("marine", "freshwater"), pch=21, pt.bg=cols.hab, bty="n", pt.cex=1.2, cex=0.8)

# define freshwater hull
fw.rows= which(grps$Habitat=="FW")
fw.hull= chull(pc.data[grps$Habitat=="FW",])

# define marine hull
mar.rows= which(grps$Habitat=="SW")
mar.hull= chull(pc.data[grps$Habitat=="SW",])

# plot hulls
polygon(pc.data[fw.rows,][fw.hull,], col= adjustcolor(cols.hab[2], alpha=0.2), border=F)
polygon(pc.data[mar.rows,][mar.hull,], col= adjustcolor(cols.hab[1], alpha=0.2), border=F)

##############################################################################
## PHYLOMORPHOSPACE BY FAMILY
##############################################################################
# Plot phylomorphospace using PC1 and PC2 scores with a family-mapped tree.
# Add convex hulls for groups with > 2 sampled species to highlight family-specific regions.
##############################################################################

# Combine PC data with Family information
pc.family.table <- cbind(pc.data, Family = grps$Family)

pc.family.table <- pc.family.table[match(ray.tree$tip.label, rownames(pc.family.table)), ]

fam <- factor(pc.family.table$Family) # Convert Family column to factor
names(Fams)= spine.all$PhyloSpp
tip.colors <- cols.fam[fam]

# Plot phylomorphospace
par(mar = c(5, 5, 2, 2))
phylomorphospace(ray.tree, pc.data, node.by.map = FALSE,
                 ftype = "off", node.size = c(0, 1), bty = "n", las = 1,
                 xlab = "PC1",
                 ylab = expression(paste("PC2")))

tiplabels(pch = 21, bg = tip.colors, cex = 1.2)

legend(x = par("usr")[1], y = par("usr")[4] + 0.5, # Adjust y higher with offset
       legend = levels(fam), cex = 0.8, pch = 21,
       pt.bg = cols.fam, pt.cex = 1.5, bty = "n", # Remove legend box
       xjust = 0, yjust = 1) # Align top left of the legend to these coordinates

# define dasyatid hull
D.rows= which(grps$Family=="Dasyatidae")
D.hull= chull(pc.data[grps$Family=="Dasyatidae",])

# define potamotrygonid hull
P.rows= which(grps$Family=="Potamotrygonidae")
P.hull= chull(pc.data[grps$Family=="Potamotrygonidae",])

# define urotrygonid hull
U.rows= which(grps$Family=="Urotrygonidae")
U.hull= chull(pc.data[grps$Family=="Urotrygonidae",])

#plot hulls
polygon(pc.data[D.rows,][D.hull,], col= adjustcolor(cols.fam[1], alpha=0.2), border=F)
polygon(pc.data[P.rows,][P.hull,], col= adjustcolor(cols.fam[3], alpha=0.2), border=F)
polygon(pc.data[U.rows,][U.hull,], col= adjustcolor(cols.fam[5], alpha=0.2), border=F)

##############################################################################
## BOXPLOT BY FAMILY
##############################################################################
# Generate a boxplot of a single trait (here, log-transformed aspect ratio) grouped by family.
# This section is modular and can be changed by assigning a different trait variable from spine.data.
##############################################################################

trait.plt= spine.all$Dorsoventral # add name of the trait you want to plot from "spine.data"
t=ggplot(grps, aes(x=Family, y=trait.plt, fill=Family))  +  
  geom_boxplot(alpha=0.6, outlier.color=NA) +
  geom_jitter(shape=21, width=0.25, size=2) +
  scale_fill_manual(values=cols.fam) + 
  labs(x="Family", y="Average Stress (Dorsoventral)")+  
  theme(panel.background=element_blank(),  
        title=element_text(family="sans"),
        axis.text.y=element_text(size=13, colour="black", family="sans"), 
        axis.text.x=element_text(size=11, colour="black", family="sans", angle=0),  
        axis.title.y=element_text(angle=90, size=14, colour="black", family="sans"), 
        axis.title.x=element_text(angle=0, size=14, colour="black", family="sans"),
        strip.background= element_blank(),
        strip.text.x= element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(colour="black", fill=NA, linewidth=1),
        axis.ticks=element_line(colour="black"),
        axis.line=element_blank(),
        legend.text=element_text(family="sans", size=10), 
        legend.title=element_blank(),
        legend.key=element_rect(fill="transparent"),
        legend.position="none"
  )
t
ggsave(filename="Dorsoventral_Plot.pdf", plot=t, device="pdf", width=8, height=6) #change name to match current trait

##############################################################################
## ESTIMATE SPINE MORPHOLOGIES ALONG PC1 FOR PHYSICAL MODELS
##############################################################################
# Estimate trait values across a range of PC1 scores to reconstruct hypothetical 
# spine morphologies along the primary axis of shape variation. 
# PC coordinates are back-transformed to trait space using PCA loadings, and logged 
# traits are returned to original scale. 
# Output is used for creating physical models.
##############################################################################

# specify coordinates along PC1 to estimate values (5 estimates)
PC1.step= (max(pc.data$PC.1) - min(pc.data$PC.1))/4 

x1= min(pc.data$PC.1) # min PC1 value
x2= x1 + PC1.step
x3= x2 + PC1.step
x4= x3 + PC1.step
x5= x4 + PC1.step # max PC1 value

sample.coords= as.data.frame(matrix(nrow=5, ncol=2))
sample.coords[,1]= c(x1, x2, x3, x4, x5)
sample.coords[,2]= c(0, 0, 0, 0, 0)
names(sample.coords)= c("PC.1","PC.2")

# Get the loadings (rotation matrix), center, and scale
loadings <- pca.funct$rotation
center <- pca.funct$center
scale <- pca.funct$scale

# constructing a data frame with estimated trait values along PC1 
PC1.est= as.data.frame(matrix(ncol=9, nrow=nrow(sample.coords)))
PC1.est[,1:2]= sample.coords
names(PC1.est)= c("PC.1", "PC.2", row.names(loadings))

for(i in 1:nrow(sample.coords)){
pc.coords= c(PC1=sample.coords[i,1], PC2=sample.coords[i,2], PC3=0, PC4=0, PC5=0, PC6=0, PC7=0)

# transform data back to original variable space
pc.pred= (pc.coords %*% t(loadings)) * scale + center 
PC1.est[i,3:9]= pc.pred
}

PC1.est[,6:9]= exp(PC1.est[,6:9]) # e^x, to convert logged traits to original scale

# table of estimated serration traits and included angle
spine.model.PC1= PC1.est[,c(1,2,4,5,6,7)]
write.csv(spine.model.PC1, "PC1 spine models.csv")

##############################################################################
## PHYLOGENETIC (M)ANOVAs
##############################################################################
# Tests for differences in spine morphology across habitats (marine vs freshwater) 
# while accounting for phylogenetic relationships. 
# A multivariate phylogenetic  MANOVA is performed on the scaled trait dataset (9999 iterations).
# Phylogenetic ANOVAS are performed for univatiate traits (9999 iterations).
##############################################################################

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
# Runs phylogenetic generalized least squares (PGLS) regressions on the summary tree 
# using a Brownian motion correlation structure. 
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
AR.Dorsoventral.pgls<-gls(Dorsoventral ~ SpineAR, correlation=corBM)
PC1.Dorsoventral.pgls<-gls(Dorsoventral ~ PC1.scr, correlation=corBM)
SM.Dorsoventral.pgls<-gls(Dorsoventral ~ SM50, correlation=corBM)
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
AR.Lateral.pgls<-gls(Lateral ~ SpineAR, correlation=corBM)
PC1.Lateral.pgls<-gls(Lateral ~ PC1.scr, correlation=corBM)
SM.Lateral.pgls<-gls(Lateral ~ SM50, correlation=corBM)
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

##############################################################################
### PGLS REGRESSION PLOTS
##############################################################################
# Generate plots for PGLS regressions.
# Plots include Dorsoventral and Lateral bending scenarios regressed against 
# both aspect ratio (SpineAR) and PC1 scores. 
# Regression lines are drawn from PGLS models, and figures are saved in the directory.
##############################################################################

### Plot 1: Dorsoventral vs AR
svg("Dorsoventral_vs_AR.svg", width=5, height=5)
par(mar=c(4,4,2,2))
buffer.x <- 0.1 * diff(range(SpineAR, na.rm = TRUE))
buffer.y <- 0.1 * diff(range(Dorsoventral, na.rm = TRUE))
plot(Dorsoventral ~ SpineAR, 
     pch=21, bg=palette()[1], cex=1.5, bty="n",
     xlab="log(aspect ratio)", ylab="Stress (Dorsoventral)",
     cex.lab=1, cex.axis=0.8, las=1,
     xlim=range(SpineAR, na.rm = TRUE) + c(-buffer.x, buffer.x),
     ylim=range(Dorsoventral, na.rm = TRUE) + c(-buffer.y, buffer.y))
abline(AR.Dorsoventral.pgls, lwd=2, col="#9E2E50", lty=2)
box(lwd=1)
dev.off()

### Plot 2: Dorsoventral vs PC1
svg("Dorsoventral_vs_PC1.svg", width=5, height=5)
par(mar=c(4,4,2,2))
buffer.x <- 0.1 * diff(range(PC1.scr, na.rm = TRUE))
buffer.y <- 0.1 * diff(range(Dorsoventral, na.rm = TRUE))
plot(Dorsoventral ~ PC1.scr, 
     pch=21, bg=palette()[1], cex=1.5, bty="n",
     xlab="PC1 - spine morphology", ylab="Stress (Dorsoventral)",
     cex.lab=1, cex.axis=0.8, las=1,
     xlim=range(PC1.scr, na.rm = TRUE) + c(-buffer.x, buffer.x),
     ylim=range(Dorsoventral, na.rm = TRUE) + c(-buffer.y, buffer.y))
abline(PC1.Dorsoventral.pgls, lwd=2, col="#9E2E50", lty=2)
box(lwd=1)
dev.off()

### Plot 3: Lateral vs AR
svg("Lateral_vs_AR.svg", width=5, height=5)
par(mar=c(4,4,2,2))
buffer.x <- 0.1 * diff(range(SpineAR, na.rm = TRUE))
buffer.y <- 0.1 * diff(range(Lateral, na.rm = TRUE))
plot(Lateral ~ SpineAR, 
     pch=21, bg=palette()[1], cex=1.5, bty="n",
     xlab="log(aspect ratio)", ylab="Stress (Lateral)",
     cex.lab=1, cex.axis=0.8, las=1,
     xlim=range(SpineAR, na.rm = TRUE) + c(-buffer.x, buffer.x),
     ylim=range(Lateral, na.rm = TRUE) + c(-buffer.y, buffer.y))
abline(AR.Lateral.pgls, lwd=2, col="#9E2E50", lty=2)
box(lwd=1)
dev.off()

### Plot 4: Lateral vs PC1
svg("Lateral_vs_PC1.svg", width=5, height=5)
par(mar=c(4,4,2,2))
buffer.x <- 0.1 * diff(range(PC1.scr, na.rm = TRUE))
buffer.y <- 0.1 * diff(range(Lateral, na.rm = TRUE))
plot(Lateral ~ PC1.scr, 
     pch=21, bg=palette()[1], cex=1.5, bty="n",
     xlab="PC1 - spine morphology", ylab="Stress (Lateral)",
     cex.lab=1, cex.axis=0.8, las=1,
     xlim=range(PC1.scr, na.rm = TRUE) + c(-buffer.x, buffer.x),
     ylim=range(Lateral, na.rm = TRUE) + c(-buffer.y, buffer.y))
abline(PC1.Lateral.pgls, lwd=2, col="#9E2E50", lty=2)
box(lwd=1)
dev.off()



##############################################################################
## LIKELIHOOD-BASED R2 CALCULATION FOR PGLS MODELS
##############################################################################
# Calculate likelihood-based R2 for PGLS (gls) models.
# Models are refit using ML and compared to intercept-only models
# with the same Brownian motion correlation structure.
# R2 values are computed for dorsoventral and lateral models.

pgls_R2_ML <- function(model) {
  model_ML <- update(model, method = "ML")
  null_ML  <- update(model_ML, . ~ 1)
  
  ll_full <- as.numeric(logLik(model_ML))
  ll_null <- as.numeric(logLik(null_ML))
  
  1 - exp((-2 / nobs(model_ML)) * (ll_full - ll_null))
}

R2_values <- c(
  AR.Dorsoventral  = pgls_R2_ML(AR.Dorsoventral.pgls),
  PC1.Dorsoventral = pgls_R2_ML(PC1.Dorsoventral.pgls),
  SA.Dorsoventral  = pgls_R2_ML(SA.Dorsoventral.pgls),
  PS.Dorsoventral  = pgls_R2_ML(PS.Dorsoventral.pgls),
  AR.Lateral       = pgls_R2_ML(AR.Lateral.pgls),
  PC1.Lateral      = pgls_R2_ML(PC1.Lateral.pgls),
  IA.Lateral       = pgls_R2_ML(IA.Lateral.pgls),
  SA.Lateral       = pgls_R2_ML(SA.Lateral.pgls)
)

print(R2_values)


sessionInfo()
