##############################################################################
## Code for "Stingray Spine Diversity Reflects Performance Trade-Offs 
## Linked to Puncture and Breakability"
## Author: Emily Poulin
##############################################################################

##############################################################################
## OPEN PACKAGES
##############################################################################
# Load required R packages for data analysis and visualization of physical 
# model puncture and removal performance.
##############################################################################
library(ggplot2)

################################################################################
## PERFORMANCE ANALYSES (PUNCTURE & REMOVAL)
################################################################################
# Analyze puncture and removal performance across physical spine models.
# Steps:
# - Read in dynamic and quasi-static puncture and removal data
# - Visualize performance across models
# - Perform ANOVA to test for differences among models
# - Run Tukey HSD post-hoc tests
# - Save ANOVA summaries (txt) and Tukey results (csv)
################################################################################

################################################################################
## DYNAMIC
################################################################################

# Set seed 
set.seed(123)

# Read data
dyn_removal_data <- read.csv(file = "DynamicRemoval.csv", header = TRUE)

# Make Model a factor
dyn_removal_data$Model <- factor(dyn_removal_data$Model, levels = 1:5)

# Plot
dyn_removal_plot <- ggplot(dyn_removal_data, aes(x = Model, y = Removal_Force, fill = Model)) +
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "grey20") +
  labs(x = "Model", y = "Removal Force (N)", fill = "Model") +
  ggtitle("Removal Force Variation among Models") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#BFE6F2", "#8FB3D9", "#6E73B8", "#4E55A8", "#1F3F9B"))
plot(dyn_removal_plot)

# ANOVA
dyn_removal_aov <- aov(Removal_Force ~ Model, data = dyn_removal_data)
sink("Dynamic_Removal_Force_ANOVA_summary.txt")
summary(dyn_removal_aov)
sink()

# Tukey HSD
dyn_removal_tukey <- TukeyHSD(dyn_removal_aov)

# Format Tukey results
dyn_removal_tukey_df <- as.data.frame(dyn_removal_tukey$Model)
dyn_removal_tukey_df$Comparison <- rownames(dyn_removal_tukey_df)
dyn_removal_tukey_df <- dyn_removal_tukey_df[, c("Comparison", "diff", "lwr", "upr", "p adj")]
colnames(dyn_removal_tukey_df) <- c("Comparison", "Difference", "Lower_CI", "Upper_CI", "P_value")
write.csv(dyn_removal_tukey_df,
          file = "Dynamic_Removal_Force_TukeyHSD.csv",
          row.names = FALSE)

# Read data
dyn_puncture_data <- read.csv(file = "DynamicPuncture.csv", header = TRUE)

# Make Model a factor
dyn_puncture_data$Model <- factor(dyn_puncture_data$Model, levels = 1:5)

# Plot
dyn_puncture_plot <- ggplot(dyn_puncture_data,
                            aes(x = Model, y = Puncture_Distance, fill = Model)) +
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "grey20") +
  labs(x = "Model", y = "Puncture Distance (mm)", fill = "Model") +
  ggtitle("Puncture Distance Variation among Models") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#BFE6F2", "#8FB3D9", "#6E73B8", "#4E55A8", "#1F3F9B"))
plot(dyn_puncture_plot)

# ANOVA
dyn_puncture_aov <- aov(Puncture_Distance ~ Model, data = dyn_puncture_data)
sink("Dynamic_Puncture_ANOVA_summary.txt")
summary(dyn_puncture_aov)
sink()

# Tukey HSD
dyn_puncture_tukey <- TukeyHSD(dyn_puncture_aov)

# Format Tukey results
dyn_puncture_tukey_df <- as.data.frame(dyn_puncture_tukey$Model)
dyn_puncture_tukey_df$Comparison <- rownames(dyn_puncture_tukey_df)
dyn_puncture_tukey_df <- dyn_puncture_tukey_df[, c("Comparison", "diff", "lwr", "upr", "p adj")]
colnames(dyn_puncture_tukey_df) <- c("Comparison", "Difference", "Lower_CI", "Upper_CI", "P_value")
write.csv(dyn_puncture_tukey_df,
          file = "Dynamic_Puncture_TukeyHSD.csv",
          row.names = FALSE)

################################################################################
## QUASI-STATIC
################################################################################

slow_removal_data <- read.csv("SlowRemoval.csv", header = TRUE)

# Make Model a factor
slow_removal_data$Model <- factor(slow_removal_data$Model, levels = 1:5)

shapiro.test(slow_removal_data$Min_Force)

# Make all forces positive
slow_removal_data$Min_Force <- abs(slow_removal_data$Min_Force)

slow_removal_plot <- ggplot(slow_removal_data,
                            aes(x = Model, y = Min_Force, fill = Model)) +
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "grey20") +
  labs(x = "Model", y = "Removal Force (N)", fill = "Model") +
  ggtitle("Minimum Force Variation among Models") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#BFE6F2", "#8FB3D9", "#6E73B8", "#4E55A8", "#1F3F9B"))
plot(slow_removal_plot)

# ANOVA
slow_removal_aov <- aov(Min_Force ~ Model, data = slow_removal_data)
sink("Slow_Removal_ANOVA_summary.txt")
summary(slow_removal_aov)
sink()

# Tukey HSD
slow_removal_tukey <- TukeyHSD(slow_removal_aov)

# Format Tukey results
slow_removal_tukey_df <- as.data.frame(slow_removal_tukey$Model)
slow_removal_tukey_df$Comparison <- rownames(slow_removal_tukey_df)
slow_removal_tukey_df <- slow_removal_tukey_df[, c("Comparison", "diff", "lwr", "upr", "p adj")]
colnames(slow_removal_tukey_df) <- c("Comparison", "Difference", "Lower_CI", "Upper_CI", "P_value")
write.csv(slow_removal_tukey_df,
          file = "Slow_Removal_TukeyHSD.csv",
          row.names = FALSE)

# Read data
slow_puncture_data <- read.csv(file = "SlowPuncture.csv", header = TRUE)

# Make Model a factor
slow_puncture_data$Model <- factor(slow_puncture_data$Model, levels = 1:5)

shapiro.test(slow_puncture_data$Max_Force)

# Plot
slow_puncture_plot <- ggplot(slow_puncture_data,
                             aes(x = Model, y = Max_Force, fill = Model)) +
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "grey20") +
  labs(x = "Model", y = "Puncture Force (N)", fill = "Model") +
  ggtitle("Maximum Force Variation among Models") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#BFE6F2", "#8FB3D9", "#6E73B8", "#4E55A8", "#1F3F9B"))
plot(slow_puncture_plot)

# ANOVA
slow_puncture_aov <- aov(Max_Force ~ Model, data = slow_puncture_data)
sink("Slow_Puncture_ANOVA_summary.txt")
summary(slow_puncture_aov)
sink()

# Tukey HSD
slow_puncture_tukey <- TukeyHSD(slow_puncture_aov)

# Format Tukey results
slow_puncture_tukey_df <- as.data.frame(slow_puncture_tukey$Model)
slow_puncture_tukey_df$Comparison <- rownames(slow_puncture_tukey_df)
slow_puncture_tukey_df <- slow_puncture_tukey_df[, c("Comparison", "diff", "lwr", "upr", "p adj")]
colnames(slow_puncture_tukey_df) <- c("Comparison", "Difference", "Lower_CI", "Upper_CI", "P_value")
write.csv(slow_puncture_tukey_df,
          file = "Slow_Puncture_TukeyHSD.csv",
          row.names = FALSE)

sessionInfo()
