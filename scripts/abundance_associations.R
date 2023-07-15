# Set the seed for reproducibility
set.seed(531)

# Load necessary libraries
library(dplyr)
library(stats)
library(car)
library(ggplot2)
library(tidyr)
library(microbiome)
library(Hmisc)
library(magrittr)

# Load data
load("p5_NPT.QMP_138_species")
load("p5_NPT.RMP_138_species")

# Function to perform correlation analysis
perform_correlation_analysis <- function(data, covariates) {
  CMB_df <- data.frame(sample_data(data)[, covariates])
  species <- data.frame(otu_table(data))
  
  # Calculate correlations and adjust p-values
  species_corr_P <- associate(CMB_df, species, method = "spearman", mode = "table", p.adj.method = "none")
  species_corr_temp <- associate(CMB_df, species, method = "spearman", mode = "table", p.adj.method = "BH")
  species_corr_adjP <- cbind(species_corr_P, species_corr_temp$p.adj)
  setnames(species_corr_adjP, c("X1", "X2", "Correlation", "P", "p.adj"))
  
  return(head(species_corr_adjP))
}

# Perform correlation analysis for QMP species
qs1 <- p5_NPT.QMP_138_species
QMP_species_corr_adjP <- perform_correlation_analysis(qs1, c("Calprotectin", "Moisture", "BMI"))
head(QMP_species_corr_adjP)

# Perform correlation analysis for RMP species
ps1 <- p5_NPT.RMP_138_species
RMP_species_corr_adjP <- perform_correlation_analysis(ps1, c("Calprotectin", "Moisture", "BMI"))
head(RMP_species_corr_adjP)
