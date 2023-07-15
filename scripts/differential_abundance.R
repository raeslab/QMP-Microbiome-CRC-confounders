# Set the seed for reproducibility
set.seed(531)

library(rstatix)
library(data.table)
library(tidyr)
library(dplyr)
library(microbiome)
library(vegan)
library(pairwiseAdonis)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(phyloseq)
library(ComplexHeatmap)
library(gridExtra)
library(scales)
library(vcd)
# Set the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary data
load("LCPM_runs_589_seq_new_SLV_nr99_v138.1_filtered")
load("QMP_ASV_LCPM_SLV_species")
load("RMP_ASV_LCPM_SLV_species")
load("RMP_ASV_LCPM_SLV_species_RA")
# Load metadata
LCMP_metadata_589 <- read.csv("LCMP_metadata_589.csv", row.names = 1)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Differencial abundance tests %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RMP @ Species SLV  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Assign variables
ps_sp <- RMP_ASV_LCPM_SLV_species_RA
df_names <- data.frame(otu_table(ps_sp))

# Remove taxa with 'unclassified' in their names
ps_sp_noUnk <- subset_taxa(ps_sp, !taxa_names(ps_sp) %in% taxa_names(ps_sp)[grep('.unclassified', taxa_names(ps_sp))])

# Remove taxa with '[' in their names
ps_sp_noUnk <- subset_taxa(ps_sp_noUnk, !taxa_names(ps_sp_noUnk) %in% taxa_names(ps_sp_noUnk)[grep("\\[", taxa_names(ps_sp_noUnk), fixed = TRUE)])

# Remove taxa with 'group.' in their names
ps_sp_noUnk <- subset_taxa(ps_sp_noUnk, !taxa_names(ps_sp_noUnk) %in% taxa_names(ps_sp_noUnk)[grep("group.", taxa_names(ps_sp_noUnk))])

# Modify taxa names by removing unnecessary characters
taxa_names(ps_sp_noUnk) <- gsub(".*\\g_", "", taxa_names(ps_sp_noUnk))
tax_table(ps_sp_noUnk)[, 7] <- gsub(".*\\g_", "", tax_table(ps_sp_noUnk)[, 7])

# Subset samples based on CRC_general_status values
LCPM_No_lesion <- subset_samples(ps_sp_noUnk, CRC_general_status == "No_lesion")
LCPM_Polyp <- subset_samples(ps_sp_noUnk, CRC_general_status == "Polyp")
LCPM_Tumor <- subset_samples(ps_sp_noUnk, CRC_general_status == "Tumor")

# Prune taxa based on prevalence and detection thresholds
p5_N <- prune_taxa(core_members(LCPM_No_lesion, prevalence = 5/100, detection = 1/1000000), LCPM_No_lesion)
p5_P <- prune_taxa(core_members(LCPM_Polyp, prevalence = 5/100, detection = 1/1000000), LCPM_Polyp)
p5_T <- prune_taxa(core_members(LCPM_Tumor, prevalence = 5/100, detection = 1/1000000), LCPM_Tumor)

# Perform set differences between taxa names
setdiff(taxa_names(p5_N), taxa_names(p5_P))
setdiff(taxa_names(p5_N), taxa_names(p5_T))
setdiff(taxa_names(p5_T), taxa_names(p5_N))
setdiff(taxa_names(p5_P), taxa_names(p5_T))

# Subset taxa based on the combined set differences
p5_NPT <- subset_taxa(ps_sp_noUnk, taxa_names(ps_sp_noUnk) %in% c(taxa_names(p5_N), taxa_names(p5_P), taxa_names(p5_T)))
p5_NPT.RMP_138_species <- p5_NPT

# Write taxa names to a CSV file
write.csv(taxa_names(p5_NPT.RMP_138_species), "sp138_names.csv")
write.csv(taxa_names(subset_taxa(ps_sp_noUnk, !taxa_names(ps_sp_noUnk) %in% taxa_names(p5_NPT.RMP_138_species))), "spNot138_names.csv")

# Perform Kruskal-Wallis tests
sub_DF <- data.frame(sample_data(p5_NPT)[, c("CRC_general_status")])
OTUdf <- as.data.frame(as(otu_table(p5_NPT), "matrix"))
dataStatus <- cbind.data.frame(sub_DF, OTUdf)

kw_chi <- lapply(dataStatus[2:length(dataStatus)], function(x) {kruskal.test(x ~ CRC_general_status, data = dataStatus)$statistic})
kw_pvl <- lapply(dataStatus[2:length(dataStatus)], function(x) {kruskal.test(x ~ CRC_general_status, data = dataStatus)$p.value})

# Create a data frame with the results of Kruskal-Wallis tests
kw_res <- as.data.frame(cbind(as.data.frame(p.adjust(kw_chi, "none")),
                              as.data.frame(p.adjust(kw_pvl, "none")),
                              as.data.frame(p.adjust(kw_pvl, "BH"))))
setnames(kw_res, c("Chi2", "P", "AdjustedP"))

RMP_138sp_SLV_sp_kw_res <- kw_res
dim(RMP_138sp_SLV_sp_kw_res)

# Write the results to a CSV file
write.csv(RMP_138sp_SLV_sp_kw_res, "RMP_138sp_SLV_sp_kw_res.csv")

# Subset significant results based on adjusted p-values
RMP_SLV_sp_kw_res_sig <- row.names(subset(kw_res, AdjustedP < 0.05))

# Create a data frame with significant results and CRC_general_status values
RMP_SLV_sp_kw_res_sig_df <- dataStatus[, c("CRC_general_status", RMP_SLV_sp_kw_res_sig)]
dim(RMP_SLV_sp_kw_res_sig_df)

# Modify row names and variable names for further analysis
rownames(RMP_138sp_SLV_sp_kw_res) <- gsub("-", ".", rownames(RMP_138sp_SLV_sp_kw_res))
rownames(RMP_138sp_SLV_sp_kw_res) <- gsub(" ", ".", rownames(RMP_138sp_SLV_sp_kw_res))
names(dataStatus) <- gsub("-", ".", names(dataStatus))
names(dataStatus) <- gsub(" ", ".", names(dataStatus))

# Assign variables for further analysis
data <- dataStatus
outcome_vars <- rownames(RMP_138sp_SLV_sp_kw_res)

# Perform Kruskal-Wallis effect size analysis
kruskal_effsize_results <- list()

for (outcome_var in outcome_vars) {
  formula <- as.formula(paste(outcome_var, "~ CRC_general_status"))
  result <- data %>%
    kruskal_effsize(formula = formula)
  result$Outcome <- outcome_var
  kruskal_effsize_results[[outcome_var]] <- result
}

# Combine the results into a single data frame
combined_results <- bind_rows(kruskal_effsize_results)
head(combined_results)

# Merge Kruskal-Wallis results with effect size results
RMP_138sp_SLV_sp_kw_res_effsize <- merge(RMP_138sp_SLV_sp_kw_res, combined_results, by.x = "row.names", by.y = "Outcome", all = TRUE)
write.csv(RMP_138sp_SLV_sp_kw_res_effsize, "RMP_138sp_SLV_sp_kw_res_effsize.csv")

# Perform Dunn's post hoc tests
colname1 <- names(RMP_SLV_sp_kw_res_sig_df)[1]
colname2 <- names(RMP_SLV_sp_kw_res_sig_df)[2:length(RMP_SLV_sp_kw_res_sig_df)]

RMP_SLV_sp_kw_res_sig_df_dunnt <- lapply(colname2, function(x) {
  rstatix::dunn_test(RMP_SLV_sp_kw_res_sig_df, reformulate(colname1, x),
                     p.adjust.method = "fdr")
})

RMP_SLV_sp_kw_res_sig_df_dunnt_r <- Reduce(full_join, RMP_SLV_sp_kw_res_sig_df_dunnt)
RMP_SLV_sp_kw_res_sig_df_dunnt_r$compari <- paste(RMP_SLV_sp_kw_res_sig_df_dunnt_r$group1, RMP_SLV_sp_kw_res_sig_df_dunnt_r$group2, sep = ".")
RMP_SLV_sp_kw_res_sig_df_dunnt_r

write.csv(RMP_SLV_sp_kw_res_sig_df_dunnt_r, "RMP_SLV_sp_kw_res_sig_df_dunnt_r.csv")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QMP @ Species SLV  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Assign variables
ps_sp <- QMP_ASV_LCPM_SLV_species
df_names <- data.frame(t(otu_table(ps_sp)))

# Remove taxa with 'unclassified' in their names
ps_sp_noUnk <- subset_taxa(ps_sp, !taxa_names(ps_sp) %in% taxa_names(ps_sp)[grep('.unclassified', taxa_names(ps_sp))])

# Modify taxa names by removing unnecessary characters
taxa_names(ps_sp_noUnk) <- gsub(".*\\g_", "", taxa_names(ps_sp_noUnk))
tax_table(ps_sp_noUnk)[, 7] <- gsub(".*\\g_", "", tax_table(ps_sp_noUnk)[, 7])

# Load saved data
load("p5_NPT.QMP_138_species")

# Assign variables for further analysis
ps1 <- p5_NPT.QMP_138_species
# Subset samples based on CRC_general_status values
sub_DF <- data.frame(sample_data(ps1)[, c("CRC_general_status")])
OTUdf <- as.data.frame(as(otu_table(ps1), "matrix"))
dataStatus <- cbind.data.frame(sub_DF, OTUdf)

# Perform Kruskal-Wallis tests
kw_chi <- lapply(dataStatus[2:length(dataStatus)], function(x) {kruskal.test(x ~ CRC_general_status, data = dataStatus)$statistic})
kw_pvl <- lapply(dataStatus[2:length(dataStatus)], function(x) {kruskal.test(x ~ CRC_general_status, data = dataStatus)$p.value})

# Create a data frame with the results of Kruskal-Wallis tests
kw_res <- as.data.frame(cbind(as.data.frame(p.adjust(kw_chi, "none")),
                              as.data.frame(p.adjust(kw_pvl, "none")),
                              as.data.frame(p.adjust(kw_pvl, "BH"))))
setnames(kw_res, c("Chi2", "P", "AdjustedP"))

# Assign variables for further analysis
QMP_138sp_SLV_sp_kw_res <- kw_res

# Subset significant results based on adjusted p-values
QMP_SLV_sp_kw_res_sig <- row.names(subset(kw_res, AdjustedP < 0.05))

# Create a data frame with significant results and CRC_general_status values
QMP_SLV_sp_kw_res_sig_df <- dataStatus[, c("CRC_general_status", QMP_SLV_sp_kw_res_sig)]

# Modify row names and variable names for further analysis
rownames(QMP_138sp_SLV_sp_kw_res) <- gsub("-", ".", rownames(QMP_138sp_SLV_sp_kw_res))
rownames(QMP_138sp_SLV_sp_kw_res) <- gsub(" ", ".", rownames(QMP_138sp_SLV_sp_kw_res))
names(dataStatus) <- gsub("-", ".", names(dataStatus))
names(dataStatus) <- gsub(" ", ".", names(dataStatus))

# Assign variables for further analysis
data <- dataStatus
outcome_vars <- rownames(QMP_138sp_SLV_sp_kw_res)

# Perform Kruskal-Wallis effect size analysis
kruskal_effsize_results <- list()

for (outcome_var in outcome_vars) {
  formula <- as.formula(paste(outcome_var, "~ CRC_general_status"))
  result <- data %>%
    kruskal_effsize(formula = formula)
  result$Outcome <- outcome_var
  kruskal_effsize_results[[outcome_var]] <- result
}

# Combine the results into a single data frame
combined_results <- bind_rows(kruskal_effsize_results)

# Merge Kruskal-Wallis results with effect size results
QMP_138sp_SLV_sp_kw_res_effsize <- merge(QMP_138sp_SLV_sp_kw_res, combined_results, by.x = "row.names", by.y = "Outcome", all = TRUE)
write.csv(QMP_138sp_SLV_sp_kw_res_effsize, "QMP_138sp_SLV_sp_kw_res_effsize.csv")

# Perform Dunn's post hoc tests
colname1 <- names(QMP_SLV_sp_kw_res_sig_df)[1]
colname2 <- names(QMP_SLV_sp_kw_res_sig_df)[2:length(QMP_SLV_sp_kw_res_sig_df)]

QMP_SLV_sp_kw_res_sig_df_dunnt <- lapply(colname2, function(x) {
  rstatix::dunn_test(QMP_SLV_sp_kw_res_sig_df, reformulate(colname1, x),
                     p.adjust.method = "fdr")
})

QMP_SLV_sp_kw_res_sig_df_dunnt_r <- Reduce(full_join, QMP_SLV_sp_kw_res_sig_df_dunnt)
QMP_SLV_sp_kw_res_sig_df_dunnt_r$compari <- paste(QMP_SLV_sp_kw_res_sig_df_dunnt_r$group1, QMP_SLV_sp_kw_res_sig_df_dunnt_r$group2, sep = ".")
QMP_SLV_sp_kw_res_sig_df_dunnt_r
write.csv(QMP_SLV_sp_kw_res_sig_df_dunnt_r, "QMP_SLV_sp_kw_res_sig_df_dunnt_r.csv")


################################################
################################################
################################################
#### 
################################################
################################################
################################################