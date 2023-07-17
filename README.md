# Microbiome Analysis Pipeline (LCPM)

This repository contains a comprehensive pipeline for analyzing microbiome data from the LCMP cohort. The pipeline covers a range of analyses including data preprocessing, quality control, taxonomic profiling, statistical analysis, and visualization. Below is a detailed guide on how to use the pipeline for each step of the analysis.

## 1. Downloading fastq files and grouping them by sequencing batches

To begin the analysis, the first step is to download the fastq files from the repository. Make sure you have the necessary permissions and access to the metadata. Once you have the fastq files, group them by sequencing batches for easy management and analysis. Use the following steps:

1. Create a directory structure to organize your data, such as `data/raw_data/` for storing the raw fastq files.
2. Place the fastq files in the corresponding sequencing batch directories within the `data/raw_data/` directory.
3. Make sure to appropriately name the directories to reflect the sequencing batches and include relevant metadata if available.

## 2. Initial analysis using dada2 to obtain ASV (controlling for sequencing batch)

After the fastq files have been grouped into sequencing batches, the next step is to perform initial analysis using the dada2 package to obtain Amplicon Sequence Variants (ASVs) for each sequencing batch. Follow these steps:

1. Install the required dependencies, including R and the dada2 package.
2. Create a script, such as `scripts/dada2_analysis.R`, and load the necessary libraries.
3. Execute the **dada2_analysis.R** script.

## 3. Filtering unclassified and non-bacterial ASVs

Once you have obtained the ASVs for each sequencing batch, it is important to filter out unclassified and non-bacterial ASVs to focus on the microbial taxa of interest. Follow these steps:

1. Create a script, such as `scripts/ASV_filtering.R`, and load the necessary libraries.
2. Read the concatenated ASV file from the `data/processed_data/` directory.
3. Implement filtering criteria to remove unclassified and non-bacterial ASVs based on taxonomic annotations.
4. Execute the **ASV_filtering.R** script.
   

## 4. Exploratory and quality control analysis

After filtering the ASVs, it is essential to perform exploratory and quality control analyses to gain insights into the dataset. Follow these steps:

1. Create a script, such as `scripts/exploratory_analysis.R`, and load the necessary libraries.
2. Read the filtered ASV file from the `data/processed_data/` directory.
3. Execute the **exploratory_analysis.R** script.
   
## 5. Quantitative Microbiome Profiling (QMP) at ASV level

Quantitative Microbiome Profiling (QMP) at the ASV level enables the quantification of microbial abundance. Follow these steps to perform QMP:

1. Create a script, such as `scripts/qmp_analysis.R`, and load the necessary libraries.
2. Read the filtered ASV file from the `data/processed_data/` directory.
3. Execute the **qmp_analysis.R** script.


## 6. Microbiota covariates identification

Identification of microbiota covariates helps in understanding the factors influencing microbial community composition. Follow these steps to identify microbiota covariates:

1. Create a script, such as `scripts/covariate_identification.R`, and load the necessary libraries.
2. Read the filtered ASV files for each sequencing batch from the `data/processed_data/` directory.
3. Prepare the necessary metadata, such as sample characteristics, clinical variables, etc.
4. Execute the **covariate_identification.R** script.

## 7. Differential abundance analysis

To identify taxa showing differential abundance between different conditions or groups, follow these steps:

1. Create a script, such as `scripts/differential_abundance.R`, and load the necessary libraries.
2. Read the filtered ASV files for each sequencing batch from the `data/processed_data/` directory.
3. Execute the **differential_abundance.R** script.
   
## 8. Taxa abundance associations

Investigating associations between taxa abundance and other variables can provide valuable insights. Follow these steps to analyze taxa abundance associations:

1. Create a script, such as `scripts/abundance_associations.R`, and load the necessary libraries.
2. Read the filtered ASV files for each sequencing batch from the `data/processed_data/` directory.
3. Prepare the necessary metadata, including variables of interest for association analysis.
4. Execute the **abundance_associations.R** script.
  
## 9. Linear model analysis

Linear model analysis helps in exploring relationships between multiple covariates and microbial abundance. Follow these steps to perform linear models analysis:

1. Create a script, such as `scripts/linear_models.R`, and load the necessary libraries.
2. Read the filtered ASV files for each sequencing batch from the `data/processed_data/Filtered_ASVs_per_batch/` directory.
3. Prepare the necessary metadata, including multiple covariates of interest.
4. Execute the **linear_models.R** script.

## 10. Enterotyping

Enterotyping is a method to categorize individuals based on their gut microbiota composition. Follow these steps to perform enterotyping analysis:

1. Create a script, such as `scripts/enterotyping.R`, and load the necessary libraries.
2. Read the filtered ASV files from the `data/processed_data/` directory.
3. Execute the **enterotyping.R** script.


