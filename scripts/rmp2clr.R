##################  from RMP to CLR ############################################
# Load required packages
library(ALDEx2)
library(CoDaSeq)
library(phyloseq)
library(ggplot2)

set.seed(531) 

# Set working directory to the location of the current script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######## CLR ########

# Load the microbiome data object
load("RMP_ASV_LCPM_SLV_species")

# Check and adjust taxonomic names if necessary
taxa_names(RMP_ASV_LCPM_SLV_species) <- gsub("-", ".", taxa_names(RMP_ASV_LCPM_SLV_species))

# Create a working copy of the data object
Obj <- RMP_ASV_LCPM_SLV_species

# Set filtering parameters
min.reads <- 10000
min.prop <- 0.001
cutoff <- 0

# Convert species table and apply filtering
matrix <- data.frame(otu_table(Obj))
matrix <- t(matrix)
matrix.f <- codaSeq.filter(matrix, min.reads = min.reads, min.occurrence = cutoff, min.prop = min.prop, samples.by.row = FALSE)

# Function to estimate zeros
estimate0.min <- function(matrix.test) { 
  print(paste("You have", ncol(matrix.test), "samples. If this is not correct, transpose matrix!"))
  matrix.test.p <- t(t(matrix.test) / rowSums(t(matrix.test)))
  samplesums <- colSums(matrix.test)
  
  matrix.f.n0 <- matrix.test
  for (i in 1:nrow(matrix.f.n0)) {
    min <- min(matrix.test.p[i, ][matrix.test.p[i, ] > 0])
    for (j in 1:ncol(matrix.f.n0)) {
      if (matrix.f.n0[i, j] == 0)
        matrix.f.n0[i, j] <- min * samplesums[j]
    }
  }
  return(matrix.f.n0)
}

# Apply zero estimation and perform CLR transformation
matrix.f.n0 <- estimate0.min(matrix.f)
matrix.f.n0.clr <- codaSeq.clr(matrix.f.n0, samples.by.row = FALSE)
hist(rowSums(matrix.f.n0.clr))

######### Create CLR-transformed phyloseq object #########

CLR_ASV_LCPM_SLV_species <- phyloseq(
  sample_data(Obj),
  otu_table(matrix.f.n0.clr, taxa_are_rows = TRUE),
  tax_table(Obj)
)


# Assign Diagnosis metadata column
sample_data(CLR_ASV_LCPM_SLV_species)$Diagnosis <- sample_data(CLR_ASV_LCPM_SLV_species)$CRC_general_status

# Perform PCoA ordination
CLRLCPM_SLV_species.ord <- ordinate(CLR_ASV_LCPM_SLV_species, "PCoA", "euclidean")

# Plot PCoA with Diagnosis as color
EU_plot_CLR_LCPM_SLV_species <- plot_ordination(
  CLR_ASV_LCPM_SLV_species, CLRLCPM_SLV_species.ord,
  type = "samples", color = "Diagnosis"
) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
  theme(panel.background = element_blank(), legend.position = c(0.85, 0.85))

# Display the plot
EU_plot_CLR_LCPM_SLV_species
