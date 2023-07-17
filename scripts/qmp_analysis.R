###########################################  
###########################################  
# QMP @ ASV (Amplicon sequence variant)
# Modified from our previous script for Quantitative Microbial Profiling (QMP) at the genus level, using RDP taxonomy and 16S rRNA copy number (below link)
https://github.com/raeslab/temporal_variability_microbiome_qmp/blob/master/QMP_optimumevensamplingdepth.R
###########################################  
########################################### 

set.seed(531) 
library(doParallel)
library(Biostrings)
library(RasperGade16S)

# Estimate GCN 16S rRNA
LCPM_runs_589_seq_new_SLV_nr99_v138.1_filtered

# Generate list of DNA sequences
list_of_DNAStrings <- refseq(LCPM_runs_589_seq_new_SLV_nr99_v138.1_filtered)

# Uncomment the line below to write the DNA sequences to a FASTA file
writeXStringSet(list_of_DNAStrings, 'LCPM_runs_589_seq_new_SLV_nr99_v138.1_filtered.fa')

# Uncomment the line below if you have previously generated the GCN predictions and want to load them
LCPM_runs_589_seq_new_SLV_nr99_v138.1_filtered.GCN <- predict_16SGCN_from_sequences(seqs = "LCPM_runs_589_seq_new_SLV_nr99_v138.1_filtered.fa")


# Uncomment the line below to save the GCN predictions
# save(LCPM_runs_589_seq_new_SLV_nr99_v138.1_filtered.GCN, file = "LCPM_runs_589_seq_new_SLV_nr99_v138.1_filtered.GCN")
load("LCPM_runs_589_seq_new_SLV_nr99_v138.1_filtered.GCN")
load("LCPM_runs_589_seq_new_SLV_nr99_v138.1_filtered")

# QMP
Phyl.Object.ASV <- LCPM_runs_589_seq_new_SLV_nr99_v138.1_filtered

# Load the CNV database
df_cnv <- LCPM_runs_589_seq_new_SLV_nr99_v138.1_filtered.GCN$tab

# Set row names for the CNV database
rownames(df_cnv) <- df_cnv$label
head(df_cnv)

# Check for taxa not present in the CNV database
taxa_names(Phyl.Object.ASV)[!(taxa_names(Phyl.Object.ASV) %in% row.names(df_cnv))]

# CNV correction
otu.copy <- data.frame(merge(t(otu_table(Phyl.Object.ASV)), df_cnv, by = "row.names", all.x = TRUE), row.names = 1)
otu.copy <- otu.copy[, !names(otu.copy) %in% c("label", "probs")]
otu.cnv <- otu.copy[, names(otu.copy) != "x"] / otu.copy[, "x"] # divide each taxon by copy number

metadata_QMP <- data.frame(sample_data(Phyl.Object.ASV))
head(metadata_QMP)

write.table(as.data.frame(metadata_QMP[, "Cell_counts"]), "Cell_counts_tempTem.txt")
counts <- read.table("Cell_counts_tempTem.txt", sep = " ", header = TRUE, row.names = 1)
head(counts)

# Rarefying to even sampling depth
rarefy_even_sampling_depth_opt <- function(cnv_corrected_abundance_table, cell_counts_table, minimum_nr_reads) {
  try(if (identical(sort(row.names(cnv_corrected_abundance_table)), sort(row.names(cell_counts_table))) == FALSE) stop("cnv_corrected_abundance_table and cell_counts_table do not have the same sample names. Please check!"))
  
  cnv_corrected_abundance_table <- ceiling(cnv_corrected_abundance_table) # data values are rounded up in order to make use of integer values during the calculations
  
  cell_counts_table <- t(cell_counts_table[order(match(row.names(cnv_corrected_abundance_table), row.names(cell_counts_table))), , drop = FALSE]) # order cell_counts_table as cnv_corrected_abundance_table
  
  sample_sizes <- rowSums(cnv_corrected_abundance_table) # sample size of each sample (total nr of reads)
  sampling_depths <- sample_sizes / cell_counts_table # sampling depth of each sample (total nr of reads divided by the cell count)
  
  minimum_sampling_depth <- min(sampling_depths) # minimum of all sampling depths
  
  rarefy_to <- cell_counts_table * minimum_sampling_depth # nr of reads to rarefy in each sample in order to get to an even sampling depth over all samples
  if (all(rarefy_to > minimum_nr_reads)) { # if all samples above min nr reads --> no problem
    samples_to_exclude <- c()
    minimum_sampling_depth_opt <- minimum_sampling_depth
  } else { # try next sampling depth (and exclude x)
    lost_max <- length(which(rarefy_to < minimum_nr_reads)) # n samples that are below nr reads
    lost_due_to_exclusion <- 1:(length(sampling_depths)) # vector with all n
    lost_after_rarefaction <- c()
    
    ## check how many samples are lost after rarefaction (because they did not reach the minimum nr of reads specified by the threshold, for each sample with a low sampling depth that is removed)
    for (i in 1:(length(sampling_depths))) {
      minimum_sampling_depth <- sort(sampling_depths, decreasing = TRUE)[length(sampling_depths) - i] # take sampling depth of next sample
      rarefy_to <- cell_counts_table * minimum_sampling_depth # nr of reads to rarefy in each sample in order to get to an even sampling depth over all samples
      lost_after_rarefaction <- c(lost_after_rarefaction, length(which(rarefy_to < minimum_nr_reads))) # for each sampling depth --> n samples that are lost
    }
    
    ## define the optimal sampling depth (depth at which the minimum of samples are lost in total)
    lost_in_total <- lost_after_rarefaction + lost_due_to_exclusion
    nr_samples_to_exclude <- min(which(lost_in_total == min(lost_in_total))) # CAUSES PROBLEMS IF SEVERAL WITH MIN
    
    plot(lost_in_total, xlim = c(0, lost_max), ylim = c(0, lost_max), pch = 16,
         ylab = "Total nr of samples lost",
         xlab = "Nr. of excluded samples",
         main = paste0("Lost using optimal sampling depth: ", min(lost_in_total), " (instead of ", lost_max, ")"))
    
    mtext(paste(nr_samples_to_exclude, "due to exclusion and ", lost_after_rarefaction[nr_samples_to_exclude], " because they didn't make the rarefaction threshold."))
    abline(h = min(lost_in_total), v = nr_samples_to_exclude)
    minimum_sampling_depth_opt <- sort(sampling_depths, decreasing = TRUE)[length(sampling_depths) - nr_samples_to_exclude] # minimum of all sampling depths at which the least samples are lost to reach the minimum nr of reads.
    # define the samples that need to be excluded for this optimal sampling depth
    samples_to_exclude <- which(sampling_depths < minimum_sampling_depth_opt)
  }
  
  # make QMP matrix with optimal sampling depth
  rarefy_to_opt <- round(cell_counts_table * minimum_sampling_depth_opt) # nr of reads to rarefy in each sample in order to get to the optimum even sampling depth over all samples
  cnv_corrected_abundance_table_phyloseq <- otu_table(cnv_corrected_abundance_table, taxa_are_rows = FALSE) # convert to phyloseq otutable
  out <- NULL
  samples_included <- c()
  
  for (i in 1:nrow(cnv_corrected_abundance_table_phyloseq)) {
    if (!(i %in% samples_to_exclude)) { # skip the samples that need to be excluded
      if (rarefy_to_opt[i] > minimum_nr_reads) { # only include the samples that pass the threshold
        print(paste("sample", row.names(cnv_corrected_abundance_table_phyloseq)[i], "rarefied to", rarefy_to_opt[i], "reads."))
        x <- rarefy_even_depth(cnv_corrected_abundance_table_phyloseq[i, ], sample.size = rarefy_to_opt[i], rngseed = 711, replace = FALSE, trimOTUs = FALSE, verbose = FALSE)
        out <- rbind(out, x)
        samples_included <- c(samples_included, i) # make a final list of samples that get included
      }
    }
  }
  
  print("_____________________________________________________________")
  print("Optimal sampling depth:")
  print(minimum_sampling_depth_opt)
  print("_____________________________________________________________")
  print(paste("The following", length(samples_to_exclude), "samples were excluded due to exclusion:"))
  print(row.names(cnv_corrected_abundance_table)[samples_to_exclude])
  #print(paste("The following", lost_after_rarefaction[nr_samples_to_exclude], "samples were excluded because of not making the rarefaction threshold:"))
  print(row.names(cnv_corrected_abundance_table)[!(row.names(cnv_corrected_abundance_table)) %in% c(row.names(cnv_corrected_abundance_table)[samples_to_exclude], row.names(cnv_corrected_abundance_table)[samples_included])])
  
  rarefied_matrix <- as.matrix(out)
  normalised_rarefied_matrix <- rarefied_matrix / rowSums(rarefied_matrix)
  QMP <- normalised_rarefied_matrix * cell_counts_table[1, samples_included]
  return(QMP)
}

cnv_corrected_abundance_table <- t(otu.cnv)
cell_counts_table <- counts[, 1, drop = FALSE]
minimum_nr_reads <- 500 / mean(otu.copy[, "x"]) 

dim(cnv_corrected_abundance_table)
dim(cell_counts_table)
rownames(cnv_corrected_abundance_table)

QMP <- rarefy_even_sampling_depth_opt(cnv_corrected_abundance_table, cell_counts_table, minimum_nr_reads) 

#[1] "Optimal sampling depth:"
#[1] 1.496569e-08


################################
QMP_ASV_LCPM_SLV <- phyloseq(otu_table(QMP, taxa_are_rows = FALSE), 
                             sample_data(Phyl.Object.ASV), 
                             tax_table(Phyl.Object.ASV), 
                             refseq(Phyl.Object.ASV))
QMP_ASV_LCPM_SLV 
# save(QMP_ASV_LCPM_SLV, file = "QMP_ASV_LCPM_SLV")
