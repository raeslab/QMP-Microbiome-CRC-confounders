# Set seed for reproducibility
set.seed(531)

# Load necessary libraries
library(phyloseq)
library(compositions)
library(ape)
library(DirichletMultinomial)
library(ggplot2)


# Get the directory path of the active document in RStudio and set it as the working directory
doc_path <- dirname(rstudioapi::getActiveDocumentContext()$path)

### ETs ###

# Load the required datasets

# FGFP data set background
load("fgfp_10Kg")
# Cohort to be enterotyped
load("lcpm_10Kg")

# Check for differences in taxa names between datasets
diff_lcpm_fgfp <- setdiff(taxa_names(lcpm_10Kg), taxa_names(fgfp_10Kg))
diff_fgfp_lcpm <- setdiff(taxa_names(fgfp_10Kg), taxa_names(lcpm_10Kg))


# Add cohort information to the sample data
sample_data(lcpm_10Kg)$Cohort <- c("LCPM")
sample_data(fgfp_10Kg)$Cohort <- c("FGFP")

# Merge the two datasets into one
lcpm_fgfp_10Kg <- merge_phyloseq(fgfp_10Kg, lcpm_10Kg)

# Perform PCoA analysis using Bray distance metric
BC.ord <- ordinate(lcpm_fgfp_10Kg, "PCoA", "bray")

# Visualize the PCoA results colored by cohort
plot_ordination(lcpm_fgfp_10Kg, BC.ord, type = "samples", color = "Cohort")

# Get the genus-level OTU table
GENUS2 <- t(lcpm_fgfp_10Kg@otu_table) ### change with 't'
GENUS2[2, 3] # columns are samples??

# Visual check to see if the merging was done correctly: if the study splits apart - bad sign
# Print dimensions and sum ranges of the OTU table
dim(GENUS2)
range(colSums(GENUS2))
range(rowSums(GENUS2))

#########################################################################################################
# DMM ENTEROTYPING
#########################################################################################################

# Set the input OTU table for enterotyping
otutable_file <- GENUS2 # columns are samples
otutable_file[2, 3]

# Add pseudocount to the OTU table
otu_table <- otutable_file + 0.000000001

# Create an empty dataframe to store enterotyping results
enterotype_df <- data.frame(matrix(nrow = ncol(otu_table), ncol = 0), row.names = colnames(otu_table))

# Compute dirichlet multinomial for different numbers of components
all_dmns <- 6 # how many DMMs should be checked?
dmn_list <- numeric(all_dmns)
for (i in 1:all_dmns) {
  print(i)
  assign(paste0("dmn_", i), dmn(as.matrix(t(otutable_file)), i, verbose = F))
}
dmn_list <- list(dmn_1, dmn_2, dmn_3, dmn_4, dmn_5, dmn_6)

# Save the DMM objects to a file
save.image(file = "fgfp_lcpm_10K_et_dmms.rdata")

# Find the DMM with the minimum number of Dirichlet components based on laplace
lplc <- sapply(dmn_list, laplace)
BIC <- sapply(dmn_list, BIC)
AIC <- sapply(dmn_list, AIC)
dmn_min_lplc <- dmn_list[[which.min(lplc)]]
dmn_min_BIC <- dmn_list[[which.min(BIC)]]
dmn_min_AIC <- dmn_list[[which.min(AIC)]]
print("DMM with minimum number of Dirichlet components based on Laplace:", dmn_min_lplc)
print("DMM with minimum number of Dirichlet components based on BIC:", dmn_min_BIC)
print("DMM with minimum number of Dirichlet components based on AIC:", dmn_min_AIC)

# Plot model fit for different numbers of Dirichlet components
pdf("fgfp_lcpm_10K_min_dirichlet_components.pdf", onefile = TRUE)
plot(lplc, type = "b", xlab = "Number of Dirichlet Components", ylab = "Model Fit, Laplace") # 4
plot(BIC, type = "b", xlab = "Number of Dirichlet Components", ylab = "Model Fit, BIC") # 2
plot(AIC, type = "b", xlab = "Number of Dirichlet Components", ylab = "Model Fit, AIC") # 4 (3 before)
dev.off()

# Load necessary library

# Extract the DMM results for different numbers of components
lplc <- sapply(dmn_list, laplace)
dmn1 <- dmn_list[[1]]
Dirichlet_multinomial_1 <- mixture(dmn1, assign = TRUE)

dmn2 <- dmn_list[[2]]
Dirichlet_multinomial_2 <- mixture(dmn2, assign = TRUE)

dmn3 <- dmn_list[[3]]
Dirichlet_multinomial_3 <- mixture(dmn3, assign = TRUE)

dmn4 <- dmn_list[[4]]
Dirichlet_multinomial_4 <- mixture(dmn4, assign = TRUE)

dmn5 <- dmn_list[[5]]
Dirichlet_multinomial_5 <- mixture(dmn5, assign = TRUE)

dmn6 <- dmn_list[[6]]
Dirichlet_multinomial_6 <- mixture(dmn6, assign = TRUE)

# Combine the DMM results into one dataframe
Dirichlet_multinomial_all <- cbind(Dirichlet_multinomial_1, Dirichlet_multinomial_2, Dirichlet_multinomial_3,
                                   Dirichlet_multinomial_4, Dirichlet_multinomial_5, Dirichlet_multinomial_6)
colnames(Dirichlet_multinomial_all) <- c("k1", "k2", "k3", "k4", "k5", "k6")


# Perform PCoA analysis using Bray distance metric
BC.ord <- ordinate(lcpm_fgfp_10Kg, "PCoA", "bray")
plot_ordination(lcpm_fgfp_10Kg, BC.ord, type = "samples", color = "K4")

# Add enterotyping results to the sample data
sample_data(lcpm_fgfp_10Kg)$K4 <- as.factor(Dirichlet_multinomial_all$k4)

### Object with enterotypes ####
#et_output = read.table("Enterotyping_Calpro_PSCFGFP_1to6ET_new.tsv",sep="\t",header=T, row.names=1)
et_output = as.data.frame(Dirichlet_multinomial_all)
table(et_output$k4)
n=row.names(et_output)
et_output=data.frame(apply(et_output,2,function(x) as.factor(x)))
row.names(et_output)=n
summary(et_output)

genus_table=data.frame(t(GENUS2)) 
taxa_list=c("g_Bacteroides","g_Faecalibacterium","g_Prevotella","g_Alistipes","g_Ruminococcus","g_Methanobrevibacter") #TAXA to plot at end ("Bacteroides","Faecalibacterium","Prevotella" are mandatory)
row.names(GENUS2)

et_column="k4" #name of your ET clusters column

#change Genus matrix to subselect relevant taxa and merge all others in "Other"
genus_table=data.frame(merge(data.frame("Other"=rowSums(genus_table[,-which(colnames(genus_table)%in%taxa_list)])),genus_table[,which(colnames(genus_table)%in%taxa_list)],by="row.names"),row.names=1)

#merge Genus matrix with enterotype definition
genus_table=data.frame(merge(genus_table,et_output[,et_column,drop=F],by="row.names"),row.names=1)

#decision to translate cluster number into ET name (B1/B2/R/P) only works if optimal is 4

genera_summary=aggregate(. ~ genus_table[,et_column], genus_table, median)
colnames(genera_summary)=gsub("genus_table[, et_column]", "clusters", colnames(genera_summary),fixed=TRUE) 

tmp=genera_summary
P=tmp[which(tmp$g_Prevotella==max(tmp$g_Prevotella,na.rm=TRUE)),"clusters"]
paste("Prev","==",P)
genera_summary[order(genera_summary[,"g_Prevotella",drop=FALSE]),c("clusters","g_Prevotella")]
tmp[which(tmp$clusters==P),]=NA

B2=tmp[which(tmp$g_Bacteroides==max(tmp$g_Bacteroides,na.rm=TRUE)),"clusters"]
paste("Bact2","==",B2)
genera_summary[order(genera_summary[,"g_Bacteroides",drop=FALSE]),c("clusters","g_Bacteroides")]
tmp[which(tmp$clusters==B2),]=NA

B1=tmp[which(tmp$g_Bacteroides==max(tmp$g_Bacteroides,na.rm=TRUE)),"clusters"]
paste("Bact1","==",B1)
genera_summary[order(genera_summary[,"g_Bacteroides",drop=FALSE]),c("clusters","g_Bacteroides","g_Faecalibacterium")]
tmp[which(tmp$clusters==B1),]=NA

R=tmp[!is.na(tmp$clusters),"clusters"]

ET_tab=genus_table[, et_column,drop=FALSE]
colnames(ET_tab)="ET_clusters"
ET_names=ET_tab[,1]
ET_names=gsub(paste0("\\b",P,"\\b"), "Prev", ET_names)
ET_names=gsub(paste0("\\b",B1,"\\b"), "Bact1", ET_names)
ET_names=gsub(paste0("\\b",B2,"\\b"), "Bact2", ET_names)
ET_names=gsub(paste0("\\b",R,"\\b"), "Rum", ET_names)
ET_tab=cbind(ET_tab,ET_names) #final list
head(ET_tab)


