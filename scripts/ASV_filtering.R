# Set the seed for reproducibility
set.seed(531)

library(phyloseq)
library(ggpubr)
library(tibble)
library(dplyr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# load dada2 files

########################################################################
########################################################################

#load files

load("seqtab_all_no_chimeras_tax_SLV138.1")
load("seqtab_all_no_chimeras")
load("LCPM_all_ids")

# make phyloseq obj
LCPM_all <- merge_phyloseq(otu_table(seqtab_all_no_chimeras),
                           tax_table(seqtab_all_no_chimeras_tax_SLV138.1), 
                           sample_data(LCPM_all_ids))
LCPM_all <- subset_taxa(LCPM_all, taxa_sums(LCPM_all)>0)
LCPM_all


# renaming taxa names
ASV_phyObj <- LCPM_all

mat=as.data.frame(tax_table(ASV_phyObj))
mat[is.na(mat)] <- "unclassified"
rank_names(ASV_phyObj)

mat=as.data.frame(tax_table(ASV_phyObj))
mat[is.na(mat)] <- "unclassified"
rank_names(ASV_phyObj)

##
mat$Genus = paste("g", mat$Genus, sep = "_")
mat$Family = paste("uc_f", mat$Family, sep = "_")
mat$Order = paste("uc_o", mat$Order, sep = "_")
mat$Class = paste("uc_c", mat$Class, sep = "_")
mat$Phylum = paste("uc_p", mat$Phylum, sep = "_")
mat[mat$Genus == "g_unclassified","Genus"] <- mat[mat$Genus == "g_unclassified","Family"]#
mat[mat$Genus == "uc_f_unclassified","Genus"] <- mat[mat$Genus == "uc_f_unclassified","Order"]#
mat[mat$Genus == "uc_o_unclassified","Genus"] <- mat[mat$Genus == "uc_o_unclassified","Class"]#
mat[mat$Genus == "uc_c_unclassified","Genus"] <- mat[mat$Genus == "uc_c_unclassified","Phylum"]#

# Combine Genus and Species with a dot separator
mat$Species <- paste(mat$Genus, mat$Species, sep = ".")
head(mat$Species)


# phyloseq obj with new taxa names

LCPM_all_names_edited <- phyloseq(otu_table(ASV_phyObj, taxa_are_rows=T),
                                  sample_data(ASV_phyObj),
                                  tax_table(as.matrix(mat)))
LCPM_all_names_edited

#####Remove non bacteria archaea taxa and homo sapiens reads
LCPM_all_names_edited_filtered <- subset_taxa(LCPM_all_names_edited, 
                                              !Order%in%c('uc_o_Chloroplast') &
                                                !Family%in%c('uc_f_Mitochondria') &
                                                !Genus%in%c("uc_p_unclassified") & 
                                                !Genus%in%c("uc_f_Mitochondria"))

save(LCPM_all_names_edited_filtered, file="LCPM_all_names_edited_filtered")


