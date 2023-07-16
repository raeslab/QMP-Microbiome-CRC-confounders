# Set the seed for reproducibility
set.seed(531)

library(phyloseq)
library(ggpubr)
library(tibble)
library(dplyr)
library(microViz)
library(vegan)
library(pairwiseAdonis)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


########################################################################
########################################################################

#load files
load("LCPM_all_names_edited_filtered")



# add type variable
sample_data(LCPM_all_names_edited_filtered)$type <- as.factor(sample_data(LCPM_all_names_edited_filtered)$type)
levels(sample_data(LCPM_all_names_edited_filtered)$type) <- c('no_RS', 'no_RS', 'no_RS', 'RS', 'no_RS')

LCPM_all_names_edited_filtered
LCPM_all_names_edited_filtered_species <- tax_glom(LCPM_all_names_edited_filtered, taxrank = "Species")
LCPM_all_names_edited_filtered_species


##### Runella reads per runs and plates
# add number of Runella reads in each sample
sample_data(LCPM_all_names_edited_filtered_species)$Runella <-
  sample_sums(subset_taxa(LCPM_all_names_edited_filtered_species, Genus=="g_Runella"))



dt_fr <- data.frame(sample_data(LCPM_all_names_edited_filtered_species))
RC_plot <- 
  ggboxplot(subset(dt_fr, type%in%c("no_RS", "RS")), x="RUN", add="jitter", y="Runella",  fill="type", color = "type", alpha=0.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(legend.position = "right") #+ scale_y_continuous(trans='log2') + 
RC_plot
ggsave("FigSup_RC_plot.pdf", width = 6, height = 4, RC_plot, useDingbats=FALSE)



##### Positive control consistency among runs and plates

LCPM_all_names_edited_filtered_species_10K <- rarefy_even_depth(LCPM_all_names_edited_filtered_species, sample.size = 10000, rngseed = 531)
LCPM_all_names_edited_filtered_species_10K_Controls <- subset_samples(LCPM_all_names_edited_filtered_species_10K, !type.1=='stool')

to_plot <- LCPM_all_names_edited_filtered_species_10K
metadata_all <- as.tibble(sample_data(to_plot))

BrayDD=phyloseq:: distance(to_plot, "bray")
distances_tax = BrayDD
distances_tax <- as.matrix(distances_tax)
namesPairs <- t(combn(colnames(distances_tax), 2))
distances_tax <- data.frame(namesPairs, dist=distances_tax[namesPairs], stringsAsFactors = F)
distances_tax <- as.tibble(distances_tax)


metadata_all$FastQ_full_ID <- sample_names(to_plot)
metadata_all$Sample_ID <- metadata_all$type.1
distances_tax <- distances_tax %>%
  left_join(metadata_all, by=c("X1"="FastQ_full_ID")) %>%
  dplyr::select(X1,X2,dist,Sample_ID) %>%
  left_join(metadata_all, by=c("X2"="FastQ_full_ID")) %>%
  dplyr::select(X1,X2,dist,Sample_ID.x,Sample_ID.y)

cat_vector <- c()
for(i in 1:nrow(distances_tax)){
  if(distances_tax$Sample_ID.x[i]!=distances_tax$Sample_ID.y[i]){
    cat_vector[i]="Intersample"
  } else {
    cat_vector[i]=paste0(distances_tax$Sample_ID.x[i],"-",distances_tax$Sample_ID.y[i])
  }
}

distances_tax <- distances_tax %>% mutate(clase=cat_vector)
distances_tax_sub = subset(distances_tax, !clase=="Intersample")
distances_tax_sub = distances_tax



distances_tax_sub$RUN <- (data.frame(do.call("rbind", strsplit(as.character(distances_tax_sub$X1), "p", fixed = TRUE)))[,1])

distances_tax_sub$run_run <-
  paste(substr(distances_tax_sub$X1 , 1, 7),
        substr(distances_tax_sub$X2 , 1, 7), sep = '_')


PC_bc_dis <- ggplot(subset(distances_tax_sub, clase == "PC-PC"), aes(x = run_run, y = dist, fill = clase, color = clase)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.7) +
  geom_hline(aes(yintercept = 0.2), linetype = 'dotted', color = 'blue') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none') +
  labs(title = "Bray-Curtis dissimilarity distances between\npositive control (samples rarefied to 10K reads)",
       x = "Run-Run",
       y = "Bray-Curtis distances")
PC_bc_dis
ggsave("FigSup_PC_bc_final_sp.pdf", width = 8, height = 4, PC_bc_dis, useDingbats=FALSE)


LCPM_all_names_edited_filtered_species_noRS <- subset_samples(LCPM_all_names_edited_filtered_species, !type=="RS")

LCPM_all_names_edited_filtered_species_noRS %>%
  tax_transform("compositional", rank = "Genus") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>%
  ord_plot(color = "type.1", shape = "type.1", size = 2) +
  scale_colour_brewer(palette = "Set1")


# Modify taxa table
# Update taxa names



LCPM_all_names_edited_filtered_species_noRS <- subset_taxa(LCPM_all_names_edited_filtered_species_noRS, taxa_sums(LCPM_all_names_edited_filtered_species_noRS)>0)
LCPM_all_names_edited_filtered_species_noRS

tax_table(LCPM_all_names_edited_filtered_species_noRS)[, 7] <- gsub(".*\\g_", "", tax_table(LCPM_all_names_edited_filtered_species_noRS)[, 7])

LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species <- subset_samples(LCPM_all_names_edited_filtered_species_noRS, type.1%in%c("NC", "NCE"))
LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species
taxa_names(LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species) <- tax_table(LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species)[,7]
LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species <- subset_taxa(LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species, taxa_sums(LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species)>0)
LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species
# Create NCE_data and export as CSV
NCE_data <- cbind(sample_data(LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species)[, c("type.1")], otu_table(LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species))
write.csv(NCE_data, "NCE_dataSp.csv")

# Filter names based on publishedTW
#names(NCE_data)[names(NCE_data) %in% publishedTW]

# Create map_QC and calculate distance
LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species

map_QC <- data.frame(sample_data(LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species))
QC_dist.bc <- phyloseq::distance(LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species, method = "bray")

# Perform adonis and assign results
adonis_QC <- adonis(QC_dist.bc ~ RUN, data = map_QC, "data.frame", perm = 9999, na.action = na.omit)
adonis_QC$aov.tab

# Load pairwiseAdonis library and perform pairwise adonis test

pairwise.adonis(QC_dist.bc, map_QC$RUN, p.adjust.m = "bonferroni")

# Update Location_type_colon as factor
sample_data(LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species)$type.1 <- as.factor(sample_data(LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species)$type.1)
levels(sample_data(LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species)$type.1) <- c("NCP", "NCE")

# Assign to psq and apply filters
psq <- LCPM_all_names_edited_filtered_noRS_sp_NCE_NC_sub_species
psq <- tax_filter(psq, min_prevalence = 1 / 10, min_sample_abundance = 1 / 10)

# Perform taxonomic aggregation
psq <- tax_agg(psq, "Species")

# Create complex heatmap and save as PDF
htmp <- psq %>%
  tax_transform("compositional", rank = "Species") %>%
  comp_heatmap(tax_anno = taxAnnotation(
    Prev. = anno_tax_prev(bar_width = 0.3, size =
                            
                            grid::unit(1, "cm"))
  ),
  sample_anno = sampleAnnotation(
    RUN = anno_sample_cat("RUN"), col = list(RUN = col), border = FALSE,
    Type = anno_sample_cat("type.1", col = c("blue", "grey"))
  ))

htmp %>% ComplexHeatmap::draw(
  annotation_legend_list = attr(htmp, "AnnoLegends"), cluster_columns = FALSE, cluster_rows = FALSE)


pdf("QC_NC_NCEn1.pdf", width = 7, height = 4)
htmp %>% ComplexHeatmap::draw(
  annotation_legend_list = attr(htmp, "AnnoLegends"), cluster_columns = FALSE, cluster_rows = FALSE)
dev.off()
