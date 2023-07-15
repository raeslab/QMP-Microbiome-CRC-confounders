# Set the seed for reproducibility
set.seed(531)

# For optimization and improvement of code
library(phyloseq)
library(vegan)
library(data.table)
library(ggplot2)


# Set the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary data

load("QMP_ASV_LCPM_SLV_species")

# Load metadata
LCMP_metadata_589 <- read.csv("LCMP_metadata_589.csv", row.names = 1)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          RDA      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     QMP n=589     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Set the metadata for QMP_ASV_LCPM_SLV_species
sample_data(QMP_ASV_LCPM_SLV_species) <- sample_data(LCMP_metadata_589)

to_test <- QMP_ASV_LCPM_SLV_species
fDSGENUS <- as.data.frame(to_test@otu_table@.Data)
m <- data.frame(sample_data(to_test))
dim(m)

m <- m[, !(names(m) %in% c("Lab_ID", "Enterotype", "Cell_counts"))]
dim(m)  # 94 variables
dim(fDSGENUS)

# Perform analysis using capscale and anova
all <- c()
for (i in 1:ncol(m)) {
  capsc <- capscale(fDSGENUS ~ m[, i], distance = "bray", na.action = na.omit)
  an <- anova.cca(capsc)
  pval <- an["Pr(>F)"][[1]][[1]]
  Fa <- an["F"][[1]][[1]]
  r2 <- RsquareAdj(capsc)[[1]]
  adjr2 <- RsquareAdj(capsc)[[2]]
  all <- rbind(all, cbind(Fa, r2, adjr2, pval))
}

colnames(all) <- c("F", "r2", "adjr2", "p-value")
row.names(all) <- colnames(m)

qval <- p.adjust(all[,"p-value"], method = "BH")
all <- data.frame(all)
allfDSall <- cbind(all, qval)

# Write the results to a CSV file
write.csv(allfDSall, "n589_RDA_94var_QMP_spececies_LCPM_SLV.csv")

# Subset the significant results
allfDS <- allfDSall
allfDS$covariate <- row.names(allfDS)
allfDS <- subset(allfDS, qval < 0.1)

#check how many variables are significantly associated to community distance
sign_cap <- row.names(allfDS)
MET=m[,sign_cap] # subset metadata to the selection
MET=na.exclude(MET) #select metadata without NAs
#prepare genus table
GEN=fDSGENUS[which(row.names(fDSGENUS)%in%row.names(MET)),] #subset Genus matrix to samples with available metadata
distmat=vegdist(GEN,method="bray")
#detach(MET)
attach(MET)
mod0=capscale(distmat ~ 1) #H0: unconstrained ordination
mod1=capscale(distmat ~ ., data=MET) #H1: full constrained ordination, all metadata in MET
step.res<-ordiR2step(mod0, scope=formula(mod1), data=MET ,direction="forward", Pin = 0.05, R2scope = TRUE, pstep = 100, perm.max = 200, trace = F) #forward stepwise dbRDA
# Perform stepwise dbRDA
step.res$call 
step.res$anova

# Save the stepwise dbRDA result as an R object
outdbRDA <- "n589_RDA_94var_QMP_spececies_LCPM_SLV.RData"
save(step.res, file = outdbRDA)

data <- read.csv("n589_RDA_94var_QMP_spececies_LCPM_SLV.csv", row.names = 1)
allfDSall <- data
ordiR2stepDF <- as.data.frame(step.res$anova)
rownames(ordiR2stepDF)
rownames(ordiR2stepDF) <- gsub(".*\\+ ", "", rownames(ordiR2stepDF))
ordiR2stepDF <- ordiR2stepDF[!rownames(ordiR2stepDF) %in% c("<All variables>"), ]
head(ordiR2stepDF)

allfDSall_ordiR2stepDF <- allfDSall[rownames(allfDSall) %in% rownames(ordiR2stepDF), ]
head(allfDSall_ordiR2stepDF)
cbind(allfDSall_ordiR2stepDF, ordiR2stepDF)
for_RDA_plot <- merge(allfDSall_ordiR2stepDF, ordiR2stepDF, by = 'row.names', all = TRUE)
rownames(for_RDA_plot) <- for_RDA_plot$Row.names
for_RDA_plot <- for_RDA_plot[, c("Row.names", "adjr2", "R2.adj")]

setnames(for_RDA_plot, c("Covariate", "individual", "cumulative"))

df.long <- melt(for_RDA_plot)
df.long$R2 <- df.long$value
orderNew <- for_RDA_plot$Covariate[order(for_RDA_plot$cumulative)]
setnames(df.long, c("variable", "R2"), c("Effect size", "Effect size (adj R2)"))

RDA_R2all <- ggplot(df.long, aes(Covariate, `Effect size (adj R2)`, fill = `Effect size`)) +
  geom_bar(stat = "identity", position = "dodge") + scale_fill_brewer() +
  scale_x_discrete(limits = rev(orderNew)) + theme(legend.position = c(0.8, 0.85)) + coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + theme(panel.background = element_blank())
RDA_R2all

# Plot PCoA with colors representing diagnosis
to_plotPCoA <- QMP_ASV_LCPM_SLV_species
sample_data(to_plotPCoA)$Diagnosis <- sample_data(to_plotPCoA)$CRC_general_status
levels(sample_data(to_plotPCoA)$Diagnosis) <- c("CTL", "ADE", "CRC")
to_plotPCoA.bc <- ordinate(to_plotPCoA, "PCoA", "bray")

BC_plot_QMP_ASV_LCPM_SLV_species <- plot_ordination(to_plotPCoA, to_plotPCoA.bc, type = "samples", color = "Diagnosis") +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
  theme(panel.background = element_blank()) + theme(legend.position = c(0.85, 0.85))
BC_plot_QMP_ASV_LCPM_SLV_species

