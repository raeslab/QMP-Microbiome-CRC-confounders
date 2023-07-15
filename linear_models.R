set.seed(531) 
# Load necessary packages
library(stats)
library(data.table)
library(phyloseq)
library(dplyr)
library(stringr)
# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dirname(rstudioapi::getActiveDocumentContext()$path)

# Load QMP and RMP data files
QMP_SLV_sp_kw_res <- read.csv("QMP_SLV_sp_kw_res.csv")
RMP_SLV_sp_kw_res <- read.csv("RMP_SLV_sp_kw_res.csv")

# Filter significant results
QMP_SLV_sp_kw_res_sig05 <- subset(QMP_SLV_sp_kw_res, AdjustedP < 0.05)
RMP_SLV_sp_kw_res_sig05 <- subset(RMP_SLV_sp_kw_res, AdjustedP < 0.05)

# Load QMP and RMP data files for final analysis
QMP_SLV_sp_kw_res <- read.csv("QMP_138sp_SLV_sp_kw_res.csv")
RMP_SLV_sp_kw_res <- read.csv("RMP_138sp_SLV_sp_kw_res.csv")

# Filter significant results for final analysis
QMP_SLV_sp_kw_res_sig05 <- subset(QMP_SLV_sp_kw_res, AdjustedP < 0.05)
RMP_SLV_sp_kw_res_sig05 <- subset(RMP_SLV_sp_kw_res, AdjustedP < 0.05)

# Load taxonomy data
load("p5_NPT.QMP_138_species")
load("p5_NPT.RMP_138_species")

# Subset RMP data for glm analysis
RMP_for_glm <- subset_taxa(p5_NPT.RMP_138_species, taxa_names(p5_NPT.RMP_138_species) %in% RMP_SLV_sp_kw_res_sig05$X)

# Subset QMP data for glm analysis
QMP_for_glm <- subset_taxa(p5_NPT.QMP_138_species, taxa_names(p5_NPT.QMP_138_species) %in% QMP_SLV_sp_kw_res_sig05$QMP)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Perform glm analysis for QMP data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj_phylo <- QMP_for_glm
species_QMP <- data.frame(otu_table(obj_phylo))
mata_crc <- data.frame(sample_data(obj_phylo))
mata_crc <- mata_crc[, c("Calprotectin", "CRC_general_status")]
setnames(mata_crc, c("Calprotectin", "CRC_general_status"), c("calpro", "CRCtype"))
mata_crc$CRCtype <- as.factor(mata_crc$CRCtype)
levels(mata_crc$CRCtype) <- c(1, 2, 3)
mata_crc$CRCtype <- as.numeric(as.character(mata_crc$CRCtype))
COHORTDATA <- mata_crc[, c("calpro", "CRCtype")]

# Check data types and dimensions
is.numeric(COHORTDATA$calpro)
is.numeric(COHORTDATA$CRCtype)

# Perform glm analysis for each feature in QMP data
subsamps <- row.names(species_QMP)
datasets <- list(species_QMP)
named <- "species_QMP"
lmstatsall_QMP <- data.frame()



# Perform glm analysis for QMP data including additional covariates (Moisture, BMI)
species_QMP <- data.frame(otu_table(obj_phylo))
mata_crc <- data.frame(sample_data(obj_phylo))
mata_crc <- mata_crc[, c("CRC_general_status", "Calprotectin", "Moisture", "BMI")]
setnames(mata_crc, c("CRC_general_status", "Calprotectin", "Moisture", "BMI"), c("CRCtype", "Calprotectin", "Moisture", "BMI"))
mata_crc$CRCtype <- as.factor(mata_crc$CRCtype)
levels(mata_crc$CRCtype) <- c(1, 2, 3)
mata_crc$CRCtype <- as.numeric(as.character(mata_crc$CRCtype))
COHORTDATA <- mata_crc

# Check data types and dimensions
is.numeric(COHORTDATA$Calprotectin)
is.numeric(COHORTDATA$Moisture)
is.numeric(COHORTDATA$BMI)
is.numeric(COHORTDATA$CRCtype)


# Perform glm analysis for each feature in QMP data including additional covariates

subsamps <- row.names(species_QMP)
datasets <- list(species_QMP)
named <- "species_QMP"
lmstatsall_CMB_QMP <- data.frame()

for (i in 1:length(datasets)) {
  MAT2 <- datasets[[i]]
  MATBIN <- MAT2[subsamps, ]
  MATBIN[MATBIN > 0] <- 1
  sub <- names(which((colSums(MATBIN) / nrow(MATBIN)) > 0.000000050))
  sub <- setdiff(sub, sub[grep("unclassified", sub)])
  MAT2 <- MAT2[, which(colnames(MAT2) %in% sub)]
  DATATAB <- MAT2
  MET <- COHORTDATA[subsamps, ]
  lmstats <- data.frame("dataset" = character(), "feature" = character(),
                        "calprotectin_estimate" = numeric(), "calprotectin_std.Error" = numeric(), "calprotectin_t.value" = numeric(), "calprotectin_P" = numeric(),
                        "Moisture_estimate" = numeric(), "Moisture_std.Error" = numeric(), "Moisture_t.value" = numeric(), "Moisture_P" = numeric(),
                        "BMI_estimate" = numeric(), "BMI_std.Error" = numeric(), "BMI_t.value" = numeric(), "BMI_P" = numeric(),
                        "anova.F" = numeric(), "anova_P" = numeric(), "N" = numeric())
  lmstats$dataset <- as.character()
  lmstats$feature <- as.character()
  flag <- 0
  for (vars in 1:ncol(DATATAB)) {
    flag <- flag + 1
    mod0 <- lm(rank(DATATAB[, vars]) ~ rank(MET$Calprotectin) + rank(MET$Moisture) + rank(MET$BMI))
    mod1 <- lm(rank(DATATAB[, vars]) ~ rank(MET$Calprotectin) + rank(MET$Moisture) + rank(MET$BMI) + (MET$CRCtype))
    summod <- summary(mod1)
    ANtest <- anova(mod0, mod1)
    lmstats[flag, "dataset"] <- as.character(named[i])
    lmstats[flag, c("dataset", "feature")] <- c(as.character(named[i]), as.character(colnames(DATATAB)[vars]))
    lmstats[flag, c("calprotectin_estimate", "calprotectin_std.Error", "calprotectin_t.value", "calprotectin_P")] <- summod$coefficients[2, c(1, 2, 3, 4)]
    lmstats[flag, c("Moisture_estimate", "Moisture_std.Error", "Moisture_t.value", "Moisture_P")] <- summod$coefficients[3, c(1, 2, 3, 4)]
    lmstats[flag, c("BMI_estimate", "BMI_std.Error", "BMI_t.value", "BMI_P")] <- summod$coefficients[4, c(1, 2, 3, 4)]
    lmstats[flag, c("CRCtype_estimate", "CRCtype_std.Error", "CRCtype_t.value", "CRCtype_P")] <- summod$coefficients[5, c(1, 2, 3, 4)]
    lmstats[flag, c("anova.F", "anova_P")] <- c(ANtest$F[2], ANtest$`Pr(>F)`[2])
    lmstats[flag, "N"] <- nrow(DATATAB)
  }
  lmstats[, "calprotectin_AdjP"] <- p.adjust(as.numeric(as.character(lmstats$calprotectin_P)), "BH")
  lmstats[, "Moisture_AdjP"] <- p.adjust(as.numeric(as.character(lmstats$Moisture_P)), "BH")
  lmstats[, "BMI_AdjP"] <- p.adjust(as.numeric(as.character(lmstats$BMI_P)), "BH")
  lmstats[, "anova_AdjP"] <- p.adjust(as.numeric(as.character(lmstats$anova_P)), "BH")
  lmstatsall_CMB_QMP <- rbind(lmstatsall_CMB_QMP, lmstats)
}

# Subset features with adjusted p-value < 0.05
subset(lmstatsall_CMB_QMP, anova_AdjP < 0.05)

QMP_LCMP_acBCM <- lmstatsall_CMB_QMP


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Perform glm analysis for RMP data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj_phylo <- RMP_for_glm
species_RMP <- data.frame(otu_table(obj_phylo))
mata_crc <- data.frame(sample_data(obj_phylo))
mata_crc <- mata_crc[, c("Calprotectin", "CRC_general_status")]
setnames(mata_crc, c("Calprotectin", "CRC_general_status"), c("calpro", "CRCtype"))
mata_crc$CRCtype <- as.factor(mata_crc$CRCtype)
levels(mata_crc$CRCtype) <- c(1, 2, 3)
mata_crc$CRCtype <- as.numeric(as.character(mata_crc$CRCtype))
COHORTDATA <- mata_crc[, c("calpro", "CRCtype")]

# Check data types and dimensions
is.numeric(COHORTDATA$calpro)
is.numeric(COHORTDATA$CRCtype)
dim(COHORTDATA)
dim(na.omit(COHORTDATA))


# Perform glm analysis for each feature in RMP data
subsamps <- row.names(species_RMP)
datasets <- list(species_RMP)
named <- "species_RMP"
lmstatsall_RMP <- data.frame()

for (i in 1:length(datasets)) {
  MAT2 <- datasets[[i]]
  MATBIN <- MAT2[subsamps, ]
  MATBIN[MATBIN > 0] <- 1
  sub <- names(which((colSums(MATBIN) / nrow(MATBIN)) > 0.000000050))
  sub <- setdiff(sub, sub[grep("unclassified", sub)])
  MAT2 <- MAT2[, which(colnames(MAT2) %in% sub)]
  DATATAB <- MAT2
  MET <- COHORTDATA[subsamps, ]
  lmstats <- data.frame("dataset" = character(), "feature" = character(),
                        "calprotectin_estimate" = numeric(), "calprotectin_std.Error" = numeric(), "calprotectin_t.value" = numeric(), "calprotectin_P" = numeric(),
                        "anova.F" = numeric(), "anova_P" = numeric(), "N" = numeric())
  lmstats$dataset <- as.character()
  lmstats$feature <- as.character()
  flag <- 0
  for (vars in 1:ncol(DATATAB)) {
    flag <- flag + 1
    mod0 <- lm(rank(DATATAB[, vars]) ~ rank(MET$calpro))
    mod1 <- lm(rank(DATATAB[, vars]) ~ rank(MET$calpro) + (MET$CRCtype))
    summod <- summary(mod1)
    ANtest <- anova(mod0, mod1)
    lmstats[flag, "dataset"] <- as.character(named[i])
    lmstats[flag, c("dataset", "feature")] <- c(as.character(named[i]), as.character(colnames(DATATAB)[vars]))
    lmstats[flag, c("calprotectin_estimate", "calprotectin_std.Error", "calprotectin_t.value", "calprotectin_P")] <- summod$coefficients[2, c(1, 2, 3, 4)]
    lmstats[flag, c("anova.F", "anova_P")] <- c(ANtest$F[2], ANtest$`Pr(>F)`[2])
    lmstats[flag, "N"] <- nrow(DATATAB)
  }
  lmstats[, "calprotectin_AdjP"] <- p.adjust(as.numeric(as.character(lmstats$calprotectin_P)), "BH")
  lmstats[, "anova_AdjP"] <- p.adjust(as.numeric(as.character(lmstats$anova_P)), "BH")
  lmstatsall_RMP <- rbind(lmstatsall_RMP, lmstats)
}

# Subset features with adjusted p-value < 0.05
subset(lmstatsall_RMP, anova_AdjP < 0.05)

RMP_LCMP_acBCM <- lmstatsall_RMP


# Perform glm analysis for RMP data including additional covariates (Moisture, BMI)
species_RMP = (data.frame(otu_table(obj_phylo)))
mata_crc = data.frame(sample_data(obj_phylo))
head(mata_crc)
mata_crc <- mata_crc[, c("CRC_general_status", "Calprotectin", "Moisture", "BMI")]
setnames(mata_crc, c("CRC_general_status", "Calprotectin", "Moisture", "BMI"),  c("CRCtype", "Calprotectin","Moisture", "BMI"))
head(mata_crc)
mata_crc$CRCtype  = as.factor(mata_crc$CRCtype)
levels(mata_crc$CRCtype) = c(1,2,3)
mata_crc$CRCtype = as.numeric(as.character(mata_crc$CRCtype))
#mata_crc$CRCtype = as.factor(as.character(mata_crc$CRCtype))
COHORTDATA <- mata_crc

###
is.numeric(COHORTDATA$Calprotectin)
is.numeric(COHORTDATA$Moisture)
is.numeric(COHORTDATA$BMI)
is.numeric(COHORTDATA$CRCtype)
dim(COHORTDATA)
dim((na.omit(COHORTDATA)))

dim(species_RMP)


row.names(species_RMP)
subsamps=row.names(species_RMP)
datasets=list(species_RMP)
named=c("species_RMP")

lmstatsall_CMB_RMP=data.frame()

for(i in 1:length(datasets)){
  MAT2=datasets[[i]]
  print(dim(MAT2))
  MATBIN=MAT2[subsamps,]
  MATBIN[MATBIN>0]=1
  sub=names(which((colSums(MATBIN)/nrow(MATBIN))>0.000000050))
  sub=setdiff(sub,sub[grep("unclassified",sub)])
  MAT2=MAT2[,which(colnames(MAT2)%in%sub)]
  DATATAB=MAT2
  MET=COHORTDATA[subsamps,]
  lmstats=data.frame("dataset"=character(),"feature"=character(),
                     "calprotectin_estimate"=numeric(),"calprotectin_std.Error"=numeric(),"calprotectin_t.value"=numeric(),"calprotectin_P"=numeric(),
                     "Moisture_estimate"=numeric(),"Moisture_std.Error"=numeric(),"Moisture_t.value"=numeric(),"Moisture_P"=numeric(),
                     "BMI_estimate"=numeric(),"BMI_std.Error"=numeric(),"BMI_t.value"=numeric(),"BMI_P"=numeric(),
                     "anova.F"=numeric(),"anova_P"=numeric(),"N"=numeric())
  lmstats$dataset=as.character()
  lmstats$feature=as.character()
  flag=0
  for(vars in 1:ncol(DATATAB)){
    flag=flag+1
    mod0=lm(rank(DATATAB[,vars]) ~ rank(MET$Calprotectin) + rank(MET$Moisture) + rank(MET$BMI))
    mod1=lm(rank(DATATAB[,vars]) ~ rank(MET$Calprotectin) + rank(MET$Moisture) + rank(MET$BMI) + (MET$CRCtype))
    summod=summary(mod1)
    ANtest=anova(mod0,mod1)
    lmstats[flag,"dataset"]=as.character(named[i])
    lmstats[flag,c("dataset","feature")]=c(as.character(named[i]),as.character(colnames(DATATAB)[vars]))
    lmstats[flag,c("calprotectin_estimate","calprotectin_std.Error","calprotectin_t.value","calprotectin_P")]=summod$coefficients[2,c(1,2,3,4)]
    lmstats[flag,c("Moisture_estimate","Moisture_std.Error","Moisture_t.value","Moisture_P")]=summod$coefficients[3,c(1,2,3,4)]
    lmstats[flag,c("BMI_estimate","BMI_std.Error","BMI_t.value","BMI_P")]=summod$coefficients[4,c(1,2,3,4)]
    lmstats[flag,c("CRCtype_estimate","CRCtype_std.Error","CRCtype_t.value","CRCtype_P")]=summod$coefficients[5,c(1,2,3,4)]
    lmstats[flag,c("anova.F","anova_P")]=c(ANtest$F[2],ANtest$`Pr(>F)`[2])
    lmstats[flag,"N"]=nrow(DATATAB)
  }
  lmstats[,"calprotectin_AdjP"]=p.adjust(as.numeric(as.character(lmstats$calprotectin_P)),"BH")
  lmstats[,"Moisture_AdjP"]=p.adjust(as.numeric(as.character(lmstats$Moisture_P)),"BH")
  lmstats[,"BMI_AdjP"]=p.adjust(as.numeric(as.character(lmstats$BMI_P)),"BH")
  lmstats[,"anova_AdjP"]=p.adjust(as.numeric(as.character(lmstats$anova_P)),"BH")
  lmstatsall_CMB_RMP=rbind(lmstatsall_CMB_RMP,lmstats)
}

print(subset(lmstatsall_CMB_RMP, anova_AdjP<0.05))

RMP_LCMP_acBCM <- lmstatsall_CMB_RMP
