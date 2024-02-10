## DADA2 Analysis Pipeline

## This code represents an analysis pipeline for processing and analyzing DNA sequencing data using the DADA2 package in R. 
## The pipeline includes steps for quality filtering, denoising, merging paired-end reads, chimera removal, and taxonomic assignment.

## Load required library
set.seed(531) 
library(dada2)
packageVersion('dada2')

## Specify the path for input FASTQ files (replace with the appropriate path after downloading the files from EGA)
path <- '/FASTQ/demultiplexed_run_DI18R24/'
fns <- list.files(path)
fastqs <- fns[grepl('.fq.gz$', fns)]
fastqs <- sort(fastqs)

## Separate forward and reverse read files and extract sample names
fnFs <- fastqs[grepl('.1.fq.gz', fastqs)]
fnRs <- fastqs[grepl('.2.fq.gz', fastqs)]
sample.names <- sapply(strsplit(fnFs, '.1.fq.gz'), `[`, 1)

## Fully specify file paths for forward and reverse reads
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

## Set up filtered file paths and perform filtering and trimming
filt_path <- file.path(path, 'filtered')
filtFs <- file.path(filt_path, paste0(sample.names, '_F_filt.fastq.gz'))
filtRs <- file.path(filt_path, paste0(sample.names, '_R_filt.fastq.gz'))

## MiSeq (LCPM cohort) and HiSeq (FGFP cohort) tested for my primer constructs
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(130, 200), trimLeft = c(30, 30), maxN = 0, maxEE = c(2, 2), truncQ = 11, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

## Create a data frame with filtering results
DFout <- data.frame(out)
## Select samples with at least 50 reads
sample.names0 <- sapply(strsplit(rownames(subset(DFout, reads.out > 50)), '.1.fq.gz'), `[`, 1)
filtFs0 <- file.path(filt_path, paste0(sample.names0, '_F_filt.fastq.gz'))
filtRs0 <- file.path(filt_path, paste0(sample.names0, '_R_filt.fastq.gz'))

## Learn forward and reverse error rates
errF <- learnErrors(filtFs0, nread = 1e6, multithread = TRUE)
errR <- learnErrors(filtRs0, nread = 1e6, multithread = TRUE)

## Dereplicate reads
derepRs <- derepFastq(filtRs0, verbose = TRUE)
derepFs <- derepFastq(filtFs0, verbose = TRUE)

## Assign names to dereplicated objects
names(derepFs) <- sample.names0
names(derepRs) <- sample.names0

## Perform denoising and merge paired-end reads
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

## Merge the denoised reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

## Create a sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
seqtabHS <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)

## Save the sequence table
saveRDS(seqtabHS, "DI18R24.seqtab.rds")

###### Perform the same analysis for the other sets of fastq files from miseq runs: DI18R36 DI18R37 DI19R04 DI19R05 DI19R06

## Merge sequence tables from the six runs 
dada2_seqtabExp1 <- readRDS('DI18R24.seqtab.rds')
dada2_seqtabExp2 <- readRDS('DI18R36.seqtab.rds')  
dada2_seqtabExp3 <- readRDS('DI18R37.seqtab.rds')
dada2_seqtabExp4 <- readRDS('DI19R04.seqtab.rds')
dada2_seqtabExp5 <- readRDS('DI19R05.seqtab.rds')
dada2_seqtabExp6 <- readRDS('DI19R06.seqtab.rds')
seqtab_all <- mergeSequenceTables(dada2_seqtabExp1, dada2_seqtabExp2, dada2_seqtabExp3, dada2_seqtabExp4, dada2_seqtabExp5, dada2_seqtabExp6)

## Remove chimeras
seqtab_all_no_chimeras <- removeBimeraDenovo(seqtab_all, method = "consensus", multithread = TRUE)
saveRDS(seqtab_all_no_chimeras, "seqtab_all_no_chimeras.rds")

## Assign taxonomy

## Assign taxonomy using the SILVA database (replace with your own reference file)
seqtab_all_no_chimeras_tax_SLV138.1 <- assignTaxonomy(seqtab_all_no_chimeras, "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread = TRUE)
seqtab_all_no_chimeras_tax_SLV138.1sp <- assignSpecies(seqtab_all_no_chimeras, "silva_species_assignment_v138.1.fa.gz")
seqtab_all_no_chimeras_tax_SLV138.1[,7] <- seqtab_all_no_chimeras_tax_SLV138.1sp[,7]
saveRDS(seqtab_all_no_chimeras_tax_SLV138.1, "seqtab_all_no_chimeras_tax_SLV138.1.rds")
