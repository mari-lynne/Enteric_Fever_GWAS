# Set up

setwd("~/RNA/oct")
plot_dir <- c("~/RNA/oct/plots") 

#load packages #
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel) 
library(tidyr)
library(tidylog)
library(limma)
library(stringi)
library(janitor)
library(stringr)
library(splines)
library(edgeR)
library(BiocManager)
library(DESeq2)
library(Glimma)
library(RColorBrewer)
library(ggfortify)
library(factoextra)
library(patchwork)

# Steps Filter D012h/24h time points from vast and Tyger datasets
# Merge using combat

# Pre-process -----------------------------------------------------------------

'%!in%' <- function(x,y)!('%in%'(x,y))

## VAST ------------------------------------------------------------------------

load("/home/mari/RNA/Filter2_VAST_STAR_HTSeq_gene_meta_autosomes_mismatch_corrected_demo_minus_rRNA_globins_autosomes_2021-04-27.R")

data <- VAST_autosomes #RNAseq data set
rm(VAST_autosomes)

#Clean names using janitor package
str(data$meta_data)
data$meta_data <- clean_names(data$meta_data)
#Rename variables
data$meta_data<- data$meta_data %>% dplyr::rename(time = days_since_challenge)
#Clean var names
colnames(data$meta_data) <-
  stri_replace_all_regex(
    colnames(data$meta_data),
    pattern = c('_x', 'x_',
                '_e2_c3', 'e3_c3', 'e5_c5'),replacement = c(''),
    vectorize = FALSE)
data$samples$batch <- rep("3",nrow(data$samples))
data$meta_data <- data$meta_data %>% select(-group)

data$samples <- bind_cols(data$samples, data$meta_data)

vast <- data

## Tyger -------------------------------------------------------------------
load(file = "Filter2_TyGER_combined_RNAseq_runs_STAR_HTSeq_autosome_gene_meta.HLA_genotype_corrected_minus_rRNA_globins_demo_autosomes_2021-05-18.R")

data <- TyGER_combined_data #VAST_autosomes #RNAseq data set
rm(TyGER_combined_data)


#Clean names using janitor package
str(data$metadata)
data$metadata<- clean_names(data$metadata)
#Rename variables
data$metadata <- data$metadata %>% dplyr::rename(time = days_since_challenge)
data$metadata <- data$metadata %>% dplyr::rename(strain = allocation)

# Mutliple entries in sample + meta for norm.factors
str(data$samples)
data$samples <- data$samples %>% dplyr::rename(samp_process_name = names_preprocessed_under)
data$metadata <- data$metadata %>% select(-group, -visit) %>%
  dplyr::rename(diagnosis = diagnosed_x)

data$metadata <- data$metadata %>%
  mutate(batch = ifelse(seq_run == "original", "1","2"))

data$samples <- bind_cols(data$samples, data$metadata)
str(data$samples)

tyger <- data

# Tidy metadata ---------------------------------------------------------------

# Intersecting metacolumns
check <- intersect(names(tyger$samples), names(vast$samples))
check

# what extra vars do we need:
# VAST - sequence_pool = seq_run
# Tyger - strain, contaminated_sample, seq_run
tyger$samples <- tyger$samples %>% dplyr::rename(sequence_pool = seq_run)
tyger$samples <- tyger$samples %>% dplyr::rename(study_arm = strain)

levels(tyger$samples$study_arm) <- c("","TN", "CTRL")
vast$samples$contaminated_sample <- rep("No", nrow(vast$samples))

vast$samples <- vast$samples[,colnames(vast$samples) %in% check]
tyger$samples <- tyger$samples[,colnames(tyger$samples) %in% check]

# Tidy countdata ---------------------------------------------------------------

## Keep only common genes ------------------------------------------------------
common <- intersect(tyger$genes$gene_id, vast$genes$gene_id)
length(common)

vast$genes <- vast$genes %>% filter(gene_id %in% common)
vast$counts <- vast$counts[row.names(vast$counts) %in% common,]
dim(vast$counts)

tyger$genes <- tyger$genes %>% filter(gene_id %in% common)

tyger$counts <- tyger$counts[row.names(tyger$counts) %in% common,]
dim(tyger$counts)

# order genes
names.use <- vast$genes$gene_id
tyger$genes <- tyger$genes[names.use,]
# order counts
tyger$counts <- tyger$counts[names.use,]

rm(data)

save.image(file = "VAST_Tyger_oct.RData")

## Merge studies ----------------------------------------------------------------

load(file = "VAST_Tyger_oct.RData")

counts <- cbind(tyger$counts,vast$counts)
samples <- rbind(tyger$samples,vast$samples)
samples$time <- as.factor(samples$time)

table(samples$sequence_pool)
table(samples$study_arm)
table(samples$time_point3)

# Filter contaminated
# should filter for contaminated samples if I'm going to just use the uncontaminated 
samples <- filter(samples, contaminated_sample == "No")

# Clean diagnosis var
table(samples$diagnosis)

no_chall <- filter(samples, diagnosis == "UNKNOWN - Not Challenged")
samples <- filter(samples, diagnosis != "UNKNOWN - Not Challenged")

# Clean diagnosis var
table(samples$diagnosis)
samples <- samples %>%
  mutate(diagnosis = ifelse(diagnosis == "1", "TD", diagnosis),
         diagnosis = ifelse(diagnosis == "0", "nTD", diagnosis))

# Update counts 
counts <- counts[ , colnames(counts) %in%  rownames(samples)]

# Make new batch sequence pool var
samples$batch_pool <- str_c(samples$batch, samples$sequence_pool, sep = "_")
table(samples$batch_pool)

# Harmonise time points -------------------------------------------------------

# recode time_point3
# D0.12h <- D0+12
# D0 <- D00
table(samples$time_point)

samples <- samples %>%
  mutate(time_point3 = ifelse(time_point3 == "D0+12h", "D0.12h", time_point3)) %>%
  mutate(time_point3 = ifelse(time_point3 == "D00", "D0", time_point3))

# Make day of challenge V28 for vaccine samples day 0
samples <- samples %>%
  mutate(time_point = ifelse(time_point3 == "D0" & (study_arm %in% c("ViTCV", "ViPS")), "V28", time_point3))


# Make new baseline time point, V0 for vaccinees and D0 for challenge
samples <- samples %>%
  mutate(time_point = ifelse(time_point %in% c("V0","D0"), "Baseline", time_point)) 
# filter vaccine timepoints and D01-14 which are study specific
samples <- filter(samples, time_point %in% c("Baseline", "D0.12h", "TD"))

tab <- table(samples$time_point, samples$batch, samples$diagnosis)
addmargins(tab)

counts <- counts[ , colnames(counts) %in%  rownames(samples)]



# Batch correction -------------------------------------------------------------

# covar_mod	
# Model matrix for multiple covariates to include in linear model (signals from these variables are kept in data after adjustment
# Can't be confounded with batch tho, therfore need to rework batch or leave out 

samples$arm_time <- str_c(samples$study_arm, samples$time_point, sep = "_")

covar <- samples %>% dplyr::select(arm_time)
table(samples$batch, samples$study_arm)

combat <- sva::ComBat_seq(counts=counts,
                          batch=samples$batch,
                          group=samples$time_point)


# It's worse with seq_pool in as unbalanced
# 


## PCA -------------------------------------------------------------------------
# Show that control samples look similar, batch independent
# Show that we can still se biological variation eg. TD vs nTD at D7

# see batch_pcas.R


### Save data ----------------------------------------------------------------
save.image(file = "VAST_Tyger_oct.RData")
load(file = "VAST_Tyger_oct.RData")

pheno_exprs <- bind_cols(samples, comb_t)

write.csv(pheno_exprs, file = "combat_vast_tyger.csv")

# So far batches look better (although not fully corrected), biological affects not changed when comparing within a batch



## Extra ------------------------------------------------------------------------

### Tyger checking ---------------------------------------------------------------
samples <- data$samples
genes <- data$genes
metadata <- data$metadata
counts <- data$counts

miss <- anti_join(samples, metadata, by = "names_preprocessed_under")
check <- bind_cols(samples$lib.size,metadata$lib_size,samples$sequence_pool,metadata$contaminated_sample)

# remove resequence for n
names(check) <- c("samp", "meta", "seq","contam")
check %>% ggplot(aes(x=samp,y=meta,colour=contam)) + geom_point()

reseq <- filter(samples, contaminated_sample != "No") # Include as an extra batch?

check <- samples[samples$lib.size >4e7,]
# they are samples with seqids (resquenced?)
