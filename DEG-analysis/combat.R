# Set up

setwd("~/RNA/oct")
plot_dir <- c("~/RNA/oct/plots") 

#load packages #
library(dplyr)
library(ggplot2)
library(ggrepel) 
library(tidyr)
library(tidylog)
library(limma)
library(stringi)
library(janitor)
library(data.table)
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
                '_e2_c3', 'e3_c3', 'e5_c5'),replacement = c(""),
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

## Merge studies (start here)---------------------------------------------------

load(file = "VAST_Tyger_oct.RData")

# Option, exclude vaccinated participants:

idx <- (vast$samples$study_arm == "CTRL")
vast <- vast[,idx]


counts <- cbind(tyger$counts,vast$counts)
samples <- rbind(tyger$samples,vast$samples)
samples$time <- as.factor(samples$time)

table(samples$sequence_pool)
table(samples$study_arm)
table(samples$time_point3)


# Filter contaminated
# should filter for contaminated samples if I'm going to just use the uncontaminated 
samples <- filter(samples, contaminated_sample == "No")

# remove baseline samples in original pool who have lab ids in check
check <- filter(samples, sequence_pool == "resequence" & time == "0")
check <- filter(samples, (time == "0" & sequence_pool == "original" & lab_id %in% check$lab_id))
samples <- filter(samples, rownames(samples) %!in% rownames(check))

# Clean diagnosis var
table(samples$diagnosis)
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
table(samples$time_point)
# D0.12h <- D0+12
# D0 <- D00
samples <- samples %>%
  mutate(time_point3 = ifelse(time_point3 == "D0+12h", "D0.12h", time_point3)) %>%
  mutate(time_point3 = ifelse(time_point3 == "D00", "D0", time_point3))

# Make day of challenge V28 for vaccine samples day 0
samples <- samples %>%
  mutate(time_point = ifelse(time_point3 == "D0" & (study_arm %in% c("ViTCV", "ViPS")), "V28", time_point3))

# Make new baseline time point, V0 for vaccinees and D0 for challenge
samples <- samples %>%
  mutate(time_point = ifelse(time_point %in% c("V0","D0"), "Baseline", time_point)) 

# Day 7 time point, around D7 for TD
samples <- samples %>% 
  mutate(time_point = ifelse(time %in% c("6","7","8") & (time_point3 == "TD"), "D7", time_point))

# samples <- samples %>% 
# mutate(time_point = ifelse(time_point == "D7" & batch == "3", "D7", time_point))

table(samples$time_point, samples$diagnosis)
# filter vaccine timepoints and D01-14 which are study specific

samples <- filter(samples, time_point %in% c("Baseline", "D0.12h", "TD"))

# Just D7
# samples <- filter(samples, time_point == "D7")
tab <- table(samples$time_point, samples$batch_pool, samples$diagnosis)
addmargins(tab)

counts <- counts[ , colnames(counts) %in%  rownames(samples)]

# Batch correction -------------------------------------------------------------

# covar_mod	
# Model matrix for multiple covariates to include in linear model (signals from these variables are kept in data after adjustment
# Can't be confounded with batch tho, therfore need to rework batch or leave out 

samples$arm_time <- str_c(samples$study_arm, samples$time_point, sep = "_")

# Merge TN and control samples
samples <- samples %>% 
  mutate(chall_vax = ifelse(str_detect(study_arm, "TN|CTRL"), "T", study_arm))

# covar <- samples %>% dplyr::select(arm_time)
covar <- samples %>% dplyr::select(time_point)

table(samples$batch_pool, samples$chall_vax)

 #combat <- sva::ComBat_seq(counts=counts,
                          #batch=samples$batch,
                          #group=samples$chall_vax)

combat <- sva::ComBat_seq(counts=counts,
                          batch=samples$batch_pool,
                          group=samples$time_point)


# 7/nov/22group is the var we want to see an effect in, so vaccine vs control, but we also need time point, unless running bc separate
# 31/oct/22
# It's worse with seq_pool in as unbalanced, therefore just correct on batch here, add seq_pool as covariate later


## PCA -------------------------------------------------------------------------
# Show that control samples look similar, batch independent
# Show that we can still se biological variation eg. TD vs nTD at D7

# see batch_pcas.R

#### Calculate TPM ------------------------------------------------------------

keep.exprs <- filterByExpr(combat, group = samples$time_point)

combat <- combat[keep.exprs,] #subset og data
dim(combat)

tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

gene_info <- vast$genes[vast$genes$gene_id %in% row.names(combat),c(1,6)]
names.use <- rownames(combat) #order
gene_info <- gene_info[names.use, 2]

tpm <- tpm3(combat, gene_info)
exprs <- as.data.frame(t(tpm))
exprs$row_names <- rownames(exprs)
samples$row_names <- rownames(samples)
pheno_exprs <- left_join(samples, exprs, by = "row_names")


### Save data ----------------------------------------------------------------
save.image(file = "VAST_Tyger_nov2.RData")

# gene info
gene_info <- vast$genes[vast$genes$gene_id %in% row.names(combat),c(1,6)]
names.use <- rownames(combat) # order
gene_info <- gene_info[names.use, 2]

vast_tyg_novax <- DGEList(counts =combat, samples=samples, genes = vast$genes)
saveRDS(vast_tyg_novax, file = "VT_novax.Rds")

dim(combat)
dim(vast$genes)

load(file = "VAST_Tyger_oct.RData")


# write.csv(pheno_exprs, file = "combat_vast_tyger_d7_tpm.csv")

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
