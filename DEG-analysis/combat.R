# Set up

setwd("~/RNA")
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

# Harmonise time points

# recode time_point3
# D0.12h <- D0+12
# D0 <- D00

samples <- samples %>%
  mutate(time_point3 = ifelse(time_point3 == "D0+12h", "D0.12h", time_point3)) %>%
  mutate(time_point3 = ifelse(time_point3 == "D00", "D0", time_point3))

samples <- samples %>%
  mutate(time_point = ifelse(time_point3 == "V0" & (study_arm == "ViPS" | study_arm == "ViTCV") , "Baseline", time_point3)) %>%
  mutate(time_point = ifelse(time_point3 == "D0" & (study_arm == "CTRL" | study_arm == "TN"), "Baseline", time_point))

table(samples$time_point)
# merge Day 7 TD time point where poss, as TD have a separate time point
# make a separate diagnoisis time var maybe, time points from 8 days plus are TD
# TD 14 is just TN and control samples

# TD_D7plus
table(samples$time_point)
# merge TD24/12hs into one timepoint as combat doesnt like single obvs
# D08,10,12 Join

samples <- samples %>%
mutate(time_point = ifelse(time == "7", "D7", time_point)) %>%
  mutate(time_point = ifelse(time_point3 == "TD.12h", "TD.24h", time_point))  %>%
  mutate(time_point = ifelse(time_point3 %in% c("D08", "D10","D12"), "D08_12",time_point))

# Clean diagnosis var
samples <- samples %>%
  mutate(diagnosis = ifelse(diagnosis == "1", "TD", "nTD"))

# Make new batch sequence pool var
samples <- samples %>%
  mutate(batch = ifelse(seq_run == "original", "1","2"))
samples$batch_pool <- str_c(samples$batch, samples$sequence_pool, sep = "_")
table(samples$batch_pool)

# should filter for contaminated samples if I'm going to just use the uncontaminated 
samples <- filter(samples, contaminated_sample == "No")
counts <- counts[ , colnames(counts) %in%  rownames(samples)]


# Batch correction -------------------------------------------------------------

combat <- sva::ComBat_seq(counts=counts,
                          batch=samples$batch_pool,
                          group=samples$time_point)

# covar_mod	
# Model matrix for multiple covariates to include in linear model (signals from these variables are kept in data after adjustment
# Can't be confounded with batch tho, therfore need to rework batch or leave out 

covar <- samples %>% dplyr::select(study_arm)
## PCA -------------------------------------------------------------------------
# Show that control samples look similar, batch independent
# Show that we can still se biological variation eg. TD vs nTD at D7

# 1) Make df with counts and vars joined for easy data subestting
count_t <- t(counts)
data <- bind_cols(samples, count_t)

comb_t <- t(combat)
data2 <- bind_cols(samples, comb_t)

# 2a) PCA of uncorrected seq data

table(data$time_point)
pca <- filter(data, time_point == "D7")
pca_1 <- prcomp(pca[,-c(1:20)], scale. = TRUE)

# 3a) Plots
fviz_eig(pca_1)

a <- autoplot(pca_1, data = pca,
         colour = "study_arm", shape="batch",
         frame = TRUE, frame.type = 'norm') +
  theme_bw() + ggtitle("Raw Data")
a

# 2b) PCA batch corrected data
pca <- filter(data2, time_point3 == "TD")
pca_2 <- prcomp(pca[,-c(1:20)], scale. = TRUE)

# 3b) Plots
fviz_eig(pca_2)
b <- autoplot(pca_2, data = pca,
         colour = "study_arm", shape="batch",
         frame = TRUE, frame.type = 'norm') +
  theme_bw() + ggtitle("Batch Corrected")
b

ggsave(filename = paste(plot_dir, "TD_PCA.png"))
(a + b) + plot_layout(guides = "collect") + plot_annotation(title = "TD")


# Batch 1 = TyGER original
# Batch 2 = TyGER resequenced
# Batch 3 = VAST

# Vaccine differences
vast_pc <- filter(data, batch == "3" & time_point3 == "TD")
pca <- prcomp(vast_pc[,-c(1:19)], scale. = TRUE)
a <- autoplot(pca, data = vast_pc,
         colour = "study_arm", shape="batch_pool",
         frame = TRUE, frame.type = 'norm') +
  theme_bw() + ggtitle("Raw Data")
a

table(vast$samples$time_point3)

#"batch_pool", shape="study_arm"
#"study_arm", shape="batch_pool"


vast_pc <- filter(data2, batch == "3" & time_point3 == "TD")
 pca <- prcomp(vast_pc[,-c(1:19)], scale. = TRUE)
b <- autoplot(pca, data = vast_pc,
         colour = "study_arm", shape="batch_pool",
         frame = TRUE, frame.type = 'norm') +
  theme_bw() + ggtitle("Batch Corrected")
b

(a + b) + plot_layout(guides = "collect") + plot_annotation(title = "VAST TD")


write.csv(pheno_exprs, file = "combat_vast_tyger.csv")




?autoplot


# Extra -------------------------------------------------------------------

#Convert counts to log to minimise the effect of small values and negatives
cpm <- cpm(data)
lcpm <- cpm(data, log=TRUE) #Used for exploratory plots/checks
L <- mean(data$samples$lib.size) * 1e-6 #average library size
M <- median(data$samples$lib.size) * 1e-6
c(L, M)

#Filter low exprs genes
keep.exprs <- filterByExpr(data)
# Warning - all samples appear to belong to the same group
data <- data[keep.exprs, keep.lib.sizes=FALSE] 

#TMM normalisation already performed
# <- calcNormFactors(data, method = "TMM")
tyger <- data
vast <- data
rm(data)
rm(cpm)
rm(lcpm)

save.image(file = "TyGER_oct.RData")


## Filter out contaminated samples, and time_point -----------------------------
# Add exprs colnames to samples for matching tables by

data <- tyger
samples <- data$samples
counts <- data$counts
genes <- data$genes

samples$exprs_cols <- colnames(counts)

table(data$samples$time_point3)
samples$visit <- str_replace_all(samples$time_point3, "D0.12h", "12h")
table(samples$visit)

time_point <- c("12h")
v2 <- c("contaminated_sample")

# 1) Filter time point and contaminated
samples <-
  samples %>%
  filter(visit == time_point & contaminated_sample != v2)

# 1b) Filter time - vast
samples <-
  samples %>%
  filter(visit == time_point)

counts <-
  counts[,colnames(counts) %in% samples$exprs_cols]
dim(counts)


#reform s3 DGElist object
tyger_12h <-DGEList(counts = counts,
                    genes = genes,
                    samples = samples)

vast_12h <- DGEList(counts = counts,
                    genes = genes,
                    samples = samples)

rm("data","samples", "counts", "keep.exprs", "genes")

save.image(file = "VAST_TYGER.RData")
load(file = "VAST_TYGER.RData")

table(vast$samples$time_point3)
table(tyger$samples$time_point3)

library(ganalyse)
test <- as.eset(tyger_12h)
test2 <- as.eset(vast_12h)

merged <- merge(test1, test2, method = "combat")

## Keep only common genes ------------------------------------------------------
common <- intersect(tyger_12h$genes$gene_id, vast_12h$genes$gene_id)
length(common)

vast_12h$genes <- vast_12h$genes %>% filter(gene_id %in% common)
vast_12h$counts <- vast_12h$counts[row.names(vast_12h$counts) %in% common,]
dim(vast_12h$counts)

tyger_12h$genes <- tyger_12h$genes %>% filter(gene_id %in% common)
tyger_12h$counts <- tyger_12h$counts[row.names(tyger_12h$counts) %in% common,]

dim(tyger_12h$samples)

tyger_12h$samples$batch <- rep("tyger",nrow(tyger_12h$samples))
vast_12h$samples$batch <- rep("vast",nrow(vast_12h$samples))

# order genes
names.use <- vast_12h$genes$gene_id
tyger_12h$genes <- tyger_12h$genes[names.use,]
# order counts
tyger_12h$counts <-t yger_12h$counts[names.use,]

# Common sample data
common <- intersect(names(tyger_12h$samples), names(vast_12h$samples))

View(data$samples)

test <- cbind

## Tyger checking ---------------------------------------------------------------
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
