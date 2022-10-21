#Library 
library(dplyr)
library(tidyr)
library(limma)
library(stringr)
library(data.table)
library(ggplot2)
library(forcats)
library(factoextra)
library(tidylog)
library(ggpubr)
library(tibble)
library(stringr)
library(MatrixEQTL)
library(Biobase)
library(janitor)
library(biomaRt)

#load(file = "~/RNA/oct/meqtl_TD.RData")

#Functions ####
'%!in%' <- function(x,y)!('%in%'(x,y))

# Directories and variables ---------------------------------------------------

setwd("~/GWAS_22/gwas_final/eQTL")
plot_dir <- c("~/GWAS_22/gwas_final/eQTL/plots") 
time_point <- c("12h")
study <- c("/T1T2_")
out_dir <- getwd()

id_data <- c("~/GWAS_22/gwas_final/meta/T1T2_ID_link.csv")
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/QC/typhoid2.IBD.fam") #fam file
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/tophits_typhoid2.txt")
covar_data <- c("~/GWAS_22/gwas_final/meta/typhoid_pcacovar.txt")

# Load micro-array experiment data
load("T1_T2_with_rsn.norm.R")
data = T1_T2_autosomes
rm(T1_T2_autosomes)

exprs <- as.data.frame(exprs(data))
pData <- clean_names(setDT(data@phenoData@data, keep.rownames = TRUE)[])
fData <- clean_names(data@featureData@data)

# Set up new vars -------------------------------------------------------------
# Make new typhoid diagnosis variable, na's = 0
pData <- pData %>%
  mutate(diagnosis = ifelse(is.na(
    pData$day_of_typhoid_diagnosis_from_challenge
  ), 0, 1))

# Rename vars
names(pData)[names(pData) == "days_since_challenge"] <- "time"
names(pData)[names(pData) == "timepoint3"] <- "visit"

# Time point
table(pData$visit)

pData$visit <- str_replace_all(pData$visit, "D0\\+12", "12h")

# Filter DGE data for time point and genotyped samples ------------------------

# Add exprs colnames to pData for matching tables by
pData$exprs_cols <- colnames(exprs)

# 1) Filter time point

pData <-
  pData %>%
  filter(visit == time_point)

exprs <-
  exprs[,colnames(exprs) %in% pData$exprs_cols]

# 2) Filter for gwas samples

# Match by intersecting participant ids

IDlink <- clean_names(read.csv(file = id_data))
#Rename ID col so can join appropriately
IDlink <- dplyr::rename(IDlink, part_number = participant_id)

# Genotyping data in DGE data
geno_ids <-
  IDlink %>%
  filter(part_number %in% pData$part_number) # Do missing samples have D1

# Filter DGE data for appropriate genotyped samples
pData <- 
  pData %>%
  filter(part_number %in% geno_ids$part_number)
exprs <-
  exprs[,colnames(exprs) %in% pData$exprs_cols]
# 50 observations , 42 Day 0

# Filter gwas data for dge/study samples --------------------------------------

# write keep list
keep <- geno_ids %>% dplyr::select(genotyping_id) %>% dplyr::rename(IID = genotyping_id)

# add in leading zero ids
keep$IID <- str_pad(keep$IID, 3, pad = "0")

fam <- fread(plink)
fam <- fam %>% dplyr::rename(FID = V1,
                      IID = V2)
# Get FIDs
keep <- inner_join(keep, fam, by = "IID") %>% dplyr::select(FID, IID)

write.table(keep,
            file = paste0(out_dir,study,time_point,"_keep.txt"),
            sep = "\t",
            row.names = F,
            quote = F,
            col.names = T)


# Filter for topsnps ----------------------------------------------------------

top <- fread(assoc_data)

snp.keep <- top$SNP

write.table(snp.keep, file = paste0(out_dir,study,time_point,"_snp_keep.txt"),
            sep = "\t",
            quote = F,
            col.names = F,
            row.names = F)

# Prep meQTL formats ----------------------------------------------------------

# Prepare genotyping data/snp location data ------------------------------------

# Run script in terminal as source ~/GWAS_22/Enteric_GWAS/DEG-analysis/eQTL/2.geno_snp_prep.sh

# Gene location table --------------------------------------------------------------

# Update array ids to gene names
row.names(exprs) <- fData$ensembl_gene_id #from fData

# Download ensembl human gene data
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
# default build is 38

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes2 <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl)

# Only assoc gene loc data with genes that are in our expression data
filtered_genes <- genes2[genes2$ensembl_gene_id %in% row.names(exprs),]

# filter autosomes
table(filtered_genes$chromosome_name)
filtered_genes <- filtered_genes[(filtered_genes$chromosome_name %in% c(1:22)),] 

# Fill in NAs
filtered_genes[filtered_genes == ""]<- NA
filtered_genes <- na.omit(filtered_genes) %>% dplyr::select(-hgnc_symbol) #rewrite table

write.table(filtered_genes, file = paste0(out_dir,study,time_point,"_gene_loc.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

# Expression table ------------------------------------------------------------

# Update Sample IDs to match geno IDs
pData <- left_join(pData, geno_ids, by = "part_number")

colnames(exprs) <- pData$genotyping_id
# add in leading zero ids
colnames(exprs) <- str_pad(colnames(exprs), 3, pad = "0")

write.table(exprs, file = paste0(out_dir,study,time_point,"_exprs.txt"),
            sep = "\t", quote = F, col.names = T, row.names = T)

# Covar table -----------------------------------------------------------------

covar <- fread("~/GWAS_22/gwas_final/meta/typhoid_pcacovar.txt")

covar <- dplyr::select(covar, IID, sex, age, chall_vax)
# filter for study participants
covar <- covar %>% filter(IID %in% keep$IID)

# recode chall_vax to numeric
covar$chall_vax <- as.factor(covar$chall_vax)
levels(covar$chall_vax)

covar$chall_vax <- ifelse(covar$chall_vax == "TxNone",0,
                        ifelse(covar$chall_vax == "TxTy21a",1,
                               ifelse(covar$chall_vax == "TxM01ZH09", 2,NA)))
# Transpose
covar <-  t(covar)
colnames(covar) <- covar[1, ]
covar <- covar[-1, ]


write.table(covar, file = paste0(out_dir,study,time_point,"_covar.txt"),
            sep = "\t", quote = F, col.names = T, row.names = T)


# Check orders -----------------------------------------------------------------
geno <- fread(paste0(out_dir,study,time_point,"_geno.txt"))
covar <- fread(paste0(out_dir,study,time_point,"_covar.txt"))
exprs <- fread(paste0(out_dir,study,time_point,"_exprs.txt"))

# Clean geno iid names
colnames(geno) <- str_sub(colnames(geno), start= -3)
geno <- geno %>% rename(V1 = SNP)

# Order other tables by geno
names.use <- names(geno)
covar <- covar[, ..names.use]
exprs <- exprs[, ..names.use]

# Final tables 
write.table(geno, file = paste0(out_dir,study,time_point,"_geno.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

write.table(covar, file = paste0(out_dir,study,time_point,"_covar.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

write.table(exprs, file = paste0(out_dir,study,time_point,"_exprs.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)


