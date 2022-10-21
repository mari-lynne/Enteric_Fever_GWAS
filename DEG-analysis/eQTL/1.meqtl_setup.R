#Library 
library(dplyr)
library(data.table)
library(tidyr)
library(limma)
library(stringr)
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

#Functions ####
'%!in%' <- function(x,y)!('%in%'(x,y))

setwd("/home/mari/GWAS_22/gwas_final/eQTL") #where data is saved
plot_dir <- c("/home/mari/GWAS_22/gwas_final/eQTL/plots") #where to save output plots

# Load micro-array experiment data
load("~/RNA/T1_T2_with_rsn.norm.R")
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
# Make new time point 0.5_1
pData$visit <- pData$visit %>% fct_collapse(D0_1 = c("D0+12", "D1"))

# Check T1 T2 IDs -------------------------------------------------------------

# Fix sample ids, T2 samples exprs names dont match
miss <- pData[sample_id %!in% colnames(exprs),]
miss2 <- colnames(exprs[ ,colnames(exprs) %!in% pData$sample_id])
miss$exprs_cols <- miss2
miss <- dplyr::select(miss, part_number, sample_id, exprs_cols, study_arm, visit)

# Add exprs colnames to pData for matching
pData$exprs_cols <- colnames(exprs)


# Filter DGE data for time point and genotyped samples ------------------------

# 1) Filter time point

pData <-
  pData %>%
  filter(visit == "TD")
exprs <-
  exprs[,colnames(exprs) %in% pData$exprs_cols]

# 2) Filter for gwas samples

# Match by intersecting participant ids
IDlink <- clean_names(read.csv(file = "~/GWAS_22/gwas_final/meta/T1T2_ID_link.csv"))
#Rename ID col so can join appropriately
IDlink <- dplyr::rename(IDlink, part_number = participant_id)

geno_ids <-
  IDlink %>%
  filter(part_number %in% pData$part_number)

pData <- 
  pData %>%
  filter(part_number %in% geno_ids$part_number)

exprs <-
  exprs[,colnames(exprs) %in% pData$exprs_cols]
# 50 observations  

# Filter gwas data for dge/study samples --------------------------------------

# write keep list
keep <- geno_ids %>% dplyr::select(genotyping_id) %>% dplyr::rename(IID = genotyping_id)
# add leading zeros
# add in leading zero ids
keep$IID <- str_pad(keep$IID, 3, pad = "0")

fam <- fread("~/GWAS_22/gwas_final/merge/typhoid/QC/typhoid2.IBD.fam")
fam <- fam %>% dplyr::rename(FID = V1,
                      IID = V2)
# Get FIDs
keep <- inner_join(keep, fam, by = "IID") %>% dplyr::select(FID, IID)
write.table(keep, file = "~/GWAS_22/gwas_final/meta/T1T2_keep.txt", row.names = F, quote = F, col.names = T)


# Filter for topsnps ----------------------------------------------------------

top <- fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/tophits_typhoid2.txt")

snp.keep <- top$SNP

write.table(snp.keep, file = "~/GWAS_22/gwas_final/eQTL/snp_keep.txt", sep = "\t", quote = F, col.names = F, row.names = F)

system("plink2 --bfile ~/GWAS_22/gwas_final/merge/typhoid/QC/typhoid2.IBD --keep ~/GWAS_22/gwas_final/meta/T1T2_keep.txt --extract ~/GWAS_22/gwas_final/eQTL/snp_keep.txt --make-bed --out ~/GWAS_22/gwas_final/eQTL/T1T2")

save.image(file = "meqtl_TD.RData")

# Prep meQTL formats ----------------------------------------------------------

# Prepare genotyping data/snp location data ------------------------------------

system("source ~/GWAS_22/Enteric_GWAS/DEG-analysis/eQTL/2.geno_snp_prep.sh")

# Gene loc table --------------------------------------------------------------

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

write.table(filtered_genes, file = "~/GWAS_22/gwas_final/eQTL/T1T2_gene_loc.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# Expression table ------------------------------------------------------------

# Update Sample IDs to match geno IDs
pData <- left_join(pData, geno_ids, by = "part_number")

colnames(exprs) <- pData$genotyping_id
#TODO fix zeros

write.table(exprs, file = "~/GWAS_22/gwas_final/eQTL/T1T2_exprs.txt", sep = "\t", quote = F, col.names = T, row.names = T)

# Covar table -----------------------------------------------------------------

load(file = "~/RNA/oct/meqtl_TD.RData")

covar <- fread("~/GWAS_22/gwas_final/meta/typhoid_pcacovar.txt")
covar <- select(covar, IID, sex, age, chall_vax)
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
write.table(covar, file = "~/GWAS_22/gwas_final/eQTL/T1T2_covar_eqtl.txt", sep = "\t", quote = F, col.names = T, row.names = T)


# Check orders -----------------------------------------------------------------
geno <- fread("T1T2_geno_file.txt")
covar <- fread("T1T2_covar_eqtl.txt",header = T)
exprs <- fread("T1T2_exprs.txt")


# Clean geno names
# geno needs iid names
colnames(geno) <- str_sub(colnames(geno), start= -3)
geno <- geno %>% rename(V1 = SNP)

# add in leading zero ids
colnames(exprs) <- str_pad(colnames(exprs), 3, pad = "0")
exprs <- exprs %>% rename(V1 = "0V1")

# Order other tables by geno
names.use <- names(geno)
exprs <- exprs[, ..names.use]
covar <- covar[, ..names.use]

write.table(geno, file = "T1T2_geno_file.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(covar, file = "T1T2_covar_eqtl.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(exprs, file = "T1T2_exprs.txt", sep = "\t", quote = F, col.names = T, row.names = F)


