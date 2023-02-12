# Libraries 
library(dplyr)
library(tidyr)
library(limma)
library(stringr)
library(tidylog)
library(data.table)
library(stringi)
library(ggplot2)
library(forcats)
library(factoextra)
library(ggpubr)
library(tibble)
library(stringr)
library(MatrixEQTL)
library(Biobase)
library(janitor)
library(biomaRt)

# load(file = "~/RNA/oct/meqtl_TD.RData")

#Functions ####


# Directories and variables ---------------------------------------------------

setwd("~/GWAS_22/gwas_final/eQTL")
plot_dir <- c("~/GWAS_22/gwas_final/eQTL/plots") 
study <- c("/T1T2")
time_point <- c("/Baseline")
time <- c("Baseline")
var <- c("_gwas_")
# cgas <- c("FCGR")
# deg_data <- c("D0_V0.csv")
# assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/fcgr_assoc_tophits.txt")

id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids_clean.csv")
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/QC/typhoid2.IBD.fam") #fam file
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus2/typhoid_tophits.txt")
covar_data <- c("~/GWAS_22/gwas_final/meta/typhoid_pcacovar.txt")
out_dir <- paste0(getwd(), study)

'%!in%' <- function(x,y)!('%in%'(x,y))

# Load micro-array experiment data
load("T1_T2_with_rsn.norm.R")
data = T1_T2_autosomes
rm(T1_T2_autosomes)

exprs <- as.data.frame(exprs(data))
pData <- clean_names(setDT(data@phenoData@data, keep.rownames = TRUE)[])
fData <- clean_names(data@featureData@data)

## Set up new vars -------------------------------------------------------------
# Make new typhoid diagnosis variable, na's = 0
pData <- pData %>%
  mutate(diagnosis = ifelse(is.na(pData$day_of_typhoid_diagnosis_from_challenge), 0, 1))
# Rename vars
names(pData)[names(pData) == "days_since_challenge"] <- "time"
names(pData)[names(pData) == "timepoint3"] <- "visit"
# Time point
table(pData$time, pData$study_arm)
pData$visit <- str_replace_all(pData$visit, "D0\\+12", "D0.12h")
# Make 12_24h time point
pData <- pData %>% mutate(visit = ifelse(visit %in% c("D0.12h","D1"), "12-24h", visit))
table(pData$visit, pData$study_arm)
# New baseline
pData <- pData %>% 
  mutate(visit = ifelse(visit == "V0", "V28", visit))
# Make new baseline time point using V0 for vaccinees and D0 for challenge
pData <- pData %>% 
  mutate(visit = ifelse(visit %in% c("V0","D0"), "Baseline", visit)) 
table(pData$visit, pData$study_arm)

# Filter DGE data for time point and genotyped samples ------------------------
# Add exprs colnames to pData for matching tables by
pData$exprs_cols <- colnames(exprs)

# 1) Filter time point
pData <- pData %>% filter(visit == "Baseline")
  #filter(time == "8" | time == "7" | time == "6.5")

# Remove duplicates
pData <- pData[!stri_duplicated(pData$part_number),]

# 2) Filter for GWAS samples

# Match by intersecting participant ids
IDlink <- clean_names(read.csv(file = id_data))
IDlink <- filter(IDlink, str_detect(meta_id, "T1|T2"))

# Participant IDs dont match T2
#Rename ID col so can join appropriately
IDlink <- dplyr::rename(IDlink, part_number = participant_id)

# Genotyping data in DGE data
geno_ids <- IDlink %>%
  filter(part_number %in% pData$part_number) # Do missing samples have D1

# Filter DGE data for appropriate genotyped samples
pData <- pData %>%
  filter(part_number %in% geno_ids$part_number)

exprs <- exprs[,colnames(exprs) %in% pData$exprs_cols]

# Filter gwas data for dge/study samples --------------------------------------

# Write keep list
keep <- geno_ids %>% dplyr::select(cat_iid) %>% dplyr::rename(IID = cat_iid)
# add in leading zero ids
keep$IID <- str_pad(keep$IID, 3, pad = "0")

fam <- fread(plink)
fam <- fam %>% dplyr::rename(FID = V1,
                      IID = V2)

# Get FIDs
keep <- inner_join(keep, fam, by = "IID") %>% dplyr::select(FID, IID)

write.table(keep,
            file = paste0(out_dir,time_point,"_keep.txt"),
            sep = "\t",
            row.names = F,
            quote = F)


# Filter for topsnps ----------------------------------------------------------

top <- fread(assoc_data)
snp.keep <- top$SNP

write.table(snp.keep, file = paste0(out_dir,time_point,"_snp_keep.txt"),
            sep = "\t",
            quote = F,
            col.names = F,
            row.names = F)

# Prep meQTL formats ----------------------------------------------------------

# Expression table ------------------------------------------------------------

# Update Sample IDs to match geno IDs
pData <- left_join(pData, geno_ids, by = "part_number")

# Add in leading zeros
pData$cat_iid <- str_pad(pData$cat_iid , 3, pad = "0")
colnames(exprs) <- pData$cat_iid
exprs$ensembl_gene_id <- rownames(exprs)

# Gene location table --------------------------------------------------------------

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl)

pathway <- c("~/GWAS_22/gwas_final/InnateDB_genes.csv")
immune <- clean_names(read.csv(file = pathway)) %>% dplyr::rename(ensembl_gene_id = ensembl) %>% inner_join(genes, by = "ensembl_gene_id") # hgnc_symbol

# Fill in NAs
immune[immune == ""]<- NA

# Filter for autsomes and immune dge samples
immune <- immune %>%
  filter(chromosome_name %in% (as.character(c(1:22))), !is.na(hgnc_symbol)) %>% inner_join(exprs, by = "ensembl_gene_id")

# Make tables
exprs <- immune[,c(27,31:ncol(immune))] # immune[,c(1,6:ncol(immune))]

filtered_genes <- dplyr::select(immune, hgnc_symbol, chromosome_name, start_position, end_position)

write.table(exprs, file = paste0(out_dir,time_point,"_exprs.txt"),
            sep = "\t", quote = F,row.names = F)
write.table(filtered_genes, file = paste0(out_dir,time_point,"_gene_loc.txt"),
            sep = "\t", quote = F, row.names = F)

# # Only assoc gene loc data with genes that are in our expression data
# filtered_genes <- genes[genes$ensembl_gene_id %in% rownames(exprs),]
# filtered_genes <- filtered_genes[(filtered_genes$chromosome_name %in% c(1:22)),] 
# # Fill in NAs
# filtered_genes[filtered_genes == ""]<- NA
# # filter na genes
# filtered_genes <- filter(filtered_genes, !is.na(hgnc_symbol))
# # Filter for gene of interest
# filtered_genes <-
#   filtered_genes[str_starts(filtered_genes$hgnc_symbol, cgas),]
# # Filter exprs
# exprs <- filter(exprs, row.names(exprs) %in% filtered_genes$ensembl_gene_id)
# # Update to symbol
# exprs$ensembl_gene_id <- rownames(exprs)
# exprs <- left_join(filtered_genes, exprs, by = "ensembl_gene_id")
# exprs <- exprs[,c(2,6:ncol(exprs))]
# filtered_genes <- filtered_genes %>% dplyr::select(-ensembl_gene_id)
# 
# write.table(filtered_genes, file = paste0(out_dir,time_point,"_gene_loc.txt"),
#             sep = "\t", quote = F, row.names = F)
# write.table(exprs, file = paste0(out_dir,time_point,"_exprs.txt"),
#             sep = "\t", quote = F, row.names = F)

# Covar table -----------------------------------------------------------------
covar <- fread("~/GWAS_22/gwas_final/meta/typhoid_pcacovar.txt")
covar <- dplyr::select(covar, IID, sex, age, dose, chall_vax)
# filter for study participants
covar <- covar %>% filter(IID %in% keep$IID)

# recode chall_vax to numeric
table(covar$chall_vax)

covar$chall_vax <- ifelse(covar$chall_vax == "TxNone",0,
                        ifelse(covar$chall_vax == "TxTy21a",1,
                               ifelse(covar$chall_vax == "TxM01ZH09", 2,NA)))

# Add extra covars
table(pData$array_experiment)
table(pData$time)

covar <- pData %>%
  dplyr::select(cat_iid, array_experiment) %>%
  dplyr::rename(IID = cat_iid) %>%
  left_join(covar) 

# covar <- dplyr::select(covar, -chall_vax)

# Transpose
covar <-  t(covar)
colnames(covar) <- covar[1, ]
covar <- covar[-1, ]

write.table(covar, file = paste0(out_dir,time_point,"_covar.txt"),
            sep = "\t", quote = F, col.names = T, row.names = T)

# Prepare genotyping data/snp location data ------------------------------------

# Run script in terminal as source ~/GWAS_22/Enteric_GWAS/DEG-analysis/eQTL/2.geno_snp_prep.sh

# Check orders -----------------------------------------------------------------
geno <- fread(paste0(out_dir,time_point,"_geno.txt"), header = T)
covar <- fread(paste0(out_dir,time_point,"_covar.txt"))
exprs <- fread(paste0(out_dir,time_point,"_exprs.txt"), header = T)

# Clean geno iid names
colnames(geno) <- str_sub(colnames(geno), start= -3)
geno <- geno %>% dplyr::rename(V1 = SNP)
exprs <- exprs %>% dplyr::rename(V1 = hgnc_symbol)

# Order other tables by geno
names.use <- names(geno)
covar <- covar[, ..names.use]
exprs <- exprs[, ..names.use]

# Final tables 
write.table(geno, file = paste0(out_dir,time_point,"_geno.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

write.table(covar, file = paste0(out_dir,time_point,"_covar.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

write.table(exprs, file = paste0(out_dir,time_point,"_exprs.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

# Run eqTL ---------------------------------------------------------------------

library(MatrixEQTL)

SNP_file_name = paste0(out_dir, time_point, "_geno.txt");
snps_location_file_name = paste0(out_dir, time_point, "_snp_loc.txt");

# Gene expression file name
expression_file_name = paste0(out_dir, time_point, "_exprs.txt");
gene_location_file_name = paste0(out_dir, time_point, "_gene_loc.txt");

# Covariates file name
# Set to character() if there are no covariates
covariates_file_name = paste0(out_dir, time_point,"_covar.txt");

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

## Thresholds -------------------------------------------------------------------

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 5e-2;
pvOutputThreshold_tra = 5e-4;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# Distance for local gene-SNP pairs
cisDist = 1e6; #1MB

# Load data --------------------------------------------------------------------

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # Space delim
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # Space delim
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # Space delim
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Load location files
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

# Run mEQTL analysis ------------------------------------------------------------

# Define model
useModel = modelLINEAR; 

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name      = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis  = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

# Results: --------------------------------------------------------------------

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');

trans <- (me$trans$eqtls)
cis <- (me$cis$eqtls)
bonferroni <- (0.05/me$trans$ntests)

#trans <- filter(trans, FDR <=0.05)
#cis <- filter(cis, FDR <=0.05)

## Merge gwas data ------------------------------------------------------------
trans$eqtl <- rep("trans", nrow(trans))
cis$eqtl <- rep("cis", nrow(cis))

eqtls <-
  bind_rows(cis, trans) 

## Merge assoc data ####
assoc <- fread(assoc_data)

eqtls <-
  eqtls %>%
  dplyr::rename(SNP = snps) %>%
  left_join(assoc, by = "SNP")
#%>%
# dplyr::select(-TEST)


## Save data -------------------------------------------------------------------
eqtl_sig <- eqtls %>% filter(FDR <= 0.05)

write.csv(eqtl_sig,
          file = paste0(out_dir, time_point, var,"sig_eqtl.csv"),
          row.names = FALSE)
write.csv(eqtls,
          file = paste0(out_dir, time_point, var,"all_eqtl.csv"),
          row.names = FALSE)



