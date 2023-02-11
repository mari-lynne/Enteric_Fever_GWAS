# Packages

#Library 
library(dplyr)
library(tidyr)
library(limma)
library(stringr)
library(data.table)
library(stringi)
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

#### Files -----------

# Typhoid:
# id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids.csv")
# plink <- c("~/GWAS_22/gwas_final/merge/typhoid/QC/typhoid2.IBD.fam")
# assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/typhoid_tophits.txt")
# covar_data <- c("~/GWAS_22/gwas_final/meta/typhoid_pcacovar.txt")
# Merged RNAseq data 
# study <- c("/vast_tyg")
# pheno_exprs <- read.csv(file = "~/RNA/oct/combat_vast_tyger.csv")
# Challenge data microarray:

# Just VAST:
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/vast/vast.fam") #fam file
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/vast/iga_tophits.txt")
covar_data <- c("~/GWAS_22/gwas_final/meta/covar_vast.txt")
pheno_exprs <- read.csv(file = "vast/vast_tcpm.csv")
# study_arm %in% c("ViTCV", "ViPS")
pathway <- c("vast/iga_genes.csv") 

# Directories and variables ---------------------------------------------------
setwd("~/GWAS_22/gwas_final/eQTL")

study <- c("/vast_tyg")
time_point <- c("/D7")
time <- c("D7")
var <- c("_gwas_")
#cgas <- c("FCGR")
out_dir <- paste0(getwd(), study)
'%!in%' <- function(x,y)!('%in%'(x,y))

deg_data <- c("V7_V0.csv")
pathway <- c("~/GWAS_22/gwas_final/InnateDB_genes.csv")

id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids.csv")
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/QC/typhoid2.IBD.fam")
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus2/typhoid_tophits.txt")
covar_data <- c("~/GWAS_22/gwas_final/meta/typhoid_pcacovar.txt")

pheno_exprs <- read.csv(file = "~/RNA/oct/combat_vast_tyger_d7_tpm.csv")

# 1) Filter dge data for time point --------------------------------------------

# Make new baseline time point, V0 for vaccinees and D0 for challenge
# pheno_exprs <- pheno_exprs %>%
#   mutate(time_point = ifelse(time_point3 == "D0" & (study_arm %in% c("ViTCV", "ViPS")), "V28", time_point3))
# pheno_exprs <- pheno_exprs %>%
#   mutate(time_point = ifelse(time_point %in% c("V0","D0"), "Baseline", time_point)) 

pheno_exprs <- pheno_exprs %>%
  dplyr::filter(time_point == "D7") 

table(pheno_exprs$study_arm, pheno_exprs$time_point)
pheno_exprs <- pheno_exprs[!stri_duplicated(pheno_exprs$lab_id),]


# 2) Filter gwas data for study samples and dge data ---------------------------

# Match by intersecting participant ids
IDlink <- clean_names(read.csv(file = id_data)) %>% dplyr::rename(lab_id = lab_id)

# Genotyping data in DGE data
geno_ids <- IDlink %>%
  filter(lab_id %in% pheno_exprs$lab_id) %>% dplyr::rename(IID = cat_iid)

# Get fam IDs
fam <- fread(plink) %>% dplyr::rename(FID = V1, IID = V2)

geno_ids <- inner_join(geno_ids, fam, by = "IID") %>% dplyr::select(-starts_with("V"))
keep <- inner_join(geno_ids, fam, by = c("FID","IID")) %>% dplyr::select(FID, IID)

write.table(keep,
            file = paste0(out_dir,time_point,"_keep.txt"),
            sep = "\t",
            row.names = F,
            quote = F,
            col.names = T)

# 3) Filter gene expression data -----------------------------------------------

# Filter for genotyped samples
pheno_exprs <- inner_join(geno_ids, pheno_exprs, by = "lab_id")

# Get expression data
exprs <- pheno_exprs[,str_detect(colnames(pheno_exprs),"ENS")]
# exprs <- pheno_exprs[,-c(1:51)]
exprs <- as.data.frame(t(exprs))
colnames(exprs) <- pheno_exprs$IID
# exprs$hgnc_symbol <- rownames(exprs)
exprs$ensembl_gene_id <- rownames(exprs)

# Filter for DEG genes
# deg <- clean_names(read.csv(file = deg_data))
# deg <- deg %>% dplyr::filter(p_value <= 0.05)
# exprs <- filter(exprs, ensembl_gene_id %in% deg$gene_id)

# Filter for genes of interest/pathway -----------------------------------------

# Biomart data
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl)

# Immune genes (for vaccination study) 
immune <- clean_names(read.csv(file = pathway)) %>% dplyr::rename(ensembl_gene_id = ensembl) %>% inner_join(genes, by = "ensembl_gene_id") # hgnc_symbol

# Filter for autsomes and immune dge samples
immune <- immune %>%
  filter(chromosome_name %in% (as.character(c(1:22)))) %>% inner_join(exprs, by = "ensembl_gene_id")

# Make tables
exprs <- immune[,c(27,31:ncol(immune))] # immune[,c(1,6:ncol(immune))]
filtered_genes <- dplyr::select(immune, hgnc_symbol, chromosome_name, start_position, end_position)

write.table(exprs, file = paste0(out_dir,time_point,"_exprs.txt"),
            sep = "\t", quote = F,row.names = F)
write.table(filtered_genes, file = paste0(out_dir,time_point,"_gene_loc.txt"),
            sep = "\t", quote = F, row.names = F)

# gwas_genes <- clean_names(read.csv(file = "~/GWAS_22/gwas_final/gwas_vaccine_traits2.csv"))
# gwas_genes <- rename(gwas_genes, hgnc_symbol = mapped_gene)
# gwas_genes <- inner_join(gwas_genes, genes2, by = "hgnc_symbol")
# gene_list <- c(immune$ensembl, gwas_genes$ensembl_gene_id)
# genes2 <- dplyr::filter(genes2, ensembl_gene_id %in% gene_list)
# Only assoc gene loc data with genes that are in our expression data
# filtered_genes <- genes2[genes2$ensembl_gene_id %in% rownames(exprs),]
# Fill in NAs
# filtered_genes[filtered_genes == ""]<- NA
# filtered_genes <- na.omit(filtered_genes) %>% dplyr::select(-ensembl_gene_id) #rewrite table
# Filter for gene of interest
# filtered_genes <-
  # filtered_genes[str_starts(filtered_genes$hgnc_symbol, cgas),]
# Filter na genes
# filtered_genes <- filter(filtered_genes, !is.na(hgnc_symbol)) 
# Filter expression data
# exprs <- filter(exprs, exprs$hgnc_symbol %in% filtered_genes$hgnc_symbol)


# Top SNPs ----------------------------------------------------------

top <- fread(assoc_data)
snp.keep <- top$SNP

write.table(snp.keep, file = paste0(out_dir,time_point,"_snp_keep.txt"),
            sep = "\t",
            quote = F,
            col.names = F,
            row.names = F)
# Covar table -----------------------------------------------------------------

covar <- fread(covar_data)
# filter for study participants
covar <- covar %>% filter(IID %in% keep$IID)
# 93-120, VAST-3104, 8280 (Vi-PS?, but control in pheno file)
# covar2 <- pheno_exprs[, !str_detect(colnames(pheno_exprs), "ENSG")]

# covar <- covar2 %>% dplyr::select(cat_iid, sequence_pool) %>%
# dplyr::rename(IID = cat_iid) %>% left_join(covar, by = "IID")

covar <- dplyr::select(covar, IID, sex, age, chall_vax) # chall_vax)
# sequence_pool, if just vast

# Recode chall_vax to numeric
table(covar$chall_vax)
covar$chall_vax <- ifelse(covar$chall_vax == "TxNone",0,
                          ifelse(covar$chall_vax == "TxTy21a",1,
                                 ifelse(covar$chall_vax == "TxM01ZH09", 2,
                                       ifelse(covar$chall_vax == "TxVi-PS", 3,
                                               ifelse(covar$chall_vax == "TxVi-TT", 4,
                                                     ifelse(covar$chall_vax == "TNxNone", 5,
                                        NA))))))

# Transpose
covar <-  t(covar)
colnames(covar) <- covar[1, ]
covar <- covar[-1, ]

write.table(covar, file = paste0(out_dir,time_point,"_covar.txt"),
            sep = "\t", quote = F, col.names = T, row.names = T)



# Prep meQTL formats ----------------------------------------------------------
# Prepare genotyping data/snp location data ------------------------------------

# Update vars in geno_prep script and run as:
# Run script as source ~/GWAS_22/Enteric_GWAS/DEG-analysis/eQTL/2.geno_snp_prep.sh

## Check orders -----------------------------------------------------------------
geno <- fread(paste0(out_dir, time_point, "_geno.txt"))
covar <- fread(paste0(out_dir, time_point, "_covar.txt"))
exprs <- fread(paste0(out_dir, time_point, "_exprs.txt"))


# Clean geno iid names
geno <- geno %>% dplyr::rename(V1 = SNP)
colnames(geno) <- str_sub(colnames(geno), start= -4)

exprs <- exprs %>% dplyr::rename(V1 = hgnc_symbol)
colnames(exprs) <- str_sub(colnames(exprs), start= -4)

colnames(covar) <- str_sub(colnames(covar), start= -4)

# Order other tables by geno
names.use <- names(geno)
covar <- covar[, ..names.use]
exprs <- exprs[, ..names.use]

# Final tables 
write.table(geno, file = paste0(out_dir,time_point,"_geno.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

write.table(covar, file = paste0(out_dir, time_point,"_covar.txt"),
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

me$cis$ntests
me$trans$ntests/2

0.05/436

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

## Final results + significance ------------------------------------------------
eqtl_sig <- eqtls %>% filter(FDR <= 0.05)

write.csv(eqtl_sig,
          file = paste0(out_dir, time_point, var,"sig_eqtl2.csv"),
          row.names = FALSE)
write.csv(eqtls,
          file = paste0(out_dir, time_point, var,"all_eqtl2.csv"),
          row.names = FALSE)

