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

# Directories and variables ---------------------------------------------------
setwd("~/GWAS_22/gwas_final/eQTL")
plot_dir <- c("~/GWAS_22/gwas_final/eQTL/plots") 
study <- c("/vast_")

time_point <- c("V0")
deg_data <- c("D0_V0.csv")
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/tophits_fcgr.txt")
var <- c("fcr")
cgas <- c("FCGR")

id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids.csv")
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/vast/vast.fam") #fam file
covar_data <- c("~/GWAS_22/gwas_final/meta/covar_vast.txt")
pheno_exprs <- read.csv(file = "vast_lcpm.csv")
out_dir <- getwd()
'%!in%' <- function(x,y)!('%in%'(x,y))

# Filter for time point and genotyped samples ----------------------------------
# 1) Filter time point
pheno_exprs <-
  pheno_exprs %>%
  filter(time_point3 == "V0", study_arm %in% c("ViTCV", "ViPS"))
table(pheno_exprs$time_point3)

# 2) Filter for gwas samples ---------------------------------------------------

# Match by intersecting participant ids

IDlink <- clean_names(read.csv(file = id_data))
#Rename ID col so can join appropriately
IDlink <- dplyr::rename(IDlink, lab_id = lab_id)

# Genotyping data in DGE data
geno_ids <-
  IDlink %>%
  filter(lab_id %in% pheno_exprs$lab_id) # Do missing samples have D1

# Filter DGE data for appropriate genotyped samples
pheno_exprs <- 
  pheno_exprs %>%
  filter(lab_id %in% geno_ids$lab_id)

# remove duplicates
test <- pheno_exprs[stri_duplicated(pheno_exprs$lab_id),]

# Filter gwas data for dge/study samples --------------------------------------

# write keep list
keep <- geno_ids %>% dplyr::select(cat_iid) %>% dplyr::rename(IID = cat_iid)

# add in leading zero ids
#keep$IID <- str_pad(keep$IID, 3, pad = "0")

fam <- fread(plink)
fam <- fam %>% dplyr::rename(FID = V1,
                             IID = V2)
# Get FIDs
keep_ids <- inner_join(keep, fam, by = "IID")
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

# Expression table ------------------------------------------------------------

# Update Sample IDs to match geno IDs
pheno_exprs <- left_join(pheno_exprs, geno_ids, by = "lab_id")

pheno_exprs <- pheno_exprs[pheno_exprs$cat_iid %in% keep$IID,]

exprs <- pheno_exprs[, str_detect(colnames(pheno_exprs), "ENSG")]
exprs <- t(exprs)

colnames(exprs) <- pheno_exprs$cat_iid


# Filter for DEG genes
# deg <- clean_names(read.csv(file = deg_data))
# deg <- deg %>% dplyr::filter(p_value <= 0.05)

exprs <- as.data.frame(exprs)
exprs$ensembl_gene_id <- rownames(exprs)

#exprs <- filter(exprs, ensembl_gene_id %in% deg$gene_id)

  # Gene location table --------------------------------------------------------------

# Download ensembl human gene data
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
# default build is 38

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes2 <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl)

# Only assoc gene loc data with genes that are in our expression data
filtered_genes <- genes2[genes2$ensembl_gene_id %in% rownames(exprs),]

exprs <- left_join(filtered_genes, exprs, by = "ensembl_gene_id")
exprs <- exprs[,c(2,6:ncol(exprs))]

# Filter autosomes
table(filtered_genes$chromosome_name)
filtered_genes <- filtered_genes[(filtered_genes$chromosome_name %in% c(1:22)),] 
# Fill in NAs
filtered_genes[filtered_genes == ""]<- NA
filtered_genes <- na.omit(filtered_genes) %>% dplyr::select(-ensembl_gene_id) #rewrite table
# Filter for gene of interest
filtered_genes <-
  filtered_genes[str_starts(filtered_genes$hgnc_symbol, cgas),]

# Filter na genes
filtered_genes <- filter(filtered_genes, !is.na(hgnc_symbol)) 

# Filter expression data
exprs <- filter(exprs, exprs$hgnc_symbol %in% filtered_genes$hgnc_symbol)


write.table(filtered_genes, file = paste0(out_dir,study,time_point,"_gene_loc.txt"),
            sep = "\t", quote = F, row.names = F)

write.table(exprs, file = paste0(out_dir,study,time_point,"_exprs.txt"),
            sep = "\t", quote = F,row.names = F)

# Covar table -----------------------------------------------------------------

covar <- fread(covar_data)
# filter for study participants
covar <- covar %>% filter(IID %in% keep$IID)
# 93-120, VAST-3104, 8280 (Vi-PS?, but control in pheno file)
covar2 <- pheno_exprs[, !str_detect(colnames(pheno_exprs), "ENSG")]

covar <- covar2 %>% dplyr::select(cat_iid, sequence_pool) %>%
  dplyr::rename(IID = cat_iid) %>% left_join(covar, by = "IID")

covar <- dplyr::select(covar, IID, sex, age, sequence_pool, chall_vax)

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

write.table(covar, file = paste0(out_dir,study,time_point,"_covar.txt"),
            sep = "\t", quote = F, col.names = T, row.names = T)



# Prep meQTL formats ----------------------------------------------------------
# Prepare genotyping data/snp location data ------------------------------------
# Update vars in geno_prep script and run as
# Run script as source ~/GWAS_22/Enteric_GWAS/DEG-analysis/eQTL/2.geno_snp_prep.sh

## Check orders -----------------------------------------------------------------
geno <- fread(paste0(out_dir,study,time_point,"_geno.txt"))
covar <- fread(paste0(out_dir,study,time_point,"_covar.txt"))
exprs <- fread(paste0(out_dir,study,time_point,"_exprs.txt"))


# Clean geno iid names
colnames(geno) <- str_sub(colnames(geno), start= -4)
colnames(covar) <- str_sub(colnames(covar), start= -4)
geno <- geno %>% dplyr::rename(V1 = SNP)
exprs <- exprs %>% dplyr::rename(V1 = hgnc_symbol)
colnames(exprs) <- str_sub(colnames(exprs), start= -4)

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

# Run eqTL ---------------------------------------------------------------------

library(MatrixEQTL)

out_dir = getwd()
SNP_file_name = paste0(out_dir, study, time_point, "_geno.txt");
snps_location_file_name = paste0(out_dir, study, time_point, "_snp_loc.txt");

# Gene expression file name
expression_file_name = paste0(out_dir, study, time_point, "_exprs.txt");
gene_location_file_name = paste0(out_dir, study, time_point, "_gene_loc.txt");

# Covariates file name
# Set to character() if there are no covariates
covariates_file_name = paste0(out_dir, study, time_point,"_covar.txt");

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

eqtl_sig <- eqtls %>% filter(FDR <= 0.05)

write.csv(eqtl_sig,
          file = paste0(out_dir, study, time_point, var,"sig_eqtl.csv"),
          row.names = FALSE)
write.csv(eqtls,
          file = paste0(out_dir, study, time_point, var,"all_eqtl.csv"),
          row.names = FALSE)

