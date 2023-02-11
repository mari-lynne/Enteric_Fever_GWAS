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

# Biomart data
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl)

study <- c("/vast")
time_point <- c("/V1")
time <- c("V1")
var <- c("WARS_iga")
out_dir <- paste0(getwd(), study, "/nov8")
deg_data <- c("V1_V0.csv") 

# VAST data 
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/vast/vast.fam") #fam file
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/vast/iga/iga_tophits.txt")
covar_data <- c("~/GWAS_22/gwas_final/meta/covar_vast.txt")
pheno_exprs <- read.csv(file = "vast/vast_tcpm.csv")

id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids.csv")

'%!in%' <- function(x,y)!('%in%'(x,y))

# 1) Filter for time point --------------------------------------------------------
pheno_exprs <- pheno_exprs %>%
  dplyr::filter(time_point3 == "V1" & study_arm != "CTRL") 
check <- pheno_exprs[!stri_duplicated(pheno_exprs$lab_id),]

# 2) Filter gwas data for study samples and dge data ---------------------------

# Match by intersecting participant ids
IDlink <- clean_names(read.csv(file = id_data))

# Genotyping data in DGE data
geno_ids <- IDlink %>%
  filter(lab_id %in% pheno_exprs$lab_id) %>% dplyr::rename(IID = cat_iid)

# Get fam IDs
fam <- fread(plink) %>% dplyr::rename(FID = V1, IID = V2)

geno_ids <- inner_join(geno_ids, fam, by = "IID") %>% dplyr::select(-starts_with("V"))
keep <- inner_join(geno_ids, fam, by = c("FID","IID")) %>% dplyr::select(FID, IID) %>% filter(IID != "125-62")

write.table(keep,
            file = paste0(out_dir,time_point,var,"_keep.txt"),
            sep = "\t",
            row.names = F,
            quote = F,
            col.names = T)

# 3) Filter gene expression data -----------------------------------------------

# Filter for genotyped samples
pheno_exprs <- inner_join(geno_ids, pheno_exprs, by = "lab_id")

# Get expression data
exprs <- pheno_exprs[,-c(1:51)]
exprs <- as.data.frame(t(exprs))
colnames(exprs) <- pheno_exprs$IID
exprs$hgnc_symbol <- rownames(exprs)

# 4) Filter for genes of interest/pathway -----------------------------------------

# Immune genes (for vaccination study) 
immune <- read.csv("~/GWAS_22/gwas_final/immune_gene_info.csv")
# Filter for autosomes and immune dge samples
#immune <- immune %>% inner_join(exprs, by = "hgnc_symbol")
#immune <- genes %>% inner_join(exprs, by = "hgnc_symbol")

immune <- filter(immune, hgnc_symbol %in% top_tab$gene_name)
exprs <- immune[,c(2:ncol(immune))]
#exprs <- immune[,c(2,11:ncol(immune))]
immune <- filter(immune, chromosome_name %in% c(1:22))

# # Filter for DEG genes
# deg <- clean_names(read.csv(file = deg_data))
# deg <- deg %>% dplyr::filter(p_value <= 0.05)
# exprs <- filter(exprs, ensembl_gene_id %in% deg$gene_id)
# exprs <- exprs[,-1]

# Make tables
# immune <- immune[immune$hgnc_symbol %in% exprs$hgnc_symbol, ]

filtered_genes <- dplyr::select(immune, hgnc_symbol, chromosome_name, start_position, end_position)
#filtered_genes <- dplyr::select(immune, hgnc_symbol, chrom_name, chrom_start, chrom_end)

write.table(exprs, file = paste0(out_dir,time_point,var,"_exprs.txt"),
            sep = "\t", quote = F,row.names = F)
write.table(filtered_genes, file = paste0(out_dir,time_point,var,"_gene_loc.txt"),
            sep = "\t", quote = F, row.names = F)

# 5) Filter association data --------------------------------------------------
top <- fread(assoc_data)

# risk loci
# top <- filter(top, CHR == "17", BP %in% c(7440148:7571159))
top <- filter(top, CHR %in% immune$chromosome_name)

snp.keep <- top$SNP
write.table(snp.keep, file = paste0(out_dir,time_point,var,"_snp_keep.txt"),
            sep = "\t",
            quote = F,
            col.names = F,
            row.names = F)

# 6) Set up covar table -----------------------------------------------------------------

covar <- fread(covar_data)
# filter for study participants
covar <- covar %>% filter(IID %in% keep$IID)
# 93-120, VAST-3104, 8280 (Vi-PS?, but control in pheno file)
covar2 <- pheno_exprs[, c(1:51)]
covar <- covar2 %>% dplyr::select(IID, sequence_pool) %>% left_join(covar, by = "IID")

covar <- dplyr::select(covar, IID, sex, age, bmi, chall_vax, sequence_pool)

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
write.table(covar, file = paste0(out_dir,time_point,var,"_covar.txt"),
            sep = "\t", quote = F, col.names = T, row.names = T)

# 7) Run bash geno prep script

# 8) Reorder files ------------------------------------------------------------
geno <- fread(paste0(out_dir, time_point, var, "_geno.txt"))
covar <- fread(paste0(out_dir, time_point, var,"_covar.txt"))
exprs <- fread(paste0(out_dir, time_point, var,"_exprs.txt"))

# Clean geno iid names
geno <- geno %>% dplyr::rename(V1 = SNP)
colnames(geno) <- str_sub(colnames(geno), start= -4)
exprs <- exprs %>% dplyr::rename(V1 = hgnc_symbol)
colnames(exprs) <- str_sub(colnames(exprs), start= -4)
colnames(covar) <- str_sub(colnames(covar), start= -4)

intersect(names(geno), names(exprs))

# Order other tables by geno
names.use <- names(geno)
covar <- covar[, ..names.use]
exprs <- exprs[, ..names.use]



# Final tables 
write.table(geno, file = paste0(out_dir,time_point,var,"_geno.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

write.table(covar, file = paste0(out_dir, time_point,var,"_covar.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

write.table(exprs, file = paste0(out_dir,time_point,var,"_exprs.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

# Run Matrix EQTL --------------------------------------------------------------
library(MatrixEQTL)

SNP_file_name = paste0(out_dir, time_point,var, "_geno.txt");
snps_location_file_name = paste0(out_dir, time_point,var, "_snp_loc.txt");

# Gene expression file name
expression_file_name = paste0(out_dir, time_point,var, "_exprs.txt");
gene_location_file_name = paste0(out_dir, time_point, var,"_gene_loc.txt");

# Covariates file name
# Set to character() if there are no covariates
covariates_file_name = paste0(out_dir, time_point,var,"_covar.txt");

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
          file = paste0(out_dir, time_point, var,"sig_eqtl.csv"),
          row.names = FALSE)
write.csv(eqtls,
          file = paste0(out_dir, time_point, var,"all_eqtl.csv"),
          row.names = FALSE)

cis <- filter(eqtls, eqtl == "cis")
write.csv(eqtls,
          file = paste0(out_dir, time_point, var,"all_eqtl_cis.csv"),
          row.names = FALSE)

tnsf12 <-filter(cis, gene == "TNFSF12")
write.csv(tnsf12,
          file = paste0(out_dir, time_point, var,"tnsf12.csv"),
          row.names = FALSE)
