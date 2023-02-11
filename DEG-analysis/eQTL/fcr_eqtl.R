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


#### Directories and variables -------------------------------------------------
setwd("~/GWAS_22/gwas_final/eQTL")
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl) %>% rename(gene_id = ensembl_gene_id)

data <- readRDS(file = "~/RNA/oct/vast_tyger/vast_tyger.Rds")
data$raw <- data$counts
data$counts <- data$tpm
genes <- left_join(data$genes, genes, by = c("gene_id", "chromosome_name"))
data$genes$chrom_start <- genes$start_position
data$genes$chrom_end <- genes$end_position
# remove rows with no gene info
data$genes <- drop_na(data$genes)
idx <- rownames(data$counts) %in% data$genes$gene_id
data <- data[idx,]

study <- c("/vast_tyg")
time_point <- c("/h12")
time <- c("h12")
var <- c("_fcgr_")
deg <- data$top_tables$D7_D0
pathway <- c("~/GWAS_22/gwas_final/InnateDB_genes.csv")

id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids.csv")
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/QC/typhoid2.IBD.fam")
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/fcgr_assoc_tophits.txt")
covar_data <- c("~/GWAS_22/gwas_final/meta/typhoid_pcacovar.txt")
out_dir <- paste0(getwd(), study, "/Nov")
'%!in%' <- function(x,y)!('%in%'(x,y))


# 1) Filter dge data for time point --------------------------------------------
# idx <- (data$samples$time_point == "D7")
idx <- (data$samples$time_point == "D0.12h")
sub <- data[,idx]
dim(sub)
check <- sub$samples[stri_duplicated(sub$samples$lab_id),]

# 2) Filter GWAS data for study samples and DGE data ---------------------------

# Match by intersecting participant ids
IDlink <- clean_names(read.csv(file = id_data)) %>% dplyr::rename(lab_id = lab_id)

# Genotyping data in DGE data
geno_ids <- IDlink %>%
  filter(lab_id %in% sub$samples$lab_id) %>% dplyr::rename(IID = cat_iid)

# Get fam IDs
fam <- fread(plink) %>% rename(FID = V1, IID = V2) %>% select(FID, IID)
geno_ids <- inner_join(geno_ids, fam, by = "IID")
keep <- inner_join(geno_ids, fam, by = c("FID","IID")) %>% select(FID, IID)

write.table(keep,
            file = paste0(out_dir,time_point,"_keep.txt"),
            sep = "\t",
            row.names = F,
            quote = F,
            col.names = T)

# 3) Filter gene expression data -----------------------------------------------

# Filter for genotyped samples
idx <- sub$samples$lab_id %in% geno_ids$lab_id
sub <- sub[,idx]
# Add these ids to the dge list
geno_ids$lab_id <- as.character(geno_ids$lab_id) 
sub$samples$geno_ids <- inner_join(sub$samples, geno_ids, by = "lab_id") %>% select(FID, IID, lab_id)

# Filter for DEG genes
deg2 <- clean_names(deg)

deg <- deg2 %>% dplyr::filter(adj_p_val <= 0.05 & log_fc > 0.25 | log_fc < -0.25)
deg <- deg[str_detect(deg$gene_name, "FCG"),]

#deg <- deg2 %>% dplyr::filter(p_value <= 0.05 & log_fc > 0.1 | log_fc < -0.1)
#f2r <- clean_names(filter(sub$top_tables$D7_D0, gene_name == "F2R"))
#deg <- bind_rows(deg, f2r)

# 4) Filter for genes of interest/pathway -----------------------------------------

# Immune genes (for vaccination study) 
immune <- clean_names(read.csv(file = pathway)) %>% rename(gene_id = ensembl) 
# Filter for autosomes and immune dge samples
# immune <- immune %>% filter(chrom_name %in% (as.character(c(1:22))) & immune$gene_id %in% sub$genes$gene_id) 
immune <- immune %>% filter(chrom_name %in% (as.character(c(1:22))) & immune$gene_id %in% deg$gene_id)
 
idx <- sub$genes$gene_id %in% immune$gene_id
sub <- sub[idx,]

# Make tables
rownames(sub$counts) <- sub$genes$gene_name
colnames(sub$counts) <- sub$samples$geno_ids$IID
exprs <- sub$counts

filtered_genes <- dplyr::select(sub$genes, gene_name, chromosome_name, chrom_start, chrom_end)

write.table(exprs, file = paste0(out_dir,time_point,"_exprs.txt"),
            sep = "\t", quote = F,row.names = T)
write.table(filtered_genes, file = paste0(out_dir,time_point,"_gene_loc.txt"),
            sep = "\t", quote = F, row.names = F)

# 5) Top SNPs ----------------------------------------------------------

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

covar <- dplyr::select(covar, IID, sex, age, bmi, vaccine)
# sequence_pool, if just vast #, vaccine)

covar$vaccine <-
  ifelse(covar$vaccine == "Vi-PS", 3,
      ifelse(covar$vaccine == "Vi-TT", 4,
              ifelse(covar$vaccine == "None", 5,
                     NA)))

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
cisDist = 2e6; #1MB

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
          file = paste0(out_dir, time_point, var,"sig_eqtl_immune_dge_nold_novax.csv"),
          row.names = FALSE)

write.csv(eqtls,
          file = paste0(out_dir, time_point, var,"all_eqtl_immune_dge_nold_novax.csv"),
          row.names = FALSE)

