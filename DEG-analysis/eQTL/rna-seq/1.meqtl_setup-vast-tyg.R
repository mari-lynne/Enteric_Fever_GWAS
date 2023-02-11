# VAST-TyGER RNA-seq eQTL analysis

# Aims: ------------------------------------------------------------------------

# Prepare RNAseq data previously preprocessed/batch corrected - counts transformed to TPM
# Alongside preparing genotyping data, also already processed in Plink
# Need to get both data sets to have the same participants, gene expression, snp data
# Also eQTL requires gene/snp location data - download from bioMart
# Then format tables as according to http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html#own

# Packages
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

# Get ensmembl data
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl) %>% rename(gene_id = ensembl_gene_id)

# study data and time point of interest
study <- c("/vast_tyg")
time_point <- c("/D0")
time <- c("D0")
var <- c("_gwas_")
deg <- data$top_tables$D7_D0 # optional associate with DEG genes, reduce MT
pathway <- c("~/GWAS_22/gwas_final/InnateDB_genes.csv") # optional, just assoc immune related genes

id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids.csv")
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/QC/typhoid2.IBD.fam")  # Fam file has participants linked 
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus2/typhoid_tophits.txt") # SNPs of interest that are sig in GWAS model 
covar_data <- c("~/GWAS_22/gwas_final/meta/typhoid_pcacovar.txt")

out_dir <- paste0(getwd(), study, "/Nov")

'%!in%' <- function(x,y)!('%in%'(x,y))


# 1) Filter dge data for time point --------------------------------------------

# Run mEQTL using tpm data instead of raw counts
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

# Get participants of interest - in this case D7, un vaccinated 
idx <- (data$samples$time_point == "D7")
idx <- (data$samples$time_point == "Baseline" & data$samples$vax_status == "no_vax")
sub <- data[,idx]
dim(sub)
check <- sub$samples[stri_duplicated(sub$samples$lab_id),] # remove dup ids if there are any

# 2) Filter GWAS data for study samples and DGE data ---------------------------

# Match by intersecting lab ids
IDlink <- clean_names(read.csv(file = id_data)) %>% dplyr::rename(lab_id = lab_id)

# Genotyping data in DGE data
geno_ids <- IDlink %>%
  filter(lab_id %in% sub$samples$lab_id) %>% dplyr::rename(IID = cat_iid)

# Get fam genoIDs
fam <- fread(plink) %>% rename(FID = V1, IID = V2) %>% select(FID, IID)
geno_ids <- inner_join(geno_ids, fam, by = "IID")
keep <- inner_join(geno_ids, fam, by = c("FID","IID")) %>% select(FID, IID)

write.table(keep,
            file = paste0(out_dir,time_point,"_keep.txt"),
            sep = "\t", row.names = F, quote = F, col.names = T)

# 3) Filter gene expression data -----------------------------------------------

# Filter for genotyped samples
idx <- sub$samples$lab_id %in% geno_ids$lab_id
sub <- sub[,idx]
# Add these ids to the dge list
geno_ids$lab_id <- as.character(geno_ids$lab_id) 
sub$samples$geno_ids <- inner_join(sub$samples, geno_ids, by = "lab_id") %>% select(FID, IID, lab_id)

# Filter for DEG genes
deg2 <- clean_names(deg)

deg <- deg2 %>% dplyr::filter(adj_p_val <= 0.05 & log_fc > 0.25 | log_fc < -0.25) #update using abs(log_fc)

# 4) Filter for genes of interest/pathway -----------------------------------------

# Immune genes (for vaccination study) 
immune <- clean_names(read.csv(file = pathway)) %>% rename(gene_id = ensembl) 
# Filter for autosomes and immune dge samples
immune <- immune %>% filter(chrom_name %in% (as.character(c(1:22))) & immune$gene_id %in% deg$gene_id)
 
idx <- sub$genes$gene_id %in% immune$gene_id
sub <- sub[idx,]

# Make gene location and gene expression tables
rownames(sub$counts) <- sub$genes$gene_name
colnames(sub$counts) <- sub$samples$geno_ids$IID
exprs <- sub$counts

filtered_genes <- dplyr::select(sub$genes, gene_name, chromosome_name, chrom_start, chrom_end)

write.table(exprs, file = paste0(out_dir,time_point,"_exprs.txt"),
            sep = "\t", quote = F,row.names = T)

write.table(filtered_genes, file = paste0(out_dir,time_point,"_gene_loc.txt"),
            sep = "\t", quote = F, row.names = F)

# 5) filter genotyping data for SNPs of interest -------------------------------

top <- fread(assoc_data) # Top GWAS SNPs or could be SNPs from region of interest
snp.keep <- top$SNP

write.table(snp.keep, file = paste0(out_dir,time_point,"_snp_keep.txt"),
            sep = "\t", quote = F, col.names = F, row.names = F)

# 6) Make covariates table -----------------------------------------------------

covar <- fread(covar_data)
# filter for study participants
covar <- covar %>% filter(IID %in% keep$IID)

covar <- dplyr::select(covar, IID, sex, age, bmi, vaccine)
# sequence_pool, if just vast no batch correction#, vaccine)

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

# 7) Prep meQTL formats --------------------------------------------------------
 
# Prepare genotyping data/snp location data using plink
# Update variables at the start of geno_prep script then run as:
# Run script as source ~/GWAS_22/Enteric_GWAS/DEG-analysis/eQTL/2.geno_snp_prep.sh

#Todo automate variable carry over

## Check data order ------------------------------------------------------------

# Matrix eqtl needs tables to be in same ID order

geno <- fread(paste0(out_dir, time_point, "_geno.txt"))
covar <- fread(paste0(out_dir, time_point, "_covar.txt"))
exprs <- fread(paste0(out_dir, time_point, "_exprs.txt"))

# Clean geno iid  for later matching
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

# 7) Run matrix-eqTL -----------------------------------------------------------

library(MatrixEQTL)
# load prepared data
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

##  Set Thresholds -------------------------------------------------------------

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
# Slice format allows for easy data compression

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

# 8) eQTL Results: -------------------------------------------------------------

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');

trans <- (me$trans$eqtls)
cis <- (me$cis$eqtls)

## Merge GWAS assoc data with new results 

trans$eqtl <- rep("trans", nrow(trans))
cis$eqtl <- rep("cis", nrow(cis))

eqtls <- bind_rows(cis, trans) 

assoc <- fread(assoc_data)
eqtls <- eqtls %>%
  dplyr::rename(SNP = snps) %>%
  left_join(assoc, by = "SNP")

### Final results + significance 

eqtl_sig <- eqtls %>% filter(FDR <= 0.05)

write.csv(eqtl_sig,
          file = paste0(out_dir, time_point, var,"sig_eqtl_immune_dge.csv"),
          row.names = FALSE)

write.csv(eqtls,
          file = paste0(out_dir, time_point, var,"all_eqtl_immune_dge.csv"),
          row.names = FALSE) #nold_novax

