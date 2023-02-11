# Directories and variables ---------------------------------------------------
# Biomart data
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','start_position','end_position'),mart = ensembl) %>% rename(gene_id = ensembl_gene_id)
# 
# data$genes <- left_join(data$genes, genes, by = c("gene_id"))
setwd("~/GWAS_22/gwas_final/eQTL/")

study <- c("/vast_tyg")
time_point <- c("/d7")
time <- c("d7")
var <- c("_TYG_")
pathway <- c("~/GWAS_22/gwas_final/InnateDB_genes.csv")

id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids.csv")
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/QC/typhoid2.IBD.fam")
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus2/typhoid_tophits.txt")
covar_data <- c("~/GWAS_22/gwas_final/meta/typhoid_pcacovar.txt")
out_dir <- paste0(getwd(), study, "/Nov/combat/redo")
'%!in%' <- function(x,y)!('%in%'(x,y))
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

# 1) Filter time point ---------------------------------------------------------
#  vast_tyger_combat.Rds
data <- readRDS(file = "combat_final_d7.RDS")

# just tyger
load(file = "~/RNA/oct/vt_bc_d7.RData")
data <- dge

# add gene info
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','start_position','end_position'),mart = ensembl) %>% rename(gene_id = ensembl_gene_id)
data$genes <- left_join(data$genes, genes, by = c("gene_id"))

# separate study bc
#data$counts <- data$combat1_counts

idx <- (data$samples$chall_vax == "T" & data$samples$time_point == "D7")
sub <- data[,idx]
dim(sub)

# TPM
tpm <- tpm3(sub$counts, sub$genes$transcript_length)
sub$counts <- tpm


# Start ------------------------------------------------------------------------
#& data$samples$chall_vax =="T"

idx <- (data$samples$time_point == "D7" & data$samples$chall_vax == "T")
sub <- data[,idx]
dim(sub)

# 2) Filter gene expression data -----------------------------------------------
# idx <- sub$genes$gene_name %in% c("TNFSF13", "TNFSF12")
idx <- sub$genes$gene_name %in% c("F2R", "F2RL1", "AGGF1", "TBCA", "TRIM69", "B2M", "ACE", "SORD")
sub <- sub[idx,]
dim(sub)

# 3) Filter geno data and dge data for intersecting samples --------------------

# Match by intersecting participant ids
IDlink <- clean_names(read.csv(file = id_data)) %>% dplyr::rename(lab_id = lab_id)

# Genotyping data in DGE data
geno_ids <- IDlink %>%
  filter(lab_id %in% sub$samples$lab_id) %>% dplyr::rename(IID = cat_iid)

# Get fam IDs
fam <- fread(plink) %>% rename(FID = V1, IID = V2) %>% select(FID, IID)
geno_ids <- inner_join(geno_ids, fam, by = "IID")
keep <- inner_join(geno_ids, fam, by = c("FID","IID")) %>% select(FID, IID)

write.table(keep, file = paste0(out_dir,time_point,var,"_keep.txt"),sep = "\t", row.names = F, quote = F, col.names = T)
# Add these ids to the dge list

# Filter DGE data for genotyped samples
idx <- sub$samples$lab_id %in% geno_ids$lab_id
sub <- sub[,idx]

sub$samples$lab_id <- as.integer(sub$samples$lab_id)
# add ID info to samples
sub$samples <- inner_join(sub$samples, geno_ids, by = c("lab_id", "participant_id"))

# Make tables expression and gene_loc tables -----------------------------------
rownames(sub$counts) <- sub$genes$gene_name
colnames(sub$counts) <- sub$samples$IID
dim(sub$counts)
exprs <- sub$counts
str(sub$genes)


target_genes <- dplyr::select(sub$genes, gene_name, chromosome_name, start_position, end_position)
#target_genes <- dplyr::select(sub$genes, gene_name, chromosome_name, chrom_start, chrom_end)

write.table(exprs, file = paste0(out_dir,time_point,var,"_exprs.txt"),
            sep = "\t", quote = F,row.names = T)
write.table(target_genes, file = paste0(out_dir,time_point,var,"_gene_loc.txt"),
            sep = "\t", quote = F, row.names = F)

# 4) Filter association data --------------------------------------------------
top <- fread(assoc_data)

# risk loci
# top <- filter(top, CHR == "5")

snp.keep <- top$SNP
write.table(snp.keep, file = paste0(out_dir,time_point,var,"_snp_keep.txt"),
            sep = "\t",
            quote = F,
            col.names = F,
            row.names = F)

# 6) Set up covar table -----------------------------------------------------------------
covar2 <- data.frame(sub$samples, stringsAsFactors=TRUE)
covar2 <- select(covar2, study_arm, IID, sequence_pool)
covar2$sequence_pool <- as.factor(covar2$sequence_pool)
levels(covar2$sequence_pool) <- c("1", "2")

covar <- fread(covar_data)
# filter for study participants
covar <- covar %>% filter(IID %in% keep$IID) %>% left_join(covar2, by = "IID")

#covar <- dplyr::select(covar, IID, sex, age, bmi, batch)
covar <- dplyr::select(covar, IID, sex, age, bmi, sequence_pool) 
# sequence_pool, if just vast #, vaccine)

covar$vaccine <-
  ifelse(covar$vaccine == "Vi-PS", 3,
         ifelse(covar$vaccine == "Vi-TT", 4,
                ifelse(covar$vaccine == "None", 5,
                       NA)))
# levels(covar$batch_pool) <- c(1:8) 

# Transpose
covar <-  t(covar)
colnames(covar) <- covar[1, ]
covar <- covar[-1, ]
write.table(covar, file = paste0(out_dir,time_point,var,"_covar.txt"),
            sep = "\t", quote = F, col.names = T, row.names = T)

# 7) Run bash geno prep script

# 8) Reorder files ------------------------------------------------------------
geno <- fread(paste0(out_dir, time_point, var, "_geno.txt"))
exprs <- fread(paste0(out_dir, time_point, var,"_exprs.txt"))
covar <- fread(paste0(out_dir, time_point, var,"_covar.txt"))

# Clean geno iid names
geno <- geno %>% dplyr::rename(V1 = SNP)
colnames(geno) <- str_sub(colnames(geno), start= -4)
colnames(exprs) <- str_sub(colnames(exprs), start= -4)
colnames(covar) <- str_sub(colnames(covar), start= -4)

# Order other tables by geno
names.use <- names(geno)
exprs <- exprs[, ..names.use]
covar <- covar[, ..names.use]

# Final tables 
write.table(geno, file = paste0(out_dir,time_point,var,"_geno.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

write.table(exprs, file = paste0(out_dir,time_point,var,"_exprs.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

write.table(covar, file = paste0(out_dir, time_point,var,"_covar.txt"),
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
output_file_name_tra = tempfile()
## Thresholds -------------------------------------------------------------------

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 5e-2;
pvOutputThreshold_tra = 1e-4;

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

cis <- (me$cis$eqtls)
trans <- (me$trans$eqtls)
bonferroni <- (0.05/me$cis$ntests)

me$cis$ntests
## Merge gwas data ------------------------------------------------------------

cis$eqtl <- rep("cis", nrow(cis))

eqtls <- cis

## Merge assoc data ####
assoc <- fread(assoc_data)

eqtls <-
  eqtls %>%
  dplyr::rename(SNP = snps) %>%
  left_join(assoc, by = "SNP")
sig_eqtl <- eqtls %>% filter(pvalue <= 0.05)
cis <- filter(eqtls, eqtl == "cis")

## Final results + significance ------------------------------------------------

# h12_vt <- sig_eqtl

write.csv(sig_eqtl,
          file = paste0(out_dir, time_point, var,"_eqtl_cis.csv"),
          row.names = FALSE)


