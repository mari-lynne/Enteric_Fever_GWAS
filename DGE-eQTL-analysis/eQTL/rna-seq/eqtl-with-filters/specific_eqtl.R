# Directories and variables ---------------------------------------------------
# Biomart data
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl)

setwd("~/GWAS_22/gwas_final/eQTL")
study <- c("/vast")
time_point <- c("/V1")
time <- c("V1")
var <- c("_tnf")
out_dir <- paste0(getwd(), study, "/nov8")

# VAST data 
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/vast/vast.fam") #fam file
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/vast/iga/iga_tophits.txt")
covar_data <- c("~/GWAS_22/gwas_final/meta/covar_vast.txt")
id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids.csv")
pathway <- c("~/GWAS_22/gwas_final/InnateDB_genes.csv")

'%!in%' <- function(x,y)!('%in%'(x,y))

# 1) Filter time point ---------------------------------------------------------

# data <- vast (from combat.r)
# 
# tpm3 <- function(counts,len) {
#   x <- counts/len
#   return(t(t(x)*1e6/colSums(x)))
# }
# 
# gene_info <- data$genes[data$genes$gene_id,c(1,6)]
# names.use <- rownames(data$counts) # order
# gene_info <- gene_info[names.use, 2]
# 
# tpm <- tpm3(data$counts, gene_info)
# data$raw <- data$counts
# data$counts <- tpm

# data <- readRDS(file = "VAST_tpm.rds")
# data <- readRDS(file = "VAST_tpm_combat.rds")

data <- readRDS(file = "VAST_tpm_nocomb.rds")

# Start ------------------------------------------------------------------------

idx <- (data$samples$time_point3 == "V1")
sub <- data[,idx]

# 2) Filter gene expression data -----------------------------------------------
idx <- sub$genes$gene_name %in% c("TNFSF13", "TNFSF12")
sub <- sub[idx,]
dim(sub)

# 3) Filter geno data and dge data for intersecting samples --------------------

# Match by intersecting participant ids
IDlink <- clean_names(read.csv(file = id_data)) %>% dplyr::rename(lab_id = lab_id)
geno_ids <- filter(IDlink, !is.na(genotyping_id)) %>% dplyr::rename(IID = cat_iid)

# Get fam IDs
fam <- fread(plink) %>% rename(FID = V1, IID = V2) %>% select(FID, IID)
geno_ids <- inner_join(geno_ids, fam, by = "IID")
keep <- left_join(fam,geno_ids, by = c("FID","IID")) %>% select(FID, IID)

write.table(keep, file = paste0(out_dir,time_point,var,"_keep.txt"),sep = "\t", row.names = F, quote = F, col.names = T)
# Add these ids to the dge list

# Filter DGE data for genotyped samples
idx <- sub$samples$lab_id %in% geno_ids$lab_id
sub <- sub[,idx]

# add ID info to samples
sub$samples <- inner_join(sub$samples, geno_ids, by = c("lab_id", "participant_id"))

# Make tables expression and gene_loc tables -----------------------------------
rownames(sub$counts) <- sub$genes$gene_name
colnames(sub$counts) <- geno_ids$IID
exprs <- sub$counts
target_genes <- dplyr::select(sub$genes, gene_name, chromosome_name, start_position, end_position)

write.table(exprs, file = paste0(out_dir,time_point,var,"_exprs.txt"),
            sep = "\t", quote = F,row.names = T)
write.table(target_genes, file = paste0(out_dir,time_point,var,"_gene_loc.txt"),
            sep = "\t", quote = F, row.names = F)

# 4) Filter association data --------------------------------------------------
top <- fread(assoc_data)

# risk loci
top <- filter(top, CHR == "17")

snp.keep <- top$SNP
write.table(snp.keep, file = paste0(out_dir,time_point,var,"_snp_keep.txt"),
            sep = "\t",
            quote = F,
            col.names = F,
            row.names = F)

# 6) Set up covar table -----------------------------------------------------------------

covar <- data.frame(sub$samples, stringsAsFactors=TRUE)
levels(covar$vaccine) <- c("3", "4")
covar$sex <- as.factor(covar$sex)
levels(covar$sex) <- c("2","1")
covar[covar$height_ == 0,24] <- 1.7 #dummy average height
covar <- covar %>% mutate(bmi =(weight_/(height_)^2))
covar <- dplyr::select(covar, IID, sex, age_at_do, vaccine, sequence_pool)
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
pvOutputThreshold_cis = 1e-1;
pvOutputThreshold_tra = 1e-10;

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

cis <- (me$cis$eqtls)
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

write.csv(eqtls,
          file = paste0(out_dir, time_point, var,"all_eqtl_cis.csv"),
          row.names = FALSE)
