# Matrix eQTL ------------------------------------------------------------------
library(MatrixEQTL)

# Update BLAS for faster matrix algrebra to atlas
# https://brettklamer.com/diversions/statistical/faster-blas-in-r/

# Set up file names and outputs ------------------------------------------------

setwd("~/GWAS_22/gwas_final/eQTL")
time_point <- c("0.5")
study <- c("/vast_tyg_")

assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/tophits_typhoid2.txt")
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
pvOutputThreshold_cis = 1e-2;
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

trans <- (me$trans$eqtls)
cis <- (me$cis$eqtls)
#trans <- filter(trans, FDR <=0.05)
#cis <- filter(cis, FDR <=0.05)

## Merge gene data ####

ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes <-
  getBM(
    attributes = c(
      'ensembl_gene_id',
      'hgnc_symbol',
      'chromosome_name',
      'start_position',
      'end_position'
    ),
    mart = ensembl
  )

genes <-
  genes %>%
  rename(gene = ensembl_gene_id)

eqtls <-
  bind_rows(cis, trans) %>%
  left_join(genes, by = "gene")

## Merge assoc data ####

assoc <- fread(assoc_data)

eqtls <-
  eqtls %>%
  rename(SNP = snps) %>%
  left_join(assoc, by = "SNP")

write.csv(eqtls,
          file = paste0(out_dir, study, time_point,"_eqtl.csv"),
          row.names = FALSE)

# Make the histogram of local and distant p-values
plot(me)
