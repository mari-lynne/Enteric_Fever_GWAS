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
time_point <- c("D0.12h")
study <- c("/vast_tyg_")
out_dir <- getwd()

id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids.csv")
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/QC/typhoid2.IBD.fam") #fam file
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/typhoid_tophits.txt")
covar_data <- c("~/GWAS_22/gwas_final/meta/typhoid_pcacovar.txt")

pheno_exprs <- read.csv(file = "~/RNA/oct/combat_vast_tyger.csv")
pheno <- pheno_exprs[, !str_detect(colnames(pheno_exprs), "ENSG")]

'%!in%' <- function(x,y)!('%in%'(x,y))

# Filter for time point and genotyped samples ----------------------------------

# Add exprs colnames to pheno_exprs for matching tables by

# 1) Filter time point
pheno_exprs <-
  pheno_exprs %>%
  filter(time_point == "D0.12h")



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

# Prep meQTL formats ----------------------------------------------------------

# Prepare genotyping data/snp location data ------------------------------------
# Update vars in geno_prep script and run as
# Run script in terminal as
#source ~/GWAS_22/Enteric_GWAS/DEG-analysis/eQTL/2.geno_snp_prep.sh

# Expression table ------------------------------------------------------------

# Update Sample IDs to match geno IDs

pheno_exprs <- left_join(pheno_exprs, geno_ids, by = "lab_id")
pheno_exprs <- pheno_exprs[pheno_exprs$cat_iid %in% keep$IID,]
exprs <- pheno_exprs[, str_detect(colnames(pheno_exprs), "ENSG")]
exprs <- t(exprs)
colnames(exprs) <- pheno_exprs$cat_iid

# Gene location table --------------------------------------------------------------

# Download ensembl human gene data
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
# default build is 38

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes2 <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl)

# Only assoc gene loc data with genes that are in our expression data
filtered_genes <- genes2[genes2$ensembl_gene_id %in% rownames(exprs),]

exprs <- as.data.frame(exprs)
exprs$ensembl_gene_id <- rownames(exprs)

exprs <- left_join(filtered_genes, exprs, by = "ensembl_gene_id")
exprs <- exprs[,c(2,6:ncol(exprs))]

# filter autosomes
table(filtered_genes$chromosome_name)
filtered_genes <- filtered_genes[(filtered_genes$chromosome_name %in% c(1:22)),] 

# Fill in NAs
filtered_genes[filtered_genes == ""]<- NA
filtered_genes <- na.omit(filtered_genes) %>% dplyr::select(-ensembl_gene_id) #rewrite table

write.table(filtered_genes, file = paste0(out_dir,study,time_point,"_gene_loc.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

write.table(exprs, file = paste0(out_dir,study,time_point,"_exprs.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

# Covar table -----------------------------------------------------------------

covar <- fread("~/GWAS_22/gwas_final/meta/typhoid_pcacovar.txt")
# filter for study participants
covar <- covar %>% filter(IID %in% keep$IID)
# 93-120, VAST-3104, 8280 (Vi-PS?, but control in pheno file)
covar <- dplyr::select(covar, IID, sex, age, chall_vax)


# recode chall_vax to numeric
covar$chall_vax <- as.factor(covar$chall_vax)
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
