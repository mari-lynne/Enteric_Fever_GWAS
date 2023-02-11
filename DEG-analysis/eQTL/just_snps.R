# geno, covar files stay the same, just change assoc

setwd("~/GWAS_22/gwas_final/eQTL")
plot_dir <- c("~/GWAS_22/gwas_final/eQTL/plots") 
time_point <- c("V7")
study <- c("/vast_")
out_dir <- getwd()

id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids.csv")
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/vast/vast.fam") #fam file
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/vast/igg_tophits.txt")
var <- c("igg")
covar_data <- c("~/GWAS_22/gwas_final/meta/covar_vast.txt")
# typhoid2.IBD, typhoid_pcacovar
# pheno_exprs <- read.csv(file = "~/RNA/oct/combat_vast_tyger.csv")
pheno_exprs <- read.csv(file = "vast_lcpm.csv")
deg_data <- c("D0_vs_V0.csv")

'%!in%' <- function(x,y)!('%in%'(x,y))


# Filter for topsnps ----------------------------------------------------------

top <- fread(assoc_data)

snp.keep <- top$SNP

write.table(snp.keep, file = paste0(out_dir,study,time_point,"_snp_keep.txt"),
            sep = "\t",
            quote = F,
            col.names = F,
            row.names = F)