# eQTL viewer
library(eQTpLot)

setwd("~/GWAS_22/gwas_final/merge/typhoid/vast/iga")
# Input GWAS
GWAS.df <- fread("~/GWAS_22/gwas_final/merge/typhoid/vast/iga/iga_gwas.txt")
GWAS.df <- select(GWAS.df, CHR, BP, SNP, P, BETA)


# Input eQTL
# eQTL.df <- read.csv("~/GWAS_22/gwas_final/eQTL/vast/nov8/V1_igatnsf12.csv")
eQTL.df <- read.csv("~/GWAS_22/gwas_final/eQTL/vast/nov8/V1_tnfall_eqtl_cis.csv")
eQTL.df <- select(eQTL.df, gene, SNP, pvalue, beta)  
eQTL.df$Tissue <- rep("whole_blood", length(eQTL.df$gene))
eQTL.df <- rename(eQTL.df, gene_name = gene)
eQTL.df$Trait <- rep("Vi-IgA titre", length(eQTL.df$gene_name))
names(eQTL.df) <- c(names(eQTL.df.example), "Trait")

# Input LD
# make SNP keep file
extract <- select(eQTL.df, SNP.Id)
write.table(extract, file = "V1_eqtl13_snp_keep.txt", row.names = F, quote = F)
system("source ~/GWAS_22/Enteric_GWAS/Post-Impute/GWAS_analysis/5.computeLD.sh")
LD.df <- fread("V1_eqtl13.ld")

# Plot data

eQTpLot(GWAS.df = GWAS.df, eQTL.df = eQTL.df, LD.df = LD.df,
        gene = "TNFSF13", gbuild = "hg38",  trait = "Vi-IgA titre", tissue =  "whole_blood", sigpvalue_GWAS = 1e-4, R2min = 0.4, LDmin = 18, range = 140)

# Public eqtls

# region
info <- read.csv("~/GWAS_22/gwas_final/eQTL/vast/nov8/V1_igatnsf12.csv")

print(paste0(min(info$BP),":", max(info$BP)))
(max(info$BP) - min(info$BP))/1000

# Typhoid F2R ------------------------------------------------------------------

GWAS.df <- fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus2/typhoid_gwas.txt")
GWAS.df <- select(GWAS.df, CHR, BP, SNP, P, OR) %>% rename(BETA = OR)

eQTL.df <- read.csv("~/GWAS_22/gwas_final/eQTL/vast_tyg/Nov/D7_f2r_sig_eqtl_cis.csv")
eQTL.df <- select(eQTL.df, gene, SNP, pvalue, beta)  
eQTL.df$Tissue <- rep("whole_blood", length(eQTL.df$gene))
eQTL.df <- rename(eQTL.df, gene_name = gene)
eQTL.df$Trait <- rep("Typhoid Fever", length(eQTL.df$gene_name))
names(eQTL.df) <- c(names(eQTL.df.example), "Trait")

extract <- select(eQTL.df, SNP.Id)
write.table(extract, file = "~/GWAS_22/gwas_final/eQTL/vast_tyg/Nov/D7_f2r_snp_keep.txt", row.names = F, quote = F)
LD.df <- fread("~/GWAS_22/gwas_final/eQTL/vast_tyg/Nov/D7_f2r.ld")

eQTpLot(GWAS.df = GWAS.df, eQTL.df = eQTL.df,
        gene = "F2R", gbuild = "hg38",  trait = "Typhoid Fever", tissue =  "whole_blood", sigpvalue_GWAS = 1e-4, R2min = 0.4, LDmin = 5, range = 1150)

#  LD.df = LD.df
info <- read.csv("~/GWAS_22/gwas_final/eQTL/vast_tyg/Nov/D7_f2r_sig_eqtl_cis.csv")
print(paste0(min(info$BP),":", max(info$BP)))
(max(info$BP) - min(info$BP))/1000

