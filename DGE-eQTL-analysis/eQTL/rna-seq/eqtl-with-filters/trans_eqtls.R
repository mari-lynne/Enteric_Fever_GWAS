# eQTL viewer
library(eQTpLot)

# Input
GWAS.df <- fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus2/typhoid_gwas.txt")
GWAS.df <- select(GWAS.df, CHR, BP, SNP, P, OR)

eQTL.df <- read.csv("~/GWAS_22/gwas_final/eQTL/vast_tyg/Nov/eqtls/D0_gwas_sig_eqtl_immune_dge_nold.csv")
eQTL.df <- select(eQTL.df, gene, SNP, pvalue, beta)  
eQTL.df$Tissue <- rep("whole_blood", length(eQTL.df))
eQTL.df <- rename(eQTL.df, gene_name = gene)
eQTL.df$Trait <- rep("Typhoid Fever", length(eQTL.df$gene_name))

# CHR    Start     Stop    Gene Build
Genes.df <- left_join(eQTL.df, genes, by = "gene_name") %>%
  select(chromosome_name, start_position, end_position, gene_name) %>%
  mutate(build = rep("hg38", length(eQTL.df$gene_name))) 
Genes.df$chromosome_name <- as.integer(Genes.df$chromosome_name)

# CHR_A     BP_A     SNP_A CHR_B     BP_B            SNP_B       R2
LD.df <- fread(file = "~/GWAS_22/gwas_final/eQTL/vast_tyg/Nov/h12.top.ld")
LD.df <- as.data.frame(LD.df)

eQTpLot(GWAS.df = GWAS.df.example, eQTL.df = eQTL.df.example, gene = "BBS1", 
        gbuild = "hg19",  trait = "LDL", tissue =  "all", CollapseMethod = "min")

# rename cols
names(GWAS.df) <- c("CHR", "BP", "SNP", "P", "BETA")
names(eQTL.df) <- c(names(eQTL.df.example), "Trait")
names(Genes.df) <- names(Genes.df.example)

# cis eqtls
eQTpLot(GWAS.df = GWAS.df, eQTL.df = eQTL.df, LD.df = LD.df, Genes.df = Genes.df,
        gene = "FCGR2C", gbuild = "hg38",  trait = "Typhoid Fever", tissue =  "whole_blood", sigpvalue_GWAS = 1e-4, R2min = 0.6, LDmin = 15)





# , 

eQTpLot(GWAS.df = GWAS.df, eQTL.df.example = eQTL.df, LD.df = LD.df,
        gene = "", gbuild = "hg38",  trait = "Typhoid Fever", tissue =  "whole_blood", sigpvalue_GWAS = 1e-4, R2min = 0.6, LDmin = 15)


str(LD.df.example)
str(GWAS.df.example)
str(LD.df)
# 