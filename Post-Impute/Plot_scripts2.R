# Genotype plotting

# Plot genotype count on x, with antibody response/expression/on y (linear associations)
# Plot genotype on x, with case control freq on y

# Run ~/GWAS_22/Enteric_GWAS/DEG-analysis/eQTL/2.geno_snp_prep.sh script to make traw file

# Packages ---------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidylog)
library(stringr)
library(stringi)
library(openxlsx)
library(data.table)
library(RColorBrewer)
library(forcats)
library(qqman)
library(janitor)
library(magrittr)

setwd("~/GWAS_22/gwas_final/merge/typhoid/assoc")

plot_dir <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/plots/")
pheno2 <- fread("~/GWAS_22/gwas_final/meta/pheno_typhoid.txt")
covar <- fread("~/GWAS_22/gwas_final/meta/covar_typhoid.txt")
id_data <- clean_names(fread("~/GWAS_22/gwas_final/meta/geno_ids.csv"))
y_var <- c("Expression (TPM)")
pathway <- c("~/GWAS_22/gwas_final/InnateDB_genes.csv")

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl)

gwas <- fread("nexus2/typhoid_tophits.txt")
geno <- fread("~/GWAS_22/gwas_final/eQTL/vast_tyg/D0.12h.traw")
pheno_exprs <- read.csv(file = "~/RNA/oct/combat_vast_tyger_tpm.csv")
pheno_exprs <- filter(pheno_exprs, time_point3 == "D0.12h") %>%
  inner_join(id_data, pheno_exprs, by = "lab_id") %>%
  dplyr::rename(IID = cat_iid)
pheno <- dplyr::select(pheno_exprs, IID, lab_id, study_arm, time, sex, ancestry, diagnosis) %>% left_join(pheno2, by = "IID") %>% mutate(geno_id = str_c(FID,IID, sep = "_"))
qtls <- c("F2R","MPL","ITGB3","MMRN1","GP9","GP1BA", "CD9", "THBS1", "ANGPT1","PTGS1")

# Gene expression data ---------------------------------------------------------
exprs <- pheno_exprs %>% dplyr::select(starts_with("ENS"))
exprs <- as.data.frame(t(exprs))
colnames(exprs) <- pheno$geno_id
exprs$ensembl_gene_id <- rownames(exprs)
exprs <- exprs[ ,!is.na(colnames(exprs))]

# Immune genes (for vaccination study) 
immune <- clean_names(read.csv(file = pathway)) %>% dplyr::rename(ensembl_gene_id = ensembl) %>% inner_join(genes, by = "ensembl_gene_id")
# Filter for autsomes and immune dge samples
immune <- immune %>%
  filter(chromosome_name %in% (as.character(c(1:22)))) %>% inner_join(exprs, by = "ensembl_gene_id")
exprs <- immune[,c(27,31:ncol(immune))] %>% filter(hgnc_symbol %in% qtls)

# Genotyping data --------------------------------------------------------------
snp <- geno[,-c(1,3:6)]
snp <-  t(snp)
colnames(snp) <- snp[1, ]
snp <- snp[-1, ]
snp <- as.data.frame(snp)
snp <- tibble::rownames_to_column(snp, var = "geno_id")
colnames(snp) <- str_replace_all(colnames(snp), "\\:", "_")

# Expression data
qtls <- as.data.frame(t(exprs))
names(qtls) <- qtls[1,]
qtls <- qtls[-1,]
qtls <- rownames_to_column(qtls, var = "geno_id")

# Join all tables
snp_pheno <-
  left_join(pheno, snp, by = "geno_id") %>%
  left_join(qtls, by = "geno_id")

# Recode factors
cols <- snp %>% dplyr::select(starts_with("chr")) %>% colnames()
snp_pheno[,cols] <- lapply(snp_pheno[,cols], factor)
snp_pheno$Diagnosed <- as.factor(snp_pheno$Diagnosed) 
snp_pheno$Diagnosed  <- fct_recode(snp_pheno$Diagnosed, nTD = "1", TD = "2")

# Make gene expression data numeric
colnames(snp_pheno) <- str_replace_all(colnames(snp_pheno), "chr", "Chr")
snp_pheno[,c(39:49)] <- lapply(snp_pheno[,c(39:49)], as.numeric)

# Set up labels ----------------------------------------------------------------

# Label P-val and OR
gwas$SNP <- str_replace_all(gwas$SNP , "chr", "Chr")
gwas$SNP <- str_replace_all(gwas$SNP, ":", "_")
# Get summary data
lab <- gwas[SNP == "Chr5_77682126_C_T", OR]
lab <- signif(lab, digits = 2)
lab_p <- gwas[SNP == "Chr5_77682126_C_T",P]
lab_p <- signif(lab_p, digits = 2)

# Get summary data
lab2 <- gwas[SNP == "Chr12_19099038_T_C", OR]
lab2 <- signif(lab2, digits = 2)
lab_p2 <- gwas[SNP == "Chr12_19099038_T_C",P]
lab_p2 <- signif(lab_p2, digits = 2)


# Plot D012h -------------------------------------------------------------------

data <- snp_pheno %>% filter(!is.na(Chr5_77682126_C_T))
compare <- list(c("0", "1"))
F2R <-
  ggplot(data, aes(x=Chr5_77682126_C_T, y=F2R, fill=Chr5_77682126_C_T)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
  geom_jitter(width = 0.1, size=1.7, alpha = 0.6) +
  scale_colour_manual(values = c("#061c23")) +
  stat_compare_means(comparisons = compare,
                     method = "wilcox.test",
                     label = "p.format", size = 4,
                     p.adjust.method = "none") +
  labs(y="Gene Expression (TPM)", title = "F2R") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size =11.5),
        axis.text.y = element_text(size =11.5)) +
  annotate("text",
           label= paste("OR:",lab,"\n","P:",lab_p, sep = " "),
           y=60, x=3, size=4, family="arial")

data2 <- snp_pheno %>% filter(!is.na(Chr12_19099038_T_C))
compare <- list(c("0", "2"))

gp9 <- 
  ggplot(data2, aes(x=Chr12_19099038_T_C, y=GP9, fill=Chr12_19099038_T_C)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
  geom_jitter(width = 0.1, size=1.7, alpha = 0.4) +
  scale_colour_manual(values = c("#061c23")) + ylim(0,300) +
  labs(y=" Gene Expression (TPM)", title = "GP9", x="Chr12_19094052_C_T") +
  theme_bw() +
  stat_compare_means(comparisons = compare,
                     method = "wilcox.test",
                     label = "p.signif", size = 5,
                     p.adjust.method = "none") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size =11.5),
        axis.text.y = element_text(size =11.5)) +
  annotate("text",
           label= paste("OR:",lab2,"\n","P:",lab_p2, sep = " "),
           y=220, x=3, size=4, family="arial")
gp9

MMRN1 <-
  ggplot(data2, aes(x=Chr12_19099038_T_C, y=MMRN1, fill=Chr12_19099038_T_C)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
  geom_jitter(width = 0.1, size=1.7, alpha = 0.4) +
  scale_colour_manual(values = c("#061c23")) +
  labs(y="Gene Expression (TPM)", title = "MMRN1", x= "Chr12_19094052_C_T") +
  theme_bw()  + ylim(0,16) +
  stat_compare_means(comparisons = compare,
                     method = "wilcox.test",
                     label = "p.signif",size = 5,
                     p.adjust.method = "none") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size =11.5),
        axis.text.y = element_text(size =11.5)) +
  annotate("text",
           label= paste("OR:",lab2,"\n","P:",lab_p2, sep = " "),
           y=12, x=3, size=4, family="arial")

write.csv(snp_pheno, file = "~/GWAS_22/gwas_final/eQTL/vast_tyg/VT_eqtl_12h.csv")

# Plot TD ----------------------------------------------------------------------

# Remake data and labs
geno <- fread("~/GWAS_22/gwas_final/eQTL/vast_tyg/TD.traw")
pheno_exprs <- read.csv(file = "~/RNA/oct/combat_vast_tyger_tpm.csv")
pheno_exprs <- filter(pheno_exprs, time_point3 == "TD") %>%
  inner_join(id_data, pheno_exprs, by = "lab_id") %>%
  dplyr::rename(IID = cat_iid)
pheno <- dplyr::select(pheno_exprs, IID, lab_id, study_arm, time, sex, ancestry, diagnosis) %>% left_join(pheno2, by = "IID") %>% mutate(geno_id = str_c(FID,IID, sep = "_"))
qtls <- c("F2R","MPL","ITGB3","MMRN1","GP9","GP1BA", "CD9", "THBS1", "ANGPT1","PTGS1")

# Gene expression data ---------------------------------------------------------
exprs <- pheno_exprs %>% dplyr::select(starts_with("ENS"))
exprs <- as.data.frame(t(exprs))
colnames(exprs) <- pheno$geno_id
exprs$ensembl_gene_id <- rownames(exprs)
exprs <- exprs[ ,!is.na(colnames(exprs))]

immune <- immune %>%
  filter(chromosome_name %in% (as.character(c(1:22)))) %>% inner_join(exprs, by = "ensembl_gene_id")
exprs <- immune[,c(27,31:ncol(immune))] %>% filter(hgnc_symbol %in% qtls)

# Genotyping data --------------------------------------------------------------
snp <- geno[,-c(1,3:6)]
snp <-  t(snp)
colnames(snp) <- snp[1, ]
snp <- snp[-1, ]
snp <- as.data.frame(snp)
snp <- tibble::rownames_to_column(snp, var = "geno_id")
colnames(snp) <- str_replace_all(colnames(snp), "\\:", "_")

# Expression data
qtls <- as.data.frame(t(exprs))
names(qtls) <- qtls[1,]
qtls <- qtls[-1,]
qtls <- rownames_to_column(qtls, var = "geno_id")

# Join all tables
snp_pheno <-
  left_join(pheno, snp, by = "geno_id") %>%
  left_join(qtls, by = "geno_id")

# Recode factors
cols <- snp %>% dplyr::select(starts_with("chr")) %>% colnames()
snp_pheno[,cols] <- lapply(snp_pheno[,cols], factor)
snp_pheno$Diagnosed <- as.factor(snp_pheno$Diagnosed) 

# Make gene expression data numeric
colnames(snp_pheno) <- str_replace_all(colnames(snp_pheno), "chr", "Chr")
snp_pheno[,c(47:56)] <- lapply(snp_pheno[,c(47:56)], as.numeric)

# Set up labels ----------------------------------------------------------------

# Label P-val and OR
gwas$SNP <- str_replace_all(gwas$SNP , "chr", "Chr")
gwas$SNP <- str_replace_all(gwas$SNP, ":", "_")
# Get summary data
lab <- gwas[SNP == "Chr5_77682126_C_T", OR]
lab <- signif(lab, digits = 2)
lab_p <- gwas[SNP == "Chr5_77682126_C_T",P]
lab_p <- signif(lab_p, digits = 2)
data <- snp_pheno %>% filter(!is.na(Chr5_77682126_C_T))

# TODO redo with the same snp (these are in high LD so results are the same and ok for now) 
lab2 <- gwas[SNP == "Chr12_19099038_T_C", OR]
lab2 <- signif(lab2, digits = 2)
lab_p2 <- gwas[SNP == "Chr12_19099038_T_C",P]
lab_p2 <- signif(lab_p2, digits = 2)
data2 <- snp_pheno %>% filter(!is.na(Chr12_19094052_C_T))

compare <- list(c("0","1"))
F2R_td <-
  ggplot(data, aes(x=Chr5_77682126_C_T, y=F2R, fill=Chr5_77682126_C_T)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
  geom_jitter(width = 0.1, size=1.7, alpha = 0.5) +
  scale_colour_manual(values = c("#061c23")) +
  stat_compare_means(comparisons = compare,
                     method = "wilcox.test",
                     label = "p.signif", size = 5,
                     p.adjust.method = "none") +
  labs(y="Gene Expression (TPM)", title = "F2R") +
  theme_bw()  +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size =11.5),
        axis.text.y = element_text(size =11.5)) +
  annotate("text",
           label= paste("OR:",lab,"\n","P:",lab_p, sep = " "),
           y=31, x=3, size=4, family="arial")

F2R_td

TD

compare <- list(c("0", "1"), c("0","2"))
gp9_td <- 
  ggplot(data2, aes(x=Chr12_19094052_C_T, y=GP9, fill=Chr12_19094052_C_T)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
  geom_jitter(width = 0.1, size=1.7, alpha = 0.4) +
  scale_colour_manual(values = c("#061c23")) +
  labs(y=" Gene Expression (TPM)", title = "GP9", x="Chr12_19094052_C_T") +
  theme_bw() +
  stat_compare_means(comparisons = compare,
                     method = "wilcox.test",
                     label = "p.signif", size = 5,
                     p.adjust.method = "none") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size =11.5),
        axis.text.y = element_text(size =11.5)) +
  annotate("text",
           label= paste("OR:",lab2,"\n","P:",lab_p2, sep = " "),
           y=150, x=3, size=4, family="arial")

gp9_td

compare <- list(c("0", "1"), c("0","2"))
MMRN1_td <-
  ggplot(data2, aes(x=Chr12_19094052_C_T, y=MMRN1, fill=Chr12_19094052_C_T)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
  geom_jitter(width = 0.1, size=1.7, alpha = 0.4) +
  scale_colour_manual(values = c("#061c23")) +
  labs(y="Gene Expression (TPM)", title = "MMRN1", x= "Chr12_19094052_C_T") +
  theme_bw() + ylim(0,16) +
  stat_compare_means(comparisons = compare,
                     method = "wilcox.test",
                     label = "p.signif",size = 5,
                     p.adjust.method = "none") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size =11.5),
        axis.text.y = element_text(size =11.5))
MMRN1_td

TD <- (F2R_td |gp9_td |MMRN1_td) + plot_annotation(tag_levels="a", title= "TD")





