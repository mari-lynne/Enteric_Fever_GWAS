# Plots IgA GWAS:

# IgA chr17:7559652:G/A:1 and IgA titre, diagnosis, gene expression

# Data set up ------------------------------------------------------------------
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

setwd("~/GWAS_22/gwas_final/merge/typhoid/vast/iga")
# snp_pheno <- read.csv(file = "vast_snp_pheno_D0.csv")
geno <- fread("~/GWAS_22/gwas_final/eQTL/vast/D0/vast_D0.12h.traw")
pheno2 <- fread("~/GWAS_22/gwas_final/meta/vast_pheno.txt") # see vast_assoc_notes.R
covar <- fread("~/GWAS_22/gwas_final/meta/covar_vast.txt")
y_var <- c("Vi-IgA")
gwas <- fread("iga_tophits.txt")
id_data <- clean_names(fread("~/GWAS_22/gwas_final/meta/geno_ids.csv"))
pheno_exprs <- read.csv(file = "~/GWAS_22/gwas_final/eQTL/vast/vast_tcpm.csv")

pheno_exprs <- pheno_exprs %>% filter(time_point3 == "V0") %>%
  inner_join(id_data, by = "lab_id") %>%
  dplyr::rename(IID = cat_iid)

pheno <- dplyr::select(pheno_exprs, IID, lab_id, study_arm, time, sex, ancestry, diagnosis) %>% left_join(pheno2, by = "IID") %>% mutate(geno_id = str_c(FID,IID, sep = "_"))
qtls <- c("IGHA2", "TNFSF13", "TNFRSF13B", "TNFSF17", "TNFSF12")

snp_interest <- c("chr17:7559652:A:G","chr20:9768452:G:A","chr18:66611261:A:G","chr2:53941744:A:G")
# chr17 - cis TNFSF13, day 0
# chr 20 - trans IgA, day 0 (FDR sig at D1)
# Chr 18 - 

# write.csv(snp_pheno, file = "all_geno_pheno.csv",row.names = FALSE)

# Gene expression data ---------------------------------------------------------
exprs <- pheno_exprs[,- c(1:51)]
exprs <- as.data.frame(t(exprs))
colnames(exprs) <- pheno$geno_id
exprs$hgnc_symbol <- rownames(exprs)
exprs <- exprs[ , !is.na(colnames(exprs))]

# Immune genes (for vaccination study) 
immune <- read.csv("~/GWAS_22/gwas_final/immune_gene_info.csv")
# Filter for autsomes and immune dge samples
immune <- immune %>% inner_join(exprs, by = "hgnc_symbol")
exprs <- immune[,c(2,10:ncol(immune))] %>% filter(hgnc_symbol %in% qtls)


snp <- geno[,-c(1,3:6)]
snp <-  t(snp)
colnames(snp) <- snp[1, ]
snp <- snp[-1, ]
snp <- as.data.frame(snp)
snp <- tibble::rownames_to_column(snp, var = "geno_id")
colnames(snp) <- str_replace_all(colnames(snp), "\\:", "_")
snp_pheno <-
  left_join(pheno, snp, by = "geno_id")

# SNP of interest 1 is Chr17_7559652_A_G (iGA pub hit)
# Recode factors
cols <- snp_pheno %>% dplyr::select(starts_with("Chr")) %>% colnames()
snp_pheno[,cols] <- lapply(snp_pheno[,cols], factor)
colnames(snp_pheno) <- str_replace_all(colnames(snp_pheno), "chr", "Chr")

# Make gene expression data numeric
snp_pheno[,c(1253:ncol(snp_pheno))] <- lapply(snp_pheno[,c(1253:ncol(snp_pheno))], as.numeric)


# 1 Plot IgA titre -------------------------------------------------------------

compare <- list(c("0", "1"))
compare <- list(c("0", "1"),c("0", "2"))

#Chr17_7559652_A_G 7544686
snp_pheno %>% drop_na(Chr17_7544686_C_T) %>%
  ggplot(aes(x = Chr17_7544686_C_T, y = vi_ig_a_titre, fill = Chr17_7544686_C_T)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#B94A22", "#BD9432", "#93913E")) +
  geom_jitter(width = 0.1,
              size = 1.7,
              alpha = 0.5) +
  scale_colour_manual(values = c("#061c23")) +
  stat_compare_means(
    comparisons = compare,
    method = "wilcox.test",
    label = "p.signif",
    size = 5,
    p.adjust.method = "none"
  ) +
  labs(y = "Vi-IgA titre", title = "rs3803800", x = "Chr17_7559652_A_G") +
  theme_bw()  +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 11.5),
    axis.text.y = element_text(size = 11.5)
  ) + scale_y_continuous(trans='log10')


# 2 Plot diagnosis -------------------------------------------------------------

snp_pheno %>% drop_na(diagnosis.x,Chr17_7559652_A_G) %>%
  ggplot(aes(x=diagnosis.x, fill=Chr17_7559652_A_G)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("#B94A22", "#BD9432", "#93913E"))  +
  labs(y = "count", x = "Diagnosis", title = "rs3803800") +
  theme_bw()  +
  theme(
    legend.position = "right", legend.title = element_text("rs3803800"),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 11.5),
    axis.text.y = element_text(size = 11.5)
  )

# 3 Plot EQTLs -----------------------------------------------------------------
# Label P-val and OR -----------------------------------------------------------
gwas$SNP <- str_replace_all(gwas$SNP , "chr", "Chr")
gwas$SNP <- str_replace_all(gwas$SNP, ":", "_")

# Get summary data
lab <- gwas[SNP == "Chr17_7559652_A_G", BETA]
lab <- signif(lab, digits = 2)
lab_p <- gwas[SNP == "Chr17_7559652_A_G", P]
lab_p <- signif(lab_p, digits = 2)
data <- snp_pheno %>% filter(!is.na(Chr17_7559652_A_G))

lab2 <- gwas[SNP == "Chr20_9768452_G_A", BETA]
lab2 <- signif(lab2, digits = 2)
lab_p2 <- gwas[SNP == "Chr20_9768452_G_A", P]
lab_p2 <- signif(lab_p2, digits = 2)
data2 <- snp_pheno %>% filter(!is.na(Chr20_9768452_G_A))

compare <- list(c("0", "1"))

TNFS <-
  ggplot(snp_pheno, aes(x = Chr17_7559652_A_G, y = TNFSF12, fill = Chr17_7559652_A_G)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
  geom_jitter(width = 0.1,
              size = 1.7,
              alpha = 0.5) +
  scale_colour_manual(values = c("#061c23")) +
  stat_compare_means(
    comparisons = compare,
    method = "wilcox.test",
    label = "p.signif",
    size = 5,
    p.adjust.method = "none"
  ) +
  labs(y = "Gene Expression (TPM)", title = "TNFSF13 (APRIL)") +
  theme_bw()  +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 11.5),
    axis.text.y = element_text(size = 11.5)
  )

TNFS


compare <- list(c("0", "1"), c("0", "2"))

TNFS12 <-
  ggplot(data, aes(x = Chr17_7525744_G_A, y = TNFSF12, fill = Chr17_7525744_G_A)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
  geom_jitter(width = 0.1,
              size = 1.7,
              alpha = 0.5) +
  scale_colour_manual(values = c("#061c23")) +
  labs(y = "Gene Expression (TPM)", title = "TNFSF12 (TWEAK)") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 11.5),
    axis.text.y = element_text(size = 11.5)
  ) +
  stat_compare_means(
    comparisons = compare,
    method = "wilcox.test",
    label = "p.signif",
    size = 5,
    p.adjust.method = "none"
  )

TNFS12
library(patchwork)

(TNFS|TNFS12) + plot_annotation(tag_levels = "a", title = "24 hours post-vaccination")