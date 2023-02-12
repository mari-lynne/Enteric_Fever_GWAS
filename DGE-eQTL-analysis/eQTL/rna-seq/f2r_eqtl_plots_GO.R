library(tidyr)

setwd("~/GWAS_22/gwas_final/merge/typhoid/vast/iga")
plot_dir <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/plots/")
geno <- fread("~/GWAS_22/gwas_final/eQTL/vast/V1.traw")  # Make with geno_prep.sh script
pheno2 <- fread("~/GWAS_22/gwas_final/meta/vast_pheno.txt") # see vast_assoc_notes.R
covar <- fread("~/GWAS_22/gwas_final/meta/covar_vast.txt")
y_var <- c("Vi-IgA")
gwas <- fread("iga_tophits.txt")
id_data <- clean_names(fread("~/GWAS_22/gwas_final/meta/geno_ids.csv"))
pheno_exprs <- read.csv(file = "~/GWAS_22/gwas_final/eQTL/vast/vast_tcpm.csv")

pheno_exprs <- pheno_exprs %>% filter(time_point3 == "V1") %>%
  inner_join(id_data, by = "lab_id") %>%
  dplyr::rename(IID = cat_iid)

pheno <- dplyr::select(pheno_exprs, IID, lab_id, study_arm, time, sex, ancestry, diagnosis) %>% left_join(pheno2, by = "IID") %>% mutate(geno_id = str_c(FID,IID, sep = "_"))
qtls <- c("IGHA2", "TNFSF13", "TNFRSF13B", "TNFSF17", "TNFSF12")

snp_interest <- c("chr17:7559652:A:G","chr20:9768452:G:A","chr18:66611261:A:G","chr2:53941744:A:G")
# chr17 - cis TNFSF13, day 0
# chr 20 - trans IgA, day 0 (FDR sig at D1)
# Chr 18 - 


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
qtls <- tibble::rownames_to_column(qtls, var = "geno_id")

# Join all tables
snp_pheno <-
  left_join(pheno, snp, by = "geno_id") %>%
  left_join(qtls, by = "geno_id")

# Recode factors
cols <- snp %>% dplyr::select(starts_with("chr")) %>% colnames()
snp_pheno[,cols] <- lapply(snp_pheno[,cols], factor)

# Make gene expression data numeric
colnames(snp_pheno) <- str_replace_all(colnames(snp_pheno), "chr", "Chr")
snp_pheno[,c(298:ncol(snp_pheno))] <- lapply(snp_pheno[,c(298:ncol(snp_pheno))], as.numeric)

# write.csv(snp_pheno, file = "vast_snp_pheno_D0.csv",row.names = FALSE)
# write.csv(snp_pheno, file = "vast_snp_pheno_D1.csv",row.names = FALSE)
# write.csv(snp_pheno, file = "vast_snp_pheno_D7.csv",row.names = FALSE)

 # snp_pheno <- read.csv(file = "vast_snp_pheno_D1.csv")
 # D1 ---------------------------------------------------------------------------
 snp_pheno <- read.csv(file = "vast_snp_pheno_D1.csv")
 
# Recode factors
cols <- snp_pheno %>% dplyr::select(starts_with("Chr")) %>% colnames()
snp_pheno[,cols] <- lapply(snp_pheno[,cols], factor)

# Make gene expression data numeric
snp_pheno[,c(298:ncol(snp_pheno))] <- lapply(snp_pheno[,c(298:ncol(snp_pheno))], as.numeric)



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

# Plot EQTLs ------------------------------------------------------------------

compare <- list(c("0", "1"))

TNFS <-
  ggplot(data, aes(x = Chr17_7559652_A_G, y = TNFSF13, fill = Chr17_7559652_A_G)) +
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


# Day 7 ------------------------------------------------------------------------

snp_pheno <- read.csv(file = "vast_snp_pheno_D7.csv")

# Recode factors
cols <- snp_pheno %>% dplyr::select(starts_with("Chr")) %>% colnames()
snp_pheno[,cols] <- lapply(snp_pheno[,cols], factor)

# Make gene expression data numeric
snp_pheno[,c(298:ncol(snp_pheno))] <- lapply(snp_pheno[,c(298:ncol(snp_pheno))], as.numeric)


data <- snp_pheno %>% filter(!is.na(Chr17_7559652_A_G))
data2 <- snp_pheno %>% filter(!is.na(Chr20_9768452_G_A))

compare <- list(c("0","2"), c("0","1"))
TNFS12_7 <-
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

TNFS12_7

# IGHA2
compare <- list(c("0","2"))
IGHA2 <-
  ggplot(data2,
         aes(x =Chr20_9768452_G_A, y = IGHA2, fill = Chr20_9768452_G_A)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
  geom_jitter(width = 0.1,
              size = 1.7,
              alpha = 0.4) +
  scale_colour_manual(values = c("#061c23")) +
  labs(y = " Gene Expression (TPM)", title = "IGHA2") +
  theme_bw() +
  stat_compare_means(
    comparisons = compare,
    method = "wilcox.test",
    label = "p.signif",
    size = 5,
    p.adjust.method = "none"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 11.5),
    axis.text.y = element_text(size = 11.5)
  ) +scale_y_continuous(trans="log10")
IGHA2


# Luminex cors -----------------------------------------------------------------

cyt <- clean_names(read.csv(file = "~/ADNP/VAST_data/cytokine/38plex.csv"))
cyt <- filter(cyt, vaccine != "Control")
cyt28 <- filter(cyt, time_point == "PV")
cyt27 <- filter(cyt, time_point == "-27")
cyt25 <- filter(cyt, time_point == "-25")

iga_snp <- dplyr::select(snp_pheno,lab_id,Chr20_9768452_G_A)

cor <- left_join(cyt25, iga_snp, by = "lab_id")
cor <- cor %>% drop_na(Chr20_9768452_G_A)

colnames <- cor %>% dplyr::select(il_4,il_7,il_6,il_10,if_na2,if_ng)

X <- as.factor(cor$Chr20_9768452_G_A)

for(i in colnames){
  plt <- ggplot(cor, aes_string(x=X, y=i, fill = X)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = c("#A92D2F", "#CE7726","#9B9759")) +
    geom_jitter(width = 0.1,
                size = 1.7,
                alpha = 0.4) +
    scale_colour_manual(values = c("#061c23")) +
    theme_bw() + labs(x = "Chr20_9768452_G_A") +
    ylab(gsub("(\\_|\\-)", " ", i)) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 11.5),
      axis.text.y = element_text(size = 11.5)
    ) 
  print(plt)
  Sys.sleep(0.2)
}

# Cibersort correlations ------------------------------------------------------


# IgG tophits - intersect

gwas_g <- fread("~/GWAS_22/gwas_final/merge/typhoid/vast/igg/igg_2_tophits.txt")
gwas_a <- fread("~/GWAS_22/gwas_final/merge/typhoid/vast/iga/iga_tophits.txt")

common <- intersect(gwas_g$SNP, gwas_a$SNP)
gwas_c <- filter(gwas_a, SNP %in% common)
main_c <- filter(main, position %in% gwas_c$BP)


extra <- fread("~/GWAS_22/gwas_final/merge/typhoid/vast/igg/near_gens.txt")
down <- filter(extra, Position %in% gwas_c$BP)

genes <- unique(c(main$overlapped_gene, down$`Nearest Downstream Gene`))
genes <- str_subset(genes, "None", negate = TRUE)
genes <- genes[!is.na(genes)]
# filter lncRNA as not in enrichment set
genes <- str_subset(genes, "AC0|AL1|LIN|AC1|MIR44", negate = TRUE)


# Geno ontology clustering ----------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
load(file = "annot_data.RData")

# 1) Build GO map ---------------------------------------------------
x <- org.Hs.egSYMBOL2EG
# Get the entrez gene identifiers that are mapped to a gene symbol
hs <- org.Hs.eg.db
select(hs, 
       keys = my.symbols,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

gomap3 <- groupGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    level = 3,
    readable = TRUE
  )

View(gomap3@result)

# 2) Set up gene universe ------------------------------------------
library(biomaRt)
library(AnnotationDbi)
library(Organism.dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DOSE)

src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
chrom <- c(as.character(unique(gwas_c$CHR)))
chrom <- str_c(rep("chr", length(chrom)), chrom)

universe <- select(src, 
                   keys = chrom,
                   columns = c("symbol","entrez"),
                   keytype = "cds_chrom")

genes_ez <- select(src, 
                   keys = genes,
                   columns = c("symbol","entrez"),
                   keytype = "symbol")

genes_ez <- distinct(genes_ez)
# 3) run GO code ---------------------------------------------------------
ego <- enrichGO(gene = genes_ez$entrez,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", readable = TRUE)
View(ego@result)


edo <- enrichDO(gene = genes_ez$entrez,
                ont = "DO",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", readable = TRUE)

View(edo@result)
