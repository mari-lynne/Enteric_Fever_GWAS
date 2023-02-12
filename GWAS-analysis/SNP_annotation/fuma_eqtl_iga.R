# Fuma output typhoid fever -----------------

library(stringr)
library(janitor)
library(data.table)
library(openxlsx)
library(readODS)
library(stringi)
library(AnnotationDbi)
library(dplyr)
library(tidylog)
library(Organism.dplyr)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DOSE)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(kableExtra)

setwd("~/GWAS_22/gwas_final/merge/typhoid/vast/fuma")

'%!in%' <- function(x,y)!('%in%'(x,y))

## EQTL fuma --------

data <- fread(file = "~/GWAS_22/gwas_final/merge/typhoid/assoc/vast/fuma/eqtl.txt")
# Get unique eqtls
DT <- data.table(data)
top_qtl <- DT[, .SD[which.min(p)], by = symbol]
#extra eqtls inc gtex
data <- fread(file = "~/GWAS_22/gwas_final/merge/typhoid/assoc/vast/fuma/FUMA_job216814/eqtl.txt")
# Get unique eqtls
DT <- data.table(data)
top_qtl2 <- DT[, .SD[which.min(p)], by = symbol]
top_qtl <- bind_rows(top_qtl, top_qtl2)
rm(top_qtl2)

# Mapped genes
data <- fread(file = "~/GWAS_22/gwas_final/merge/typhoid/assoc/vast/fuma/genes.txt")
# mapped genes
DT <- data.table(data)
top_qene <- DT[, .SD[which.min(minGwasP)], by = symbol]

## Magma -------
magma <- fread(file = "~/GWAS_22/gwas_final/merge/typhoid/assoc/vast/fuma/magma.genes.out")
magma <- filter(magma, P <0.05)
DT <- data.table(magma)
magma <- DT[, .SD[which.min(P)], by = GENE] # just keep one of each gene
magma <- clean_names(magma)

gene_list <- unique(c(magma$gene, top_qene$gene, top_qtl$gene))
snp_list <-unique(c(magma$symbol, top_qene$symbol, top_qtl$symbol))

# Top table --------------------------------------------------------------------
# Filter TD top table or TD nTD at day 7 (V7_V0)
# adj.P.Val <= 0.0, h12 - relax to 0.0 unadjust

load(file = "/home/mari/RNA/VAST_deg.RData")
sig_v7 <- V7_V0  %>% filter(adj.P.Val <= 0.05) #2667


top_tab <- V1_V0 %>% mutate(sig_pos = row_number()) %>% filter(adj.P.Val <= 0.05 & gene_id %in% gene_list & (logFC < -0.1 | logFC >0.1))

top_tab2 <- V7_V0 %>% mutate(sig_pos = row_number()) %>% filter(adj.P.Val <= 0.05 & gene_id %in% gene_list & (logFC < -0.1 | logFC >0.1))

top_tab <- bind_rows(top_tab, top_tab2)

top_tab <- top_tab %>% mutate(eqtl = ifelse(gene_id %in% top_qtl$gene, "y", "n"),
                              magma = ifelse(gene_id %in% magma$gene, "y", "n"),
                              inter = ifelse((eqtl == "y" & magma == "y"), "n", NA))

# Genes are sig lower in TD or higher in nTD
protec <- filter(top_tab, logFC < -0.1)
sus <- filter(top_tab, logFC > 0.1)

# STAT1(98% of iterations),SLAMF8(76%),PSME2(39%),WARS(37%) and ALDH1A1(36%; Fig 4A)
#write.csv(top_tab, file = "final_top.csv", row.names = F)

## Background -------------------------------------------------------------------
# check what wasn't detected

# keep intersecting genes only
gene_list <-  gene_list[gene_list %in% V7_V0$gene_id]

check1 <- top_qtl[top_qtl$symbol %in% V7_V0$gene_name] #31
check2 <- magma[magma$symbol %!in% V7_V0$gene_name]#36 - 310 = 346 of og genes not measured
miss <- bind_rows(check1, check2)
background_snps <- gene_list[gene_list %!in% miss$gene]
background_mag <- magma[magma$gene %in% background_snps]
background_qtls <- top_qtl[top_qtl$gene %in% background_snps]
background_genes <-  V7_V0 %>% filter(adj.P.Val <= 0.05 & (logFC < -0.1 | logFC >0.1 ))

# so 415 out of 759 annotated genes could be measured in top table - 251 of 415 snps in toptable
# top table 5722 diff expressed genes
# 7.5~ top genes were also in our annotated snps

spread <- list('final_top' = top_tab,
               'back_snps' = background_mag,
               'back_qtl' = background_qtls,
               'back_topgene' = background_genes,
               'top_magma' = magma,
               'top_qtls' = top_qene)

write.xlsx(spread, file = "snp_deg_all.xlsx", rowNames = F)

top_snp <- unique(magma$gene) # unique snps (magma has multiple gene pairs)

## Sig values ---------
top_tab <- rename(top_tab, gene = gene_id)
sig_qtl <- tidylog::inner_join(top_qtl, top_tab, by = "gene")
sig_mag <- tidylog::inner_join(magma, top_tab, by = "gene")

spread <- list('sig_qtl' = sig_qtl,
               'sig_mag' = sig_mag)
write.xlsx(spread, file = "snp_deg_sig.xlsx", rowNames = F)

sig <- sig_qtl %>%
  select(gene, uniqID, db, tissue, p, signed_stats, alignedDirection) %>%
  rename(qtl_p = p)
sig_m <- sig_mag %>% select(gene, nsnps, p) %>% rename(gene_p = p)

rank <- left_join(top_tab, sig, by = "gene")
rank <- left_join(rank, sig_m, by = "gene")

gsea <- rank %>% select(geneID, logFC, P.Value)
write.csv(gsea, file = "gsea.csv", row.names = F)

# get gwas summary statistics (diff builds)
gwas <- clean_names(fread("typhoid_37_gwas.txt"))
names(rank)
sum <- rank %>% select(-magma, -eqtl, -inter)


# VOlCANO VIEW -----------------------------------------------------------------

deg <- V1_V0
snp_genes <- unique(top_tab$gene_name)

plot_vol(deg,
              lab = 'top_lab',
              genes = snp_genes,
         title = "7-days post vaccination", top = 6)


qtls <- unique(top_qtl$symbol)
vol <- plot_vol(deg,
         lab = 'top_lab',
         genes = qtls,
         title = "24h post-vaccination", top = 6)

vol
save.image(file = "iga_fuma.RData")

# use v1 in thesis

deg <- V7_V0
qtls <- unique(top_qene$symbol)
plot_vol(deg,
         lab = 'sig',
         genes = qtls,
         title = "7-days post-vaccination")


# Gene ontology ----------------------------------------------------------------

# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html
# https://isglobal-brge.github.io/Master_Bioinformatics/enrichment-analysis.html

## Common genes  --------------
# from final top_tab snp_genes <- unique(top_tab$gene_name)

go_snpG <- groupGO(
  gene = snp_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  level = 3,
  readable = TRUE
)

View(go_snpG@result)

## Just SNP genes ---------------
go_snps <-
  groupGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  level = 3,
  readable = TRUE
)

## Just DGE genes --------------
ghe_genes <- V7_V0 %>% filter(adj.P.Val <= 0.05 & (logFC < -0.1 | logFC >0.1))
dge_genes <- unique(ghe_genes$gene_name)

go_dge <-
  groupGO(
    gene = dge_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    level = 3,
    readable = TRUE
  )

View(go_snpG@result)
View(go_dge@result)
View(go_snps@result)

# Enrichment start -------------------------------------------------------------------

## Get enterez ids -----------
src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")

snp_genes <- unique(top_tab$gene_name)
snp_genes <- sig_all

genes_ez <- select(src, 
                   keys = snp_genes,
                   columns = c("symbol","entrez"),
                   keytype = "symbol")

genes_ez <- distinct(genes_ez)

## ) Run ernrichGO code -------

ego <- enrichGO(gene = genes_ez$entrez,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", readable = TRUE)

View(ego@result)

snp_ego <- ego
dge_ego <- ego
cross_ego <-ego


View(snp_ego@result)

# could filter by significance then compare clusters
library(enrichplot)
barplot(dge_ego, showCategory=20)
dotplot(cross_ego, showCategory=20) + ggtitle("dotplot for GSEA")


## MsigDB -------

c7 <- read.gmt("~/Downloads/public_data/c7.immunesigdb.v2022.1.Hs.entrez.gmt")
c5 <- read.gmt("~/Downloads/public_data/c5.go.bp.v2022.1.Hs.entrez.gmt")

msig_ego <- enricher(genes_ez$entrez, TERM2GENE=c7) #c5
msig <- subset(msig_ego@result, p.adjust <0.05)

# update descriptions
msig_ego@result$Description <- str_replace_all(msig_ego@result$Description, "_", " ")
msig_ego@result$Description <- str_remove_all(msig_ego@result$Description, "^[GSE\\d]*")
msig_ego@result$Description <- str_to_title(msig_ego@result$Description)

msig_ego@result <- subset(msig_ego@result, str_detect(msig_ego@result$Description, "Newcastle", negate = TRUE))
msig_ego@result$Description <-
  str_replace_all(msig_ego@result$Description, "Tcell", "T-cell") %>%
  str_replace_all("Cd", "CD") %>%
  str_replace_all("Kaech|Dn", "") %>%
  str_replace_all("Vs", "vs") %>%
  str_replace_all("Gammadelta", "yD") %>%
  str_replace_all("A Phagocytophilum", "Phagocyte") %>%
  str_replace_all("Lps", "LPS") %>%
  str_replace_all("Pma", "PMA") %>%
  str_replace_all("Pbmc", "PBMC") %>%
  str_replace_all("Up", "(upreg)") %>%
  str_replace_all("Id2", "ID2") %>%
  str_replace_all("Ko", "KO") %>%
  str_replace_all("With", "-") %>%
  str_replace_all("Pdc", "PDc")  %>%
  str_replace_all("Foxp3", "FoxP3") %>%
  str_replace_all("Transduced", "transduced")

# filter weird rows
msig_ego@result <- filter(msig_ego@result, str_detect(Description, "Lupus", negate = TRUE))

imm <- dotplot(msig_ego, showCategory=15) + ggtitle("Immune GSEA: GWAS-DEG")

c5 <- dotplot(msig_ego, showCategory=15) + ggtitle("BP: GWAS-DEG")



# Biological function c5

# plot with volcano

vol|imm + plot_annotation(tag_levels = "a")


# save gene ontology genes ------------
write.csv(genes_ez, file = "cross_genes_ez.csv", row.names = F)

result <- msig_ego@result

neut_genes <- result[str_detect(result$Description, "Neutrophil"),]
neut_genes <- neut_genes[1:4,8]

neut_ez <- unlist(strsplit(neut_genes, "/"))

msig_ez <- unlist(strsplit(msig$geneID, "/"))

neut_ez <- select(src, 
                   keys = neut_ez,
                   columns = c("symbol","entrez"),
                   keytype = "entrez")

neut_ez <- distinct(neut_ez) %>% rename(gene_name = symbol)

all_genes <- final_top$gene_name
# gene ontology
go_all <- groupGO(
  gene = all_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "MF",
  level = 3,
  readable = TRUE
)

View(go_all@result)
write.csv(go_all@result, "typhoid-gene-ont-toptab-MF.csv", row.names = F)



msig_ez <- distinct(msig_ez) %>% rename(gene_name = symbol)

View(msig_ego@result)

# get joining snp data too
msig_rank <- rank %>% dplyr::select(-inter) %>% left_join(msig_ez, by = "gene_name")
write.csv(rank, file = "cross_rank_ez.csv", row.names = F)
ont_genes <- subset(rank, !is.na(entrez))
write.csv(ont_genes, file = "ont_genes.csv", row.names = F)

## link2 loci ----

# For eQTL analysis, test all GWAS sig snps in loci, with msig genes
# so we need a gene name list to test
# and then also filter top snps for snps in magma (but could just try without first)

loci <- fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/FUMA_extra2/GenomicRiskLoci.txt")

GenomicRiskLoci.txt

ont_genes <- left_join(msig_ez, sig_mag, by = "gene_name") %>% subset(!is.na(entrez))

sig_loci <- left_join(ont_genes, loci, by = "chr") 

sig_loci$chr<- as.factor(sig_loci$chr)
sig_loci <- subset(sig_loci, !is.na(GenomicLocus))
write.csv(sig_loci, file = "sig_loci.csv", row.names = F)

# Plots --------------------------------------------------------------

b <-
  sig_loci %>% group_by(chr) %>%
  ggplot(aes(x=chr, y =nIndSigSNPs, fill = chr)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE, option = "inferno") + theme_bw() +
  labs(y="Ind Sig SNPs (n)", x = "chromosome locus") + theme(legend.position = "none")


c <- dotplot(msig_ego, showCategory=15) + ggtitle("Immune GSEA: GWAS-DGE genes") +
  theme(axis.text.y=element_text(size=11))
c

(a|b) # add volcano

deg <- V1_V0
loci_genes <- unique(sig_loci$gene_name)

vol <- plot_vol(deg,
         lab = 'top_lab',
         genes = loci_genes,
         p = 0.01,
         title = "7-days post challenge (TD vs nTD)", top = 6,
         colours = c("#dd513a", "#932667", "#420a68"))

vol/(b|c) + plot_annotation(tag_levels = "a")

write.csv(V7_V0, file = "deg_V7_V0.csv", row.names = F)

# RNAi screen ----------------------------------------
load(file = "iga_fuma.RData")

si_up <- clean_names(read.xlsx("~/Downloads/public_data/siRNAup.xlsx", sheet = 1))
si_down <- clean_names(read.xlsx("~/Downloads/public_data/siRNAdown.xlsx", sheet = 1))

si_up$reg <- rep("up",nrow(si_up))
si_down$reg <- rep("down",nrow(si_down))
siRNA <- bind_rows(si_up, si_down)

gsea <- rename(gsea, gene_symbol = geneID)
inter <- inner_join(gsea, siRNA)
