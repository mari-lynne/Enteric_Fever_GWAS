# Fuma output typhoid fever -----------------
setwd("~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma")

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

'%!in%' <- function(x,y)!('%in%'(x,y))

## EQTL fuma --------

dat1 <- fread(file = "~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/FUMA_extra2/eqtl.txt")
dat2 <- fread(file = "~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/FUMA_job/eqtl.txt")
dat3 <- fread(file = "~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/FUMA_extra/eqtl.txt")
# Analysis just had tweaked search, do from the get go next time
all_qtl <- bind_rows(dat1, dat2, dat3)

table(all_qtl$tissue)
# Filter for non-relavant qtl tissues
data <- all_qtl %>% filter(tissue %!in% c("TwinsUK_ge_skin", "TwinsUK_ge_fat", "GEUVADIS_ge_LCL", "TwinsUK_ge_LCL", "BrainSeq_ge_brain"))

# Get unique eqtl-mapped genes
data <- data.table(data)
top_qtl <- data[, .SD[which.min(p)], by = symbol] # 62
# TODO add clean files from env that end in rm or a num
rm(dat1, dat2, dat3)

# Salmonella qtls
stim_all <- all_qtl %>% filter(str_detect(tissue, "Salmonella|Listeria|LPS|R848|IFN|IAV|Pam3CSK4")) #102
stim <- data.table(stim_all)
stim <- stim[, .SD[which.min(p)], by = symbol] 

# Mapped genes
# top_gs <- bind_rows(top_gs, data)
# dat3 <- fread(file = "~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/FUMA_typhoid_og/genes.txt") 
# DT <- data.table(top_gs)
# top_qene <- DT[, .SD[which.min(minGwasP)], by = symbol]

## Magma -------
magma1 <- fread(file = "~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/FUMA_typhoid_og/magma.genes.out")
magma2 <- fread(file = "~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/FUMA_extra/magma.genes.out")
magma <- bind_rows(magma1, magma2)
# Filter sig
magma <- filter(magma, P <0.05) 

magma <- data.table(magma)
magma <- magma[, .SD[which.min(P)], by = GENE] # just keep one of each gene
magma <- clean_names(magma)

gene_list <- unique(c(magma$gene, top_qtl$gene))
snp_list <-unique(c(magma$symbol, top_qtl$symbol))

# Top table --------------------------------------------------------------------
# Filter TD top table or TD nTD at day 7 (D7_TDnTD)
# adj.P.Val <= 0.01, h12 - relax to 0.0 unadjust

#h12_TDnTD <- read.csv(file = "final_top.csv")
D7_TDnTD <- read.csv(file = "deg_D7_TDnTD.csv")

top_tab_all <- D7_TDnTD %>% mutate(sig_pos = row_number()) %>% filter(P.Value <= 0.05 & (logFC < -0.2 | logFC >0.2))
top_tab <- D7_TDnTD %>% mutate(sig_pos = row_number()) %>%
  filter(adj.P.Val <= 0.05 & gene_id %in% gene_list & (logFC < -0.2 | logFC >0.2))

top_tab <- top_tab %>% mutate(eqtl = ifelse(gene_id %in% top_qtl$gene, "qtl", "pos"),
                              reg = ifelse(sign(top_tab$logFC) == "1", "TD", "nTD")) # pos indicates higher Fc in TD compared
# early_top <- top_tab
# final_top <- top_tab
# top_tab <- bind_rows(final_top, early_top) #n=204

# STAT1(98% of iterations),SLAMF8(76%),PSME2(39%),WARS(37%) and ALDH1A1(36%; Fig 4A)
#write.csv(top_tab, file = "final_top.csv", row.names = F)

## Background -------------------------------------------------------------------
# check what wasn't detected

# keep intersecting genes only
gene_list <-  gene_list[gene_list %in% D7_TDnTD$gene_id]

check1 <- top_qtl[top_qtl$symbol %in% D7_TDnTD$gene_name] #31
check2 <- magma[magma$symbol %!in% D7_TDnTD$gene_name]#36 - 310 = 346 of og genes not measured
miss <- bind_rows(check, check2)
background_snps <- gene_list[gene_list %!in% miss$gene]
background_mag <- magma[magma$gene %in% background_snps]
background_qtls <- top_qtl[top_qtl$gene %in% background_snps]

background_genes <-  D7_TDnTD %>% filter(adj.P.Val <= 0.05 & (logFC < -0.2 | logFC >0.2 ))

# so 415 out of 759 annotated genes could be measured in top table - 251 of 415 snps in toptable
# top table 5722 diff expressed genes
# 7.5~ top genes were also in our annotated snps

spread <- list('final_top' = final_top,
               'back_snps' = background_mag,
               'back_qtl' = background_qtls,
               'back_topgene' = background_genes,
               'back_allgene' = all,
               'top_magma' = magma,
               'top_qtls' = top_qene)

# write.xlsx(spread, file = "snp_deg_all.xlsx", rowNames = F)

#top_snp <- unique(magma$gene) # unique snps (magma has multiple gene pairs)

## Sig values ---------
top_tab <- rename(top_tab, symbol = gene_name)
sig_qtl <- tidylog::inner_join(top_qtl, top_tab, by = "symbol")
sig_mag <- tidylog::inner_join(magma, top_tab, by = "symbol")
sig_all <- unique(c(sig_qtl$symbol, sig_mag$symbol))


spread <- list('sig_qtl' = sig_qtl,
               'sig_mag' = sig_mag)
write.xlsx(spread, file = "snp_deg_sig2.xlsx", rowNames = F)

sig <- sig_qtl %>%
  select(gene, uniqID, db, tissue, p, signed_stats, alignedDirection) %>%
  rename(qtl_p = p)

sig_m <- sig_mag %>% select(gene, nsnps, p) %>% rename(gene_p = p)
rank <- left_join(top_tab, sig, by = "gene") %>% left_join(sig_m, by = "gene")
gsea <- rank %>% select(geneID, logFC, P.Value)
write.csv(gsea, file = "gsea.csv", row.names = F)

# GWAS sum stats -------------------------------------------
gwas <- clean_names(fread("typhoid_37_gwas.txt"))
gwas <- filter(gwas, p <= 1e-4) #161
gwas <- rename(gwas, chr = number_chrom)

names(rank)
rank <- rank %>% select(-c("gene_biotype","GC_content","geneID","t","signed_stats","alignedDirection","db","gene_p"))

# Get top variant magma and match to table
genes4 <- fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/FUMA_typhoid_og/genes.txt")
genes1 <- fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/FUMA_extra/genes.txt")
genes2 <- fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/FUMA_extra2/genes.txt")
genes3 <- fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/FUMA_job/genes.txt")

# Analysis just had tweaked search, do from the get go next time
all_gen <- data.table(bind_rows(genes1, genes2, genes3, genes4))
all_gen <- rename(all_gen, gene_name = symbol)
all_gen <- filter(all_gen, !stri_duplicated(gene_name)) %>% filter(type == "protein_coding")#1186 snps

# Join GWAS - then rank
#all_gen <- left_join(gwas, all_gen, by = "gene_name") #161 #by = c("chr", "pos"
merge <- inner_join(rank, all_gen, by = "gene_name")
write.csv(merge, file = "merge-data.csv")
merge$db_snp <- merge$IndSigSNPs
merge <- merge %>% mutate(db_snp  = sub(";.*", "",merge$IndSigSNPs))

# Get ancestry: rsID - 
# SNP Nexus ------
top <-clean_names(fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus2/typhoid_gwas.txt"))
top <- top %>% filter(p <= 1e-3) %>% rename(position = bp)

near <- clean_names(fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus2/near_gens.txt"))
kg <- clean_names(fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus2/1KGen.txt"))
#kg <- left_join(near, kg, by = "variation_id") %>% select(!contains(".y"))
kg <- left_join(kg, near, by = c("chromosome", "position", "variation_id")) %>% select(!contains(".y")) %>% rename(gene_name = overlapped_gene) %>% left_join(top, by = "position") %>% filter(!stri_duplicated(db_snp))

merge2 <- inner_join(merge, kg, by = "db_snp") %>% rename(symbol = gene_name.x)

extra <- inner_join(stim_all, merge2, by = "symbol")
write.csv(merge2, file = "sum2.csv")
write.csv(extra, file = "extra2.csv")

sum <- clean_names(read.csv(file = "final_snp_gene_sum.csv"))
num_cols <- names(select_if(sum, is.numeric))
sum[,num_cols]
sum[,num_cols]<- unlist(lapply(sum[,num_cols], function(x) {signif(x, digits = 2)}))
sum <- data.frame(lapply(sum, function(x) {gsub("_", " ", x)}))
sum <- sum %>% rename(Ref = ref_allele,
                      Alt = minor_allele,
                      Chr = chr_loci)
sum$tissue <- str_replace_all(sum$tissue ,"Alasoo 2018 ge macrophage", "Alasoo Macrophage")

write.csv(sum, file = "clean_sum.csv", row.names = F)


# ov <- kg %>% inner_join(rank, by = "gene_name") #8
# down <- kg %>% rename(gene_name = nearest_downstream_gene) %>% inner_join(rank, by = "gene_name") #3
# up  <- kg %>% rename(gene_name = nearest_upstream_gene) %>% inner_join(rank, by = "gene_name")

top_gen <- all_gen[, .SD[which.min(minGwasP)], by = symbol] # 62
# TODO add clean files from env that end in rm or a num
rm(genes1, genes2, genes3)
top_gen <- filter(top_gen, type == "protein_coding") %>% rename(gene_name = symbol)#52

# Then bind min gwas sig or min eqtl sig with 
merge <- left_join(rank, top_gen, by = "gene_name")
merge <- rank %>% rename(chr = chromosome_name) %>% left_join(all_gen, by = c("chr", "pos"))

# VOlCANO VIEW -----------------------------------------------------------------
deg <- TD_vs_D0
snp_genes <- unique(top_tab$gene_name)

plot_vol(deg,
              lab = 'top_lab',
              genes = snp_genes,
              title = "TD vs nTD 7-days post infection")

qtls <- unique(top_qtl$symbol)
plot_vol(deg,
         lab = 'top_lab',
         genes = qtls,
         title = "TD vs nTD 7-days post infection", top = 5)


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
ghe_genes <- D7_TDnTD %>% filter(adj.P.Val <= 0.05 & (logFC < -0.1 | logFC >0.1))
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

# Enrichment -------------------------------------------------------------------
load(file = "snp_dge.RData")

## Get enterez ids -----------
# data <-  dge_genes,
# data <- snp_list
# data <- unique(rank$gene_name) / top_qtl
# snp_genes = common -> snp_cross
# keys = sig_all$symbol

snp_genes <- sig_all

genes_ez <- select(src, 
                   keys = snp_genes,
                   columns = c("symbol","entrez"),
                   keytype = "symbol")

genes_ez <- distinct(genes_ez)

## ) Run ernrichGO code -------

ego <- enrichGO(gene = genes_ez$entrez,
                OrgDb = org.Hs.eg.db,
                universe = genes_bg$entrez,
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
dotplot(dge_ego, showCategory=20) + ggtitle("dotplot for GSEA")


## MsigDB -------

c7 <- read.gmt("~/Downs/public_data/c7.all.v2022.1.Hs.entrez.gmt")
immune_ego  <- enricher(genes_ez$entrez, TERM2GENE=c7)
barplot(immune_ego, showCategory=20)
msig <- subset(msig_ego@result, p.adjust <0.05)

# update descriptions
msig_ego@result$Description <- str_replace_all(msig_ego@result$Description,
                                              "_", " ")
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

dotplot(msig_ego, showCategory=15) + ggtitle("Immune GSEA: GWAS-DGE genes")

barplot(msig_ego, showCategory=15)


# save genes
write.csv(genes_ez, file = "cross_genes_ez.csv", row.names = F)

msig_ez <- unlist(strsplit(msig$geneID, "/"))

msig_ez <- select(src, 
                   keys = msig_ez,
                   columns = c("symbol","entrez"),
                   keytype = "entrez")

msig_ez <- distinct(msig_ez) %>% rename(gene_name = symbol)
# get joining snp data too
rank <- rank %>% dplyr::select(-inter) %>% left_join(msig_ez, by = "gene_name")
write.csv(rank, file = "cross_rank_ez.csv", row.names = F)

ont_genes <- subset(rank, !is.na(entrez))
write.csv(ont_genes, file = "ont_genes.csv", row.names = F)

## link2 loci ----

# For eQTL analysis, test all GWAS sig snps in loci, with msig genes
# so we need a gene name list to test
# and then also filter top snps for snps in magma (but could just try without first)

loci <- fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/FUMA_extra2/GenomicRiskLoci.txt")

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

deg <- D7_TDnTD
loci_genes <- unique(sig_loci$gene_name)

vol <- plot_vol(deg,
         lab = 'top_lab',
         genes = loci_genes,
         p = 0.01,
         title = "7-days post challenge (TD vs nTD)", top = 6,
         colours = c("#dd513a", "#932667", "#420a68"))

vol/(b|c) + plot_annotation(tag_levels = "a")

write.csv(D7_TDnTD, file = "deg_D7_TDnTD.csv", row.names = F)







# old ---------------------------------
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens")
em <- enricher(genes_ez$entrez, TERM2GENE=m_df)
head(em)



str(eqtl3)
table(eqtl3$tissue)
sigeqtl <- filter(eqtl1, FDR <= 0.05)
unique(eqtl3$symbol)
platelet <- sigeqtl %>% filter(tissue == "CEDAR_platelet") # SORD gene
pubmed <-   fread(file = "~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/FUMA_extra2/gwascatalog.txt")
blood <- pubmed[str_detect(pubmed$Trait, "platelet|Platelet|blood|count|volume|cell"), ]
blood_p <- pubmed[str_detect(pubmed$Trait, "platelet|Platelet"), ]

genes <- unique(sigeqtl$symbol)
genes

# Plan
# what are immmune eqtls rs ids
# investigate those snps and public snps for eqtl analysis in my data

immune <- read.csv
