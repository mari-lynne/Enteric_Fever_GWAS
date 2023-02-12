# fisher enrichment test

# Compare vax-gwas-degs with vax-degs
# Also are they enriched in infection degs at D1 or D7

ego <- enrichGO(gene = genes_ez$entrez,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", readable = TRUE)

View(ego@result)

snp_ego <- ego

# get ez ids of gwas-dge  genes
src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
snp_genes <- unique(top_tab$gene_name) #33
genes_ez <- select(src, 
                   keys = snp_genes,
                   columns = c("symbol","entrez"),
                   keytype = "symbol")

genes_ez <- distinct(genes_ez)

# get ez ids of background genes
sig_v1 <- filter(sig_v1, abs(logFC) >= 0.1)
sig_v7 <- filter(sig_v7, abs(logFC) >= 0.1)
sig_bkg <- bind_rows(sig_v1, sig_v7)

bkg_genes <- unique(sig_bkg$gene_name)
genes_bk <- select(src, 
                   keys = bkg_genes,
                   columns = c("symbol","entrez"),
                   keytype = "symbol")

genes_bk <- distinct(genes_bk)

library(DOSE)
enrich_gwas_dge <- DOSE::enricher_internal(
  gene = genes_ez$entrez,
  pvalueCutoff,
  pAdjustMethod = "BH",
  universe = genes_bk$entrez,
  minGSSize = 3,
  maxGSSize = 500,
  qvalueCutoff = 0.2
  )








