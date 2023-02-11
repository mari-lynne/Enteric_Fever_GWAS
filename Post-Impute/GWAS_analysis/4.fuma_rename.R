setwd('~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma')

library(data.table)
library(tidylog)
library(dplyr)
library(stringr)
library(openxlsx)
library(biomaRt)

fuma <- fread("typhoid_37_gwas.txt")
colnames(fuma) <- c("CHROM","POS","A2","A1","P","OR","SE")
write.table(fuma, "typhoid_fuma.txt", row.names = F, quote = F, sep = "\t")

system("gzip ~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma/typhoid_fuma.txt")

DT <- data.table(filter(fuma, P>1e-5 & P<1e-4))
top_snps <- DT[, .SD[which.min(P)], by = CHROM]
top_snps <- select(top_snps, CHR, POS)

write.table(lead, "typhoid_extra_lead.txt", row.names = F, quote = F, sep = "\t")

snp_mart <- useEnsembl(biomart="snps", 
                       host="grch37.ensembl.org", 
                       dataset="hsapiens_snp")


## combine the positions in to a single vector
position <- str_c(top_snps$CHROM, top_snps$POS, top_snps$POS, sep = ":")
position <- c(position, "5:154144110:154144110")
position

getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start'), 
      filters = 'chromosomal_region', 
      values = position, 
      mart = snp_mart)


library(biomaRt)
library(AnnotationDbi)
library(Organism.dplyr)
library(TxDb.Hsapiens.UCSC.hg37.knownGene)

src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
chrom <- c(as.character(unique(gwas_c$CHR)))
chrom <- str_c(rep("chr", length(chrom)), chrom)
genes_ez <- select(src, 
                   keys = genes,
                   columns = c("symbol","entrez"),
                   keytype = "symbol")









