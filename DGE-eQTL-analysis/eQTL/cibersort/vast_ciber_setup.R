# Cibersort ####

library(sva)


# Filter time expression data -------------------------------------------------

setwd("~/GWAS_22/gwas_final/eQTL")
study <- c("/vast")
out_dir <- paste0(getwd(), study)
'%!in%' <- function(x,y)!('%in%'(x,y))

pathway <- c("vast/iga_genes.csv") # "~/GWAS_22/gwas_final/InnateDB_genes.csv"
id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids.csv")
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/vast/vast.fam")
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/vast/iga/iga_tophits.txt")
covar_data <- c("~/GWAS_22/gwas_final/meta/covar_vast.txt")
load(file = "vast_filtered_raw.RData")
str(data$samples)
data$meta_data <- select(data$meta_data, -x, -group, -y, -vaccine_y, -screening_id_y)
data$samples <- bind_cols(data$samples, data$meta_data)

# filter just vaccine samples
idx <- data$samples$time_point3 %in% c("V0", "V1", "V7")
sub <- data[,idx]
dim(sub)

# Combat batch adjustment -----------------------------------------------------
levels(covar$vaccine)
sub$samples <- droplevels(sub$samples)
covar <- model.matrix(~(study_arm), data =sub$samples)

# Get expression data
combat <- sva::ComBat_seq(counts=sub$counts,
                          batch=sub$samples$sequence_pool,
                          group=sub$samples$time_point3,
                          covar_mod = covar)

# TPM --------------------------------------------------------------------------

# get correct gene fData info

tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

tpm <- tpm3(sub$counts, sub$genes$transcript_length)
exprs <- as.data.frame(tpm)

sub$raw_counts <- sub$counts
sub$counts <- tpm

# Biomart data
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','start_position','end_position'), mart = ensembl) %>% rename(gene_id = ensembl_gene_id)

sub$genes <- inner_join(sub$genes, genes, by = "gene_id")
idx <- rownames(sub) %in% sub$genes$gene_id
sub <- sub[idx,]
dim(sub$genes)

saveRDS(sub, file = "VAST_tpm_nocomb.rds")





# Immune genes (for vaccination study)
exprs <- genes %>% inner_join(exprs, by = "ensembl_gene_id") %>% filter(chromosome_name %in% (as.character(c(1:22))))

# Filter NAs
exprs[exprs == ""]<- NA
exprs <- filter(exprs, !is.na(hgnc_symbol))

# Prep mixture file ------------------------------------------------------------

pheno <- data$meta_data
pheno <- pheno %>%
  mutate(time_point = ifelse(time_point3 == "D0" & (study_arm %in% c("ViTCV", "ViPS")), "V28", time_point3))
pheno <- pheno %>%
  mutate(time_point = ifelse(time_point %in% c("V0","D0"), "Baseline", time_point)) 

pheno_exprs <- pheno_exprs %>%
  dplyr::filter(time_point == "Baseline")

pheno$time_id <- str_c(pheno$lab_id, pheno$arm_time, sep = "_")
colnames(exprs) <- c("Gene", pheno$time_id)

# Filter time point
V0 <- exprs[ ,grep("-28|Gene", colnames(exprs))]
V1 <- exprs[ ,grep("-27|Gene", colnames(exprs))]
V7 <- exprs[ ,grep("-21|Gene", colnames(exprs))]

# Make tables
write.table(V0, file = "vast/V0_cibersort.txt", sep = "\t", quote = FALSE, row.names = F)
write.table(V1, file = "vast/V1_cibersort.txt", sep = "\t", quote = FALSE, row.names = F)
write.table(V7, file = "vast/V7_cibersort.txt", sep = "\t", quote = FALSE, row.names = F)

# Pheno exprs with combat adjusted, tpm data
exprs <- as.data.frame(t(exprs))
names(exprs) <- exprs[1,]
exprs <- exprs[-1,]
exprs$row_names <- pheno$row_names
pheno_exprs <- left_join(pheno, exprs, by = "row_names")

write.csv(pheno_exprs, file = "vast/vast_tcpm.csv")

# Start Analysis ---------------------------------------------------------------
pheno_exprs <- read.csv(file = "vast/vast_tcpm.csv")

# Make new baseline time point, V0 for vaccinees and D0 for challenge
pheno_exprs <- pheno_exprs %>%
  mutate(time_point = ifelse(time_point3 == "D0" & (study_arm %in% c("ViTCV", "ViPS")), "V28", time_point3))
pheno_exprs <- pheno_exprs %>%
  mutate(time_point = ifelse(time_point %in% c("V0","D0"), "Baseline", time_point)) 
pheno_exprs <- pheno_exprs %>%
  dplyr::filter(time_point == "Baseline")

exprs <- pheno_exprs[,-c(1:46)]
exprs <- as.data.frame(t(exprs))

pheno_exprs$time_id <- str_c(pheno_exprs$lab_id, pheno_exprs$study_arm, pheno_exprs$time_point, sep = "_")
colnames(exprs) <- c(pheno_exprs$time_id)
write.table(pheno_exprs, file = "vast/baseline_cibersort.txt", sep = "\t", quote = FALSE, row.names = F)

gene_list <- c("IGHA2", "TNFRSF13B", "TNFSF13", "TNFRSF13")

write.table(gene_list, file = "vast/iga_gene_list.txt", sep = "\t", quote = FALSE, row.names = F, col.names = F)

# Results ---------------------------------------------------------------------

# Cell -specific --------------------------------------------------------------

setwd("~/GWAS_22/gwas_final/eQTL/vast/final_eqtls")

# Day 7
cell7 <- fread("V7_cibersort/CIBERSORTxHiRes_Job22_Plasmacells_Window34.txt")
time <- c("V7")
type <- c("B-cell")
# Join tables
cell7 <- melt(cell7, id.vars = c("GeneSymbol"), measure.vars = c(2:ncol(cell)))
colnames(cell7) <- c("Gene","time_id", "Expression")
cell7$cell_type <- rep(type, nrow(cell7))
cell7$time_point <- rep(time, nrow(cell7))

# Day 1
cell1 <- melt(cell1, id.vars = c("GeneSymbol"), measure.vars = c(2:ncol(cell1)))
colnames(cell1) <- c("Gene","time_id", "Expression")
cell1$cell_type <- rep("cellrophil", nrow(cell1))
cell1$time_point <- rep("D1", nrow(cell1))


# Baseline
cell7 <- fread("V0_cibersort/CIBERSORTxHiRes_Job23_Plasmacells_Window34.txt")
time <- c("V0")
type <- c("B-cell")
cell0 <- melt(cell0, id.vars = c("GeneSymbol"), measure.vars = c(2:ncol(cell0)))
colnames(cell0) <- c("Gene","time_id", "Expression")
cell0$cell_type <- rep("cellrophil", nrow(cell0))
cell0$time_point <- rep("Baseline", nrow(cell0))

cell <- bind_rows(cell0, cell1, cell)

# Plot 
cell %>% group_by(Gene) %>% summarise(avg = mean(Expression)) %>%
  ggplot(aes(x = Gene, y = avg, fill = Gene)) +
  geom_bar(width = 1, stat = "identity") +
  labs(title = "D7 - Plasma Cell", x = "\nFc-Receptor", y = "Gene Expression (average)\n") +
  scale_fill_viridis_d() + theme_pubr() +
  theme(legend.position = "none")

cell  %>% group_by(Gene, time_point) %>% summarise(avg = mean(Expression)) %>%
  ggplot(aes(x = Gene, y = avg, fill = Gene)) +
  geom_bar(width = 1, stat = "identity") +
  labs(title = "cellrophil FcR-Expression", x = "\nFc-Receptor", y = "Gene Expression (average)\n") +
  scale_fill_viridis_d() +
  theme(legend.position = "none") + facet_wrap(~time_point)