setwd("~/GWAS_22/gwas_final/eQTL")

# snp pheno ----------

gwas <- fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus2/typhoid_tophits.txt")
geno <- fread("~/GWAS_22/gwas_final/eQTL/vast_tyg/Nov/combat/D0_f2r_.traw") # chr5_loci.traw
# Prep geno 
id_data <- clean_names(fread("~/GWAS_22/gwas_final/meta/geno_ids.csv"))
id_data$geno_id <- str_sub(id_data$cat_iid, start= -4)

plate <- read.ods(file = "~/GWAS_22/gwas_final/meta/vast_tyger_blood.ods", sheet = 1)
names(plate) <- plate[1,]
plate <- clean_names(plate)
plate <- plate[-1,]
cols <- plate[, c(6:ncol(plate))] %>% colnames
plate[, cols] <- lapply(plate[, cols], as.numeric)

snp <- geno[,-c(1,3:6)]
snp <-  t(snp)
colnames(snp) <- snp[1, ]
snp <- snp[-1, ]
snp <- as.data.frame(snp)
snp <- tibble::rownames_to_column(snp, var = "geno_id")
colnames(snp) <- str_replace_all(colnames(snp), "\\:", "_")
snp$geno_id <- str_sub(snp$geno_id, start= -4)
snp <- left_join(snp, id_data, by = "geno_id")
# Add pheno
plate$lab_id <- as.numeric(plate$id)

# Restart here if needed
snp_pheno <- left_join(snp, plate, by = "lab_id")
colnames(snp_pheno) <- str_replace_all(colnames(snp_pheno), "chr", "Chr")
cols <- snp_pheno %>% dplyr::select(starts_with("Chr")) %>% colnames()
snp_pheno[,cols] <- lapply(snp_pheno[,cols], factor)

# Add expression data -----------------
# from eqtl prep script
data <- readRDS(file = "combat_final.RDS")
# add iids to data #
ids <-  id_data %>% select(geno_id, lab_id, meta_id)
ids$lab_id <- as.character(ids$lab_id)
data$samples <- left_join(data$samples, ids, by = "lab_id")

idx <- (data$samples$time_point == "TD") # & data$samples$chall_vax =="T"
sub <- data[,idx]
dim(sub)

# idx <- sub$genes$gene_name %in% c("TNFSF13", "TNFSF12")
idx <- sub$genes$gene_name %in% c("F2R", "F2RL1", "AGGF1", "TBCA")
sub <- sub[idx,]
dim(sub)

exprs <- as.data.frame(t(sub$counts))
names(exprs) <- sub$genes$gene_name
# exprs <- tibble::rownames_to_column(exprs, var = "geno_id")
exprs$geno_id <- str_sub(sub$samples$geno_id)

# Join all tables
snp_pheno <-
  left_join(snp_pheno, exprs, by = "geno_id")

# Set up labels ----------------------------------------------------------------

# Label P-val and OR
gwas$SNP <- str_replace_all(gwas$SNP , "chr", "Chr")
gwas$SNP <- str_replace_all(gwas$SNP, ":", "_")
# Get summary data
lab <- gwas[SNP == "Chr5_77682126_C_T", OR]
lab <- signif(lab, digits = 2)
lab_p <- gwas[SNP == "Chr5_77682126_C_T", P]
lab_p <- signif(lab_p, digits = 2)
data <- snp_pheno %>% filter(!is.na(Chr5_77682126_C_T))

# TODO redo with the same snp (these are in high LD so results are the same and ok for now)
# lab2 <- gwas[SNP == "Chr12_19094052_C_T", OR]
# lab2 <- signif(lab2, digits = 2)
# lab_p2 <- gwas[SNP == "Chr12_19094052_C_T", P]
# lab_p2 <- signif(lab_p2, digits = 2)
# data2 <- snp_pheno %>% filter(!is.na(Chr12_19094052_C_T))

# Plot eQTLs ------------------------------------------------------------------

snp_names <-
  c("Chr5_77682126_C_T", "Chr5_77685638_A_G", "Chr5_77775879_G_T", "Chr5_77795103_G_A")

## F2R plot -----------------

y = snp_pheno$F2R
for (i in snp_names) {
  plt <-
    ggplot(snp_pheno, aes_string(x = i, y = y)) +
    geom_violin(outlier.shape = NA, aes_string(fill = i)) +
    scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
    geom_jitter(width = 0.1,
                size = 1.7,
                alpha = 0.5) +
    scale_colour_manual(values = c("#061c23"))  +
    labs(y = "Gene Expression\n", title = "TD") +
    theme_bw()  +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 11.5),
      axis.text.y = element_text(size = 11.5)
    )
  print(plt)
  Sys.sleep(0.5)
}

compare <- list(c("0", "1"))

## Diagnosis expression -----------------

# get 12 hour data or baseline data

data <- readRDS(file = "combat_final.RDS")
# add iids to data #
idx <- (data$samples$time_point == "D0.12h") # & data$samples$chall_vax =="T"
sub <- data[,idx]
dim(sub)

# idx <- sub$genes$gene_name %in% c("TNFSF13", "TNFSF12")
idx <- sub$genes$gene_name %in% c("F2R", "F2RL1", "AGGF1", "TBCA")
sub <- sub[idx,]
dim(sub)

exprs <- as.data.frame(t(sub$counts))
names(exprs) <- sub$genes$gene_name
exprs$lab_id <- str_sub(sub$samples$lab_id)
exprs$lab_id <- as.numeric(exprs$lab_id)
# Join all tables
pheno_exprs <-
  left_join(plate, exprs, by = "lab_id")

pheno_exprs  %>% filter(!is.na(outcome)) %>%
  ggplot(aes(x = outcome, y = F2R, fill = outcome)) +
  geom_violin()  +
  geom_jitter(width = 0.1,
              size = 1.7,
              alpha = 0.5) +
  scale_colour_manual(values = c("#061c23")) + theme_bw() +
  labs(y = "Expression (TPM)", title = "D0.12h", x = "Diagnosis") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 11.5),
    axis.text.y = element_text(size = 11.5)
  ) + scale_fill_viridis(discrete = TRUE) 

# 12h
# lower f2r gene expression assoc with protective alleles
# lower f2r expression at baseline/TD, less likely to have typhoid
# slight trend towards lower platelets

# TD
# one copy of protective allele, lower platelet gene xpressioin at day of diagnosis (check symptom cors)
# risk allele - higher platelet gene expression


# Diag by genotype

for(i in snp_names){
  plt <-
    snp_pheno %>% ggplot(aes_string(x="outcome", fill=i)) +
    geom_bar(stat = "count") + scale_fill_viridis(discrete = TRUE)+
    labs(y="Counts\n", x="Diagnosis", title = gsub("_", ":", i)) +
    theme_bw() +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(face="bold")) 
  print(plt)
  Sys.sleep(0.5)
}

# Risk allele (add OR data)
#Chr5_77682126_C_T", "Chr5_77685638_A_G

risk <-
  snp_pheno %>% ggplot(aes_string(x="outcome", fill= "Chr5_77682126_C_T")) +
  geom_bar(stat = "count") + scale_fill_viridis(discrete = TRUE)+
  labs(y="Counts\n", x="Diagnosis", title = gsub("_", ":", "Chr5_77682126_C_T")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(face="bold")) +
  annotate("text",
           label= paste("OR:",lab,"\n","P:",lab_p, sep = " "),
           y=65, x="No Typhoid", size=4, family="arial")

lab2 <- gwas[SNP == "Chr5_77685638_A_G", OR]
lab2 <- signif(lab2, digits = 2)
lab_p2 <- gwas[SNP == "Chr5_77685638_A_G", P]
lab_p2 <- signif(lab_p2, digits = 2)

protec <-
  snp_pheno %>% ggplot(aes_string(x="outcome", fill= "Chr5_77685638_A_G")) +
  geom_bar(stat = "count") + scale_fill_viridis(discrete = TRUE)+
  labs(y="", x="Diagnosis", title = gsub("_", ":", i)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(face="bold")) +
  annotate("text",
           label= paste("OR:",lab2,"\n","P:",lab_p2, sep = " "),
           y=65, x="No Typhoid", size=4, family="arial")

(risk|protec) + plot_annotation(tag_levels = "A", title = "Chr5 Risk Loci",) +
  plot_layout(guides = "collect")


save.image(file = "chr5_platelets.RData")



## Platelets plot ---------

y = snp_pheno$platelets_x10_9_l
for (i in snp_names) {
  plt <-
    ggplot(snp_pheno, aes_string(x = i, y = y)) +
    geom_boxplot(outlier.shape = NA, aes_string(fill = i)) +
    scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
    geom_jitter(width = 0.1,
                size = 1.7,
                alpha = 0.5) +
    scale_colour_manual(values = c("#061c23"))  +
    labs(y = "Platelet count", title = "TD") +
    theme_bw()  +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 11.5),
      axis.text.y = element_text(size = 11.5)
    )
  print(plt)
  Sys.sleep(0.5)
}


# Expression td ntd
exprs2 <- as.data.frame(t(data$counts))
names(exprs2) <- data$genes$gene_name
rownames(exprs2) <- data$samples$lab_id
exprs2$lab_id <- rownames(exprs2)
exprs2$lab_id <- as.integer(exprs2$lab_id)
exprs2 <- clean_names(exprs2)
exprs2 <- left_join(data$samples, exprs2, by = "lab_id")


exprs2 %>% filter(!is.na(diagnosis)) %>%
  ggplot(aes(x = diagnosis, y = f2r, fill = diagnosis)) +
  geom_violin()  +
  geom_jitter(width = 0.1,
              size = 1.7,
              alpha = 0.5) +
  scale_fill_manual(values = c("seagreen3","sandybrown")) + theme_bw() +
  labs(y = "Gene expression", title = "", x = "Diagnosis") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 11.5),
    axis.text.y = element_text(size = 11.5))



#+ stat_compare_means(
    comparisons = compare,
    method = "wilcox.test",
    label = "p.format",
    size = 5, 
    p.adjust.method = "none")






# old data ---------------------------------------------------------------------
  
  setwd("~/GWAS_22/gwas_final/eQTL/vast_tyg")
  snp_pheno <- read.csv(file = "VT_eqtl_TD.csv")
  gwas <- fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus2/typhoid_tophits.txt")
  geno <- fread("~/GWAS_22/gwas_final/eQTL/vast_tyg/TD.traw")
  
  cols <- snp_pheno %>% dplyr::select(starts_with("Chr")) %>% colnames()
  snp_pheno[,cols] <- lapply(snp_pheno[,cols], factor)
  snp_pheno$Diagnosed <- as.factor(snp_pheno$Diagnosed) 
  
  
  
  # Label P-val and OR
  gwas$SNP <- str_replace_all(gwas$SNP , "chr", "Chr")
  gwas$SNP <- str_replace_all(gwas$SNP, ":", "_")
  # Get summary data
  lab <- gwas[SNP == "Chr5_77682126_C_T", OR]
  lab <- signif(lab, digits = 2)
  lab_p <- gwas[SNP == "Chr5_77682126_C_T", P]
  lab_p <- signif(lab_p, digits = 2)
  data <- snp_pheno %>% filter(!is.na(Chr5_77682126_C_T))
  
  # TODO redo with the same snp (these are in high LD so results are the same and ok for now)
  lab2 <- gwas[SNP == "Chr12_19094052_C_T", OR]
  lab2 <- signif(lab2, digits = 2)
  lab_p2 <- gwas[SNP == "Chr12_19094052_C_T", P]
  lab_p2 <- signif(lab_p2, digits = 2)
  data2 <- snp_pheno %>% filter(!is.na(Chr12_19094052_C_T))
  

  
  compare <- list(c("0", "1"))
  
  F2R_td <-
    ggplot(data, aes(x = Chr5_77682126_C_T, y = F2R, fill = Chr5_77682126_C_T)) +
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
    labs(y = "Gene Expression (TPM)", title = "F2R") +
    theme_bw()  +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 11.5),
      axis.text.y = element_text(size = 11.5)
    ) +
    annotate(
      "text",
      label = paste("OR:", lab, "\n", "P:", lab_p, sep = " "),
      y = 31,
      x = 3,
      size = 4,
      family = "arial"
    )
  
  compare <- list(c("0", "1"), c("0", "2"))
  THBS1 <-
    ggplot(data2,
           aes(x = Chr12_19094052_C_T, y = THBS1, fill = Chr12_19094052_C_T)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
    geom_jitter(width = 0.1,
                size = 1.7,
                alpha = 0.4) +
    scale_colour_manual(values = c("#061c23")) +
    labs(y = " Gene Expression (TPM)", title = "THBS1", x = "Chr12_19094052_C_T") +
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
    ) +
    annotate(
      "text",
      label = paste("OR:", lab2, "\n", "P:", lab_p2, sep = " "),
      y = 35,
      x = 3,
      size = 4,
      family = "arial"
    )
  THBS1
  
  TD <- (vol)/(F2R_td |THBS1) + plot_annotation(tag_levels="a")
  plot(TD)
