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

# VAST directories -------------------------------------------------------------
setwd("~/GWAS_22/gwas_final/merge/typhoid/vast/iga")
plot_dir <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/plots/")
geno <- fread("~/GWAS_22/gwas_final/eQTL/vast/V0.traw")  # Make with geno_prep.sh script
pheno <- fread("~/GWAS_22/gwas_final/meta/vast_pheno.txt") # see vast_assoc_notes.R
covar <- fread("~/GWAS_22/gwas_final/meta/covar_vast.txt")
y_var <- c("Vi-IgA")
gwas <- fread("iga_tophits.txt")

# Set up variables and directories ---------------------------------------------

setwd("~/GWAS_22/gwas_final/merge/typhoid/assoc")
plot_dir <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/plots/")
geno <- fread("~/GWAS_22/gwas_final/eQTL/vast_tyg/TD.traw")  # Make with geno_prep.sh script
pheno <- fread("~/GWAS_22/gwas_final/meta/pheno_typhoid.txt") # see vast_assoc_notes.R
covar <- fread("~/GWAS_22/gwas_final/meta/covar_typhoid.txt")
y_var <- c("Expression (TPM)")
gwas <- fread("nexus2/typhoid_tophits.txt")
# ciber <- fread()
pheno_exprs <- read.csv(file = "~/RNA/oct/combat_vast_tyger_tpm.csv")
pheno_exprs <- filter(pheno_exprs, time_point3 == "TD")
exprs <- pheno_exprs[,str_detect(colnames(pheno_exprs),"ENS")]
exprs <- as.data.frame(t(exprs))
colnames(exprs) <- pheno_exprs$lab_id
exprs$ensembl_gene_id <- rownames(exprs)
id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids.csv")

# Genotyping data --------------------------------------------------------------
# capture charcters after the second instance of _ starting with

# clean IDs from colnames
# capture all the characters .*, at the start of ^ string unitl first instance ? of -, then replace with "" (i.e remove)
# colnames(geno) <- gsub(("^.*?-"), "", colnames(geno))

pheno$IID2 <- strc(pheno$FID, pheno$IID, sep = "_")

# transpose and just keep snp info
snp <- geno[,-c(1,3:6)]
snp <-  t(snp)
colnames(snp) <- snp[1, ]
snp <- snp[-1, ]
snp <- as.data.frame(snp)
snp <- tibble::rownames_to_column(snp, var = "IID2")
colnames(snp) <- str_replace_all(colnames(snp), "\\:", "_")

# Expression data
# qtls <- c("IL10RA", "TNFSF12", "IGHA2")
qtls <- c("F2R","MPL","ITGB3","MMRN1","GP9","GP1BA", "CD9", "THBS1", "ANGPT1","PTGS1")

qtls <- dplyr::filter(exprs, V1 %in% qtls)

qtls <- as.data.frame(t(qtls))
names(qtls) <- qtls[1,]
# 
qtls<- rownames_to_column(qtls, var = "FID")
qtls<- qtls[-1,]

pheno2 <- pheno
pheno2$FID <- stringr::str_sub(pheno$FID, start= -4)
snp2 <- snp 
snp2$FID <- stringr::str_sub(snp2$FID, start= -4)

# left join with pheno data
snp_pheno <-
  left_join(pheno2, snp2, by = "FID") %>%
  left_join(qtls, by = "FID")

# Make SNPs factor
cols <- colnames(dplyr::select(snp, -FID))

snp_pheno[,cols] <- lapply(snp_pheno[, ..cols], factor)
snp_pheno <- snp_pheno %>% filter(!is.na(snp_pheno[,30]))

snp_pheno$Diagnosed <-
  as.factor(snp_pheno$Diagnosed) 
snp_pheno$Diagnosed  <- fct_recode(snp_pheno$Diagnosed, nTD = "1", TD = "2")


# Plot -------------------------------------------------------------------------

# Basic Plot
snp_pheno %>%
  ggplot(aes(x=chr9_14597725_T_C, y=adnob)) +
  geom_boxplot(outlier.shape = NA, aes(fill =`chr9_14597725_T_C`)) +
  scale_fill_manual(values = c("#114B5F","#5C9EAD", "#E4FDE1")) +
  geom_jitter(aes(colour=Diagnosed), width = 0.1, size=1.5) +
  scale_colour_manual(values = c("seagreen3", "sandybrown")) +
  labs(y="Oxidative Burst Score") +
  theme_bw()  

# Without diagnosis colouring
snp_pheno %>%
  ggplot(aes(x=chr9_14597725_T_C, y=adnob)) +
  geom_boxplot(outlier.shape = NA, aes(fill =`chr9:14597725:T:C`)) +
  scale_fill_manual(values = c("#114B5F","#5C9EAD", "#E4FDE1")) +
  geom_jitter(width = 0.2, size=1.7) +
  scale_colour_manual(values = c("#061c23")) +
  labs(y="Phagocytic Score\n") +
  theme_bw() +
  theme(legend.position = "none")

# Loop snp v phenotype plots ---------------------------------------------------

# Set up vars
colnames(snp_pheno) <- str_replace_all(colnames(snp_pheno), "chr", "Chr")
snp_names <- str_subset(colnames(snp_pheno), "Chr4") 
y=snp_pheno$vi_ig_a_titre
y=as.numeric(snp_pheno$IGHA2)
# filter if sum of snp$pheno chr <4 for factor 2 
geno2 <- snp_pheno %>% dplyr::select(starts_with("Chr"))
melted = tbl_df(melt(df, id=c("owner")))
geno2 <- as.data.frame(t(geno2))
summary <- lapply(geno2, fct_count)

# Run loop
for(i in snp_names){
  plt <-
    ggplot(snp_pheno, aes_string(x=i, y=y)) +
    geom_boxplot(outlier.shape = NA, aes_string(fill =i)) +
    scale_fill_manual(values = c("#114B5F", "#5C9EAD", "#E4FDE1")) +
    geom_jitter(width = 0.2, size=1.7) +
    scale_colour_manual(values = c("#061c23")) +
    labs(y="IGHA2 Expression") +
    theme_bw() +
    theme(legend.position = "none")
  print(plt)
  Sys.sleep(0.5)
}

save.image(file = "VAST_iga.RData")


# Gene expression
qtls<- dplyr::filter(exprs, V1 == "TNFSF12")


eqtlsnp <- 


# Diagnosis - categorical plots -----------------------------------------------

snp_pheno$diagnosis <- as.factor(snp_pheno$diagnosis)

# Basic plot
snp_pheno %>% ggplot(aes(x=diagnosis, fill=Chr1_161607867_G_A)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("#114B5F", "#5C9EAD", "#E4FDE1")) +
  labs(y="Counts\n", x="Diagnosis", title = gsub("_", ":", "Chr1_161607867_G_A")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(face="bold"))

# Loop plot
for(i in snp_names){
  plt <-
    snp_pheno %>% ggplot(aes_string(x="diagnosis", fill=i)) +
    geom_bar(stat = "count") +
    scale_fill_manual(values = c("#114B5F", "#5C9EAD", "#E4FDE1")) +
    labs(y="Counts\n", x="Diagnosis", title = gsub("_", ":", i)) +
    theme_bw() +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(face="bold")) 
  print(plt)
  Sys.sleep(0.5)
}


# Label P-val and OR
gwas$ID <- str_replace_all(gwas$ID , "chr", "Chr")
gwas$ID <- str_replace_all(gwas$ID, ":", "_")
gwas <- rename(gwas, OR = BETA)

#test <- snp_names[1:2]
snp_risk <- gwas[(OR>2),ID]

for(i in snp_risk){
  
  # Get summary data
  lab <- gwas[ID == i,OR]
  lab <- signif(lab, digits = 2)
  lab_p <- gwas[ID == i,P]
  lab_p <- signif(lab_p, digits = 2)
  # Plot bar chart
  plt <-
    snp_pheno %>% ggplot(aes_string(x="diagnosis", fill=i)) +
    geom_bar(stat = "count") +
    scale_fill_manual(name = "Genotype",
                      values = c("#114B5F", "#5C9EAD", "#E4FDE1")) +
    labs(y="Counts\n", x="Diagnosis", title = gsub("_", ":", i)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12),
          axis.title = element_text(size =11.5),
          axis.text.y = element_text(size =11.5)) +
    annotate("text",
             label= paste("OR:",lab,"\n","P:",lab_p, sep = " "),
             y=35, x="TD", size=4, family="arial")
  print(plt)
  Sys.sleep(0.5)
}


# With pvals --------------------------------------

compare <- list(c("0","1"), c("1","2"), c("0", "2"))

for(i in snp_names){
  plt <-
    ggplot(snp_pheno, aes_string(x=i, y=y)) +
    geom_boxplot(outlier.shape = NA, aes_string(fill =i)) +
    scale_fill_manual(values = c("#114B5F", "#5C9EAD", "#E4FDE1")) +
    geom_jitter(width = 0.2, size=1.7) +
    scale_colour_manual(values = c("#061c23")) +
    labs(y="Phagocytic Score\n") +
    theme_bw() +
    theme(legend.position = "none") +
    stat_compare_means(comparisons = compare,
                       method = "wilcox.test",
                       label = "p.format",
                       p.adjust.method = "none")
  print(plt)
  Sys.sleep(0.5)
}

save.image(file="iga_eqtl-snp.RData")


# Save to PDF -----------------------------------------------------------------
#TODO work this out

# Open PDF
pdf(file="diagnosis_plots.pdf")
# Set up plot grid
par(mfrow = c(2,2))
# Run loop
for(i in snp_names){
  plt <-
    ggplot(snp_pheno, aes_string(x=i, y=y)) +
    geom_boxplot(outlier.shape = NA, aes_string(fill =i)) +
    scale_fill_manual(values = c("#114B5F", "#5C9EAD", "#E4FDE1")) +
    geom_jitter(width = 0.2, size=1.7) +
    scale_colour_manual(values = c("#061c23")) +
    labs(y="Phagocytic Score\n") +
    theme_bw() +
    theme(legend.position = "none")
  print(plt)
}


# Make plots.
plot_list = list()
for(i in snp_names){
  p = ggplot(snp_pheno, aes_string(x=i, y=y)) +
    geom_boxplot(outlier.shape = NA, aes_string(fill =i)) +
    scale_fill_manual(values = c("#114B5F", "#5C9EAD", "#E4FDE1")) +
    geom_jitter(width = 0.2, size=1.7) +
    scale_colour_manual(values = c("#061c23")) +
    labs(y="Phagocytic Score\n") +
    theme_bw() +
    theme(legend.position = "none")
  plot_list[[i]] = p
}


# Another option: create pdf where each page is a separate plot.
pdf("diagnosis_plots.pdf")
for (i in 1:3) {
  print(plot_list[[i]])
}
dev.off()






