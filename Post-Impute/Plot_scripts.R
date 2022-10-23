# Genotype plotting

# Plot genotype count on x, with antibody response/expression/on y (linear associations)
# Plot genotype on x, with case control freq on y

# Run ~/GWAS_22/Enteric_GWAS/DEG-analysis/eQTL/2.geno_snp_prep.sh script to make traw file

# Packages ---------------------------------------------------------------------
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(tidylog)
library(stringr)
library(stringi)
library(openxlsx)
library(RColorBrewer)
library(forcats)
library(qqman)
library(janitor)
library(magrittr)

# Set up variables and directories ---------------------------------------------

setwd("~/GWAS_22/gwas_final/merge/typhoid/assoc/vast/diagnosis")
plot_dir <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/plots/")
geno <- fread("vast_diagnosis.traw")  # Make with geno_prep.sh script
pheno <- fread("~/GWAS_22/gwas_final/meta/vast_pheno.txt") # see vast_assoc_notes.R
covar <- fread("~/GWAS_22/gwas_final/meta/covar_typhoid.txt")
y_var <- c("Vi-adnp")
gwas <- fread("tophits_diagnosis.txt")


# Genotyping data --------------------------------------------------------------
#capture charcters after the second instance of _ starting with

# clean IDs from colnames
# capture all the characters .*, at the start of ^ string unitl first instance ? of -, then replace with "" (i.e remove)
colnames(geno) <- gsub(("^.*?-"), "", colnames(geno))

# transpose and just keep snp info
snp <- geno[,-c(1,3:6)]
snp <-  t(snp)
colnames(snp) <- snp[1, ]
snp <- snp[-1, ]
snp <- as.data.frame(snp)
snp$FID <- row.names(snp)
colnames(snp) <- str_replace_all(colnames(snp), "\\:", "_")

# left join with pheno data
snp_pheno <-
  left_join(pheno, snp, by = "FID") %>%
  left_join(covar, by = "IID")

# Make SNPs factor
cols <- colnames(select(snp, -FID))

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
snp_names <- str_subset(colnames(snp_pheno), "Chr") 
y=snp_pheno$adnob

# Run loop
for(i in snp_names){
  plt <-
    ggplot(snp_pheno, aes_string(x=i, y=y)) +
    geom_boxplot(outlier.shape = NA, aes_string(fill =i)) +
    scale_fill_manual(values = c("#114B5F", "#5C9EAD", "#E4FDE1")) +
    geom_jitter(width = 0.2, size=1.7) +
    scale_colour_manual(values = c("#061c23")) +
    labs(y="Oxidative Burst Score\n") +
    theme_bw() +
    theme(legend.position = "none")
  print(plt)
  Sys.sleep(0.5)
}


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

save.image(file="diagnosis-snp.RData")


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






