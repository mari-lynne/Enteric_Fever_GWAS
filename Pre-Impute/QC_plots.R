
#Plotting code

library(data.table)
library(stringr)
library(dplyr)
library(plinkQC)
library(ggplot2)
library(patchwork)

#Aims
#Generate QC plots or summaries,for pre-imputation and post-imputation data

#File names 
#P1Tyger og files
#P1Tyger_remap QC - pre imputed files
#- P1Tyger, T1T2 merged > P1T1. Upload separately from VAST_PATCH study to topmed

#Pre-imp thresholds ####
#./plink --bfile P1Tyger_remap --geno 0.1 --mind 0.1 --make-bed --out P1Tyger_remap
#./plink --bfile P1Tyger_remap --geno 0.05 --make-bed --out P1Tyger_remap
#./plink --bfile T1T2_remap --geno 0.1 --mind 0.1 --make-bed --out T1T2_remap
#./plink --bfile T1T2_remap --geno 0.05 --make-bed --out T1T2_remap
#./plink --bfile VAST_remap --geno 0.1 --mind 0.1 --make-bed --out VAST_remap #2 people removed due to --geno #should have done this after second filtering step
#./plink --bfile VAST_remap --geno 0.05 --make-bed --out VAST_remap


#Imputation thresholds ####
##QC imputed data in plink 
#Details at https://github.com/mari-lynne/Enteric_Fever_GWAS/blob/main/Post-Impute/Post-imputation_QC/Post_Impute_QC.Rmd

#plink2 --bfile VAST_concat --geno 0.1 --mind 0.1 --make-bed --out VAST_imp_QC (10% cut off)
#plink2 --bfile VAST_imp_QC --geno 0.05 --mind 0.05 --maf 0.01 --make-bed --out VAST_imp_QC2

#VAST Imp results 
#16873950 variants loaded from VAST_imp_QC.bim.
#7396147 variants removed due to allele frequency threshold(s)
#947,7803 variants remaining after main filters.

#P1_T1 results #
#7944509 variants removed due to allele frequency threshold(s)
#895,0259 variants remaining after main filters.
#((16894768 - 7944509)/16894768)*100 = 52% removed

#Plans ####

#Pre-impute plots/summaries 
#1) Per marker quality control
#2) Individual missingness (2 people removed from VAST_PATCH study)
#remove snps with low MAF and HWE

#Post-Imputation plots/summaries
#QC stats, geno, maf, hwe, R2 from top med 
#Sex check, ancesttry, IBD


#Pre-impute plots ####

#Start with T1T2

Need plots for MIND/Geno thresholds, IBD, sex_check

# Summary:
#Using Plink we ran an initial filtering step for snps of low quality. These included SNPS

#2) Individual missingness SNP plot 





#1) Plot MAF ####
#Get frequency file
#T1T2_remap.frq

setwd('/home/mari/GWAS_22/new_gwas/QC')
#Generate frequency file
system("plink bfile P1T1 --freq")

maf_freq <- read.table("all_enteric.frq", header =TRUE, as.is=T)





setwd('/home/mari/GWAS_22/new_gwas/Post-imp/Merged/Assoc')
#PLINKQC Set up
indir<-"~/GWAS_22/new_gwas/Post-imp/Merged"
qcdir<-"~/GWAS_22/new_gwas/Post-imp/Merged"
name<-"VAST_concat"
path2plink <- "/home/mari/bin/plink"








#Post-impute plots ####
#Set up wd ##
setwd('/home/mari/GWAS_22/new_gwas/Post-imp/Merged/Assoc')
#PLINKQC Set up
indir<-"~/GWAS_22/new_gwas/Post-imp/Merged"
qcdir<-"~/GWAS_22/new_gwas/Post-imp/Merged"
name<-"VAST_concat"
path2plink <- "/home/mari/bin/plink"


maf_freq <- read.table("all_enteric.frq", header =TRUE, as.is=T)
pdf("MAF_distribution.pdf")

hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")

str(maf_freq)

b <- maf_freq %>% filter(MAF > 0.05) %>%
  ggplot(aes(x = MAF, y = ..count..)) + 
  geom_histogram(bins = 24, fill = "purple", colour = "grey") +ggtitle("Filtered SNPs") +ylab("Count\n") + xlab("Minor allele frequency") +theme_minimal()

a <- maf_freq %>%
  ggplot(aes(x = MAF, y = ..count..)) + 
  geom_histogram(fill = "purple", colour = "grey") + geom_vline(xintercept = 0.05, linetype="dashed", 
                                                                color = "blue", size=1) +ggtitle("Total Imputed SNPs") +ylab("Count\n") + xlab("Minor allele frequency") +theme_minimal()

a

a + b + plot_annotation(tag_levels = "A")

