#Background ####

### Have QCd data for VAST_PATCH, (P1,Tyger,T1,T2)- these were merged
#rename P1TygerT1T2 to P1T1 for short
#Imputed snps for VAST_PATCH and P1T1 study cohorts separately
#downloaded chr data, unzipped and concatanated chrs
#rezipped then merged x2 .vcf.gz files using bcftools merge
#merged data saved as all_enteric
#convert to plink format

#Thresholds for pre-imputation were not as stringent as final GWAS QC
#Redo so to include snps just for assoc studies
#Visualise using PlinkQC
#run thresholds to remove poor qual/rare snps



#Plink_qc.R ####
install.packages("devtools")
library(devtools)
install_github("meyer-lab-cshl/plinkQC")

library(plinkQC)

cd 
#set directories
indir<-"/home/mari/GWAS/all_enteric/all_data"
qcdir<-"/home/mari/GWAS/all_enteric/all_data"
name<-"all_enteric"
path2plink <- "/home/mari/GWAS/all_enteric/all_data/plink"

#filter snps with poor call rate/rare ####

fail_maf <- check_maf(indir=indir, qcdir=qcdir, name=name, interactive=TRUE,
                      path2plink=path2plink, showPlinkOutput=FALSE)

fail_maf

?check_maf

failed

#big grouped plot of MAF, HQW and MAF distribution
fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name,
                            path2plink=path2plink, mafTh = 0.05,
                            verbose=TRUE, interactive=TRUE,
                            showPlinkOutput=FALSE)
fail_markers$p_markerQC



#Filter individuals ####

fail_het_imiss <- check_het_and_miss(indir=indir, qcdir=qcdir, name=name,
                                     interactive=TRUE, path2plink=path2plink)

#relatedness, dont know if this one works
exclude_relatedness <- check_relatedness(indir=indir, qcdir=qcdir, name=name,
                                         interactive=TRUE,
                                         path2plink=path2plink)

#PLINK filtering ####

./plink --bfile all_enteric --geno 0.01 --maf 0.05 --hwe 0.000001 --make-bed --out all_enteric_QC

#plot MAF ####

maf_freq <- read.table("all_enteric.frq", header =TRUE, as.is=T)
pdf("MAF_distribution.pdf")

hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")

str(maf_freq)

library(ggplot2)
b <- maf_freq %>% filter(MAF > 0.05) %>%
  ggplot(aes(x = MAF, y = ..count..)) + 
  geom_histogram(bins = 24, fill = "purple", colour = "grey") +ggtitle("Filtered SNPs") +ylab("Count\n") + xlab("Minor allele frequency") +theme_minimal()

a <- maf_freq %>%
  ggplot(aes(x = MAF, y = ..count..)) + 
  geom_histogram(fill = "purple", colour = "grey") + geom_vline(xintercept = 0.05, linetype="dashed", 
                                                                color = "blue", size=1) +ggtitle("Total Imputed SNPs") +ylab("Count\n") + xlab("Minor allele frequency") +theme_minimal()

a
library(patchwork)
#install.packages("patchwork")

a + b + plot_annotation(tag_levels = "A")

#dev.off()

#filter individuals
#
./plink --bfile all_enteric_QC --mind 0.01 --make-bed --out all_enteric_QC

#sort phenotype file (I believe it has all the genetic IDs so make file just with genetic ID and outcome, the other can have the covar information, rename to IID :))
#ent

#check ancestry/genomic inflation factor pca

#Pca redo ######

./plink2 --bfile all_enteric_QC --pheno pheno_rn.txt --pheno-name Diagnosed --covar covar_rn.txt --covar-name Sex, Challenge, Dose, Study, Vaccine, Age --pca

#recode categorical pheno to none 

./plink2 --bfile all_enteric_QC --pheno pheno_rn2.txt --pheno-name Diagnosed --covar covar_rn2.txt --covar-name Sex, Challenge, Dose, Vaccine, Age --pca


library(data.table)

pca <- fread("./plink2.eigenvec")
eigenval <- scan("./plink2.eigenval")

#Make table:
  Sex <- rep(NA, length(pca$IID))
Sex[grep("1", pca$IID)] <- "Male"
Sex[grep("2", pca$IID)] <- "Female"
# location
Challenge <- rep(NA, length(pca$IID))
Challenge[grep("1", pca$IID)] <- "Paratyphoid"
Challenge[grep("2", pca$IID)] <- "Typhoid"
Challenge[grep("3", pca$IID)] <- "Toxin-Neg"


# outcome
Outcome <- rep(NA, length(pca$IID))
Outcome[grep("2", pca$IID)] <- "Typhoid"
Outcome[grep("1", pca$IID)] <- "No Typhoid"


# remake data.frame
pca2 <- as_tibble(data.frame(pca, Outcome, Challenge, Sex))

# Convert to percentage variance explained
pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

library(ggplot2)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

#use PCA 1-5 as covars

# plot pca
b <- pca2 %>% filter(Sex != "NA") %>% ggplot(aes(PC1, PC2, col = Outcome, shape = Challenge)) + geom_point(size = 1.4)
b
b <- b + scale_colour_manual(values = c("#40d5c4", "#d54051")) + ylim(-0.1,0.1)

b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

#divergent ancestry - pre correction

#merge PCA with covar file
#remove hashtag from IID
covar <- fread("covar_update.txt")


pca3 <- pca[,1:7]


covar_pca <- left_join(covar, pca, by = "IID")

covar$IID <-as.character(covar$IID)

library(tidylog)

write.table(covar_pca, file = "covar_pca.txt", row.names = FALSE, col.names = TRUE, sep="\t", quote=F)

#run plink

#update NAs ####

covar_pca <- str_replace_all()

covar_pca[] <- lapply(covar_pca, gsub, pattern = "FALSE", replacement = "NA", fixed = TRUE)

#recode challenge, vaccine, sex ####

str(covar$Vaccine)
#NONE = 0
#Ty21a = 1
#M01ZH09= 2
#Vi-PS = 3
#Vi-TT= 4

covar_pca$Vaccine <- ifelse(covar_pca$Vaccine == "None",0,
                                 ifelse(covar_pca$Vaccine == "Ty21a",2,
                                        ifelse(covar_pca$Vaccine == "M01ZH09", 2,
                                               ifelse(covar_pca$Vaccine == "Vi-PS", 3,
                                                      ifelse(covar_pca$Vaccine == "Vi-TT", 4,NA)))))

#Recode strain

covar_pca$Challenge <- as.factor(covar_pca$Challenge)

levels(covar_pca$Challenge)

covar_pca <- filter(covar_pca, Challenge != "C")

#C, P = 3, T = 1, TN = 2
levels(covar_pca$Challenge) <- c("0", "3", "1", "2")

#sex M=1, F=2
covar_pca$Sex <- ifelse(covar_pca$Sex == "M", 1, ifelse(covar_pca$Sex == "F", 2, NA))

#update table
write.table(covar_pca, file = "covar_pca.txt", row.names = FALSE, col.names = TRUE, sep="\t", quote=F)

#covar_pca still
Error: 'PC1' entry on line 4 of covar_pca.txt is categorical, while earlier
entries are not.
(Case/control and quantitative phenotypes must all be numeric/'NA'.
  Categorical phenotypes cannot be 'NA'--use e.g. 'NONE' to represent missing
  categorical values instead--or start with a number.)
End time: Sat Nov 27 01:20:52 2021


#replace to NONE - note i am concenered it's being interpreted as categorical

covar_pca[,19:28] <- lapply(covar_pca[,19:28], gsub, pattern = "NA", replacement = "NONE", fixed = TRUE)

glimpse(covar_pca)

covar_pca <- dplyr::select(covar_pca, -dup_check)

write.table(covar_pca, file = "covar_pca2.txt", row.names = FALSE, col.names = T, sep = "\t", quote=F)

#other ethnic group tabbed into the next col was the problem xx
#fix in R #remove col for now #add _delimiter

covar_pca <- select(covar_pca, -Ethnicity)

?as.numeric

  covar_pca$PC1 <- as.numeric(covar_pca$PC1)

#still no joy, convert to plink bfile for now and use 1.9?
  
pheno <- fread("pheno_update.txt")

pheno <- mutate(pheno, FID = rep(0))

pheno = select(pheno, FID, IID, Diagnosed)

write.table(pheno, file = "pheno2.txt", row.names = FALSE, col.names = T, sep = "\t", quote=F)

covar_pca <- mutate(covar_pca, FID = rep(0))

covar_pca <- covar_pca %>%
  select(FID, everything())

write.table(covar_pca, file = "covar2.txt", row.names = FALSE, col.names = T, sep = "\t", quote=F)

#Error: Duplicate sample ID in --covar file. ####

covar2 <- covar_pca %>% mutate(dup_check = stri_duplicated(covar_pca$IID))

stri_duplicated(covar_pca$IID)

#VAST IDs massive issue, go back to ID_updating script ####


#LD #####
./plink2 --bfile all_enteric_QC --indep-pairwise 50 5 0.2 --out all_enteric_LD


./plink2 --bfile all_enteric_QC --extract all_enteric_LD.prune.in --make-bed --out  all_enteric_LD 

#redo pca

./plink2 --bfile all_enteric_LD --pheno pheno_rn2.txt --pheno-name Diagnosed --covar covar_rn2.txt --covar-name Sex, Challenge, Dose, Vaccine, Age --pca

#run GWAS of enteric fever exclude Vi-PS/TCV participants ####
./plink2 --bfile all_enteric_LD --pheno pheno_rn2.txt --pheno-name Diagnosed --covar covar_pca_LD.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --glm --adjust --out enteric_assoc

./plink2 --bfile all_enteric_LD --pheno pheno_rn2.txt --pheno-name Diagnosed --covar cplink2.eigenvec --covar-name PC1, PC2, PC3, PC4 --glm --adjust --out enteric_assoc

./plink2 --bfile all_enteric_LD --pheno pheno_rn2.txt --pheno-name Diagnosed --covar covar_rn2.txt --covar-variance-standardize --covar-name Age, Sex, Challenge, Vaccine, Dose --glm --adjust --out enteric_assoc2

./plink2 --bfile all_enteric_LD --chr 1 --from-bp 161284166 --to-bp 161697933 --pheno pheno_rn2.txt --pheno-name Diagnosed --covar covar_rn2.txt --covar-variance-standardize --covar-name Age, Sex, Challenge, Vaccine, Dose --glm --adjust --out enteric_assoc2

#./plink2 --bfile FCGRv2merge --exclude  FCGRv2_LD.prune.out --pheno fullmerge2.txt --pheno-name Outcome --covar plink2.eigenvec --covar-name PC1, PC2 --glm --adjust --out pcatestld

filter(covar_pca, is.na(PC1))

#Visualise with qqman

Needs fields   $ SNP: chr  "rs1" "rs2" "rs3" "rs4" ...
$ CHR: int  1 1 1 1 1 1 1 1 1 1 ...
$ BP : int  1 2 3 4 5 6 7 8 9 10 ...
$ P

install.packages("qqman")
library(qqman)

gwasResults <- fread("enteric_assoc")
gwasResults <- fread("enteric_assoc2")

colnames(gwasResults) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1","FIRTH?","TEST","OBS_CT","OR",        
                         "LOG(OR)_SE","Z_STAT","P","ERRCODE")

gwasResults$TEST <- as.factor(gwasResults$TEST)

gwasResults <- filter(gwasResults, TEST == "ADD")

gwasResults <- dplyr::select(gwasResults, SNP, CHR, BP, P)

manhattan(gwasResults, ylim = c(0, 8), cex = 0.6, cex.axis = 0.9, col = c("#34b368","#b3347"))

manhattan(gwasResults, main = "Manhattan Plot", ylim = c(0, 8), cex = 0.6, cex.axis = 0.9, 
          col = c("#34b368","#b3347"), suggestiveline = F, genomewideline = F, chrlabs = c(1:20, 
                                                                                           "P", "Q"))

manhattan(gwasResults, main = "Manhattan Plot", ylim = c(0, 8), cex = 0.6, cex.axis = 0.9, 
          col = c("#34b368", "orange3"), suggestiveline = F, genomewideline = F, chrlabs = c(1:20, 
                                                                                           "P", "Q"))

?manhattan #F

#Check for FcR associations and pray


./plink2 --bfile all_enteric_LD --chr 1 --from-bp 161505430 --to-bp 161678022 --pheno pheno_rn2.txt --pheno-name Diagnosed --covar covar_rn2.txt --covar-variance-standardize --covar-name Age, Sex, Challenge, Vaccine, Dose --logistic --freq case-control --out Fc_assoc

library(data.table)
library(dplyr)
setwd('/home/mari/GWAS/all_enteric/all_data')

#Plink script

#gwasResults <- fread("ennteric_assoc_glm_Fcy")
gwasResults <- fread("Fc_assoc.Diagnosed.glm.logistic.hybrid")


str(gwasResults$TEST)

colnames(gwasResults) <- c("CHR", "BP", "ID", "REF", "ALT", "A1","FIRTH?","TEST","OBS_CT","OR",        
                           "LOG(OR)_SE","Z_STAT","P","ERRCODE")

gwasResults$TEST <- as.factor(gwasResults$TEST)

gwasResults <- filter(gwasResults, TEST == "ADD")

gwasResults <- dplyr::select(gwasResults, ID, REF, ALT, CHR, BP, OR, P)

topsnps <- filter(gwasResults, P <0.05)

write.csv(topsnps, file = "topfcsnps_2.csv")

#read in counts table 

counts <- fread(Fc_assoc_counts)

table <- left_join(topsnps, Fc_assoc_counts, by = "ID")

View(table)

./plink --bfile all_enteric_LD --chr 1 --from-bp 161505430 --to-bp 161678022 --pheno pheno_rn2.txt --pheno-name Diagnosed --covar covar_rn2.txt --covar-name Age, Sex, Challenge, Vaccine, Dose --allow-no-sex --logistic --freq counts --out Fc_assoc2
#get frequecy counts using old plink
#freq rate using freq case control

frq <- select(Fc_frq, SNP, MAF_A, MAF_U)

names(frq) <- c("ID", "MAF_A", "MAF_U")

table <- left_join(table, frq, by = "ID")

table2 <- table %>% dplyr::select(-c(REF.x, ALT.x, CHR, OBS_CT))

77/0.4*100

table2 <- table %>% mutate(count_alt = (ALT_CTS/MAF_U)*100)

                 #TFCAR gene 54874231..54874231
#./plink2 --bfile all_enteric_LD --chr 19 --from-bp 54874231 --to-bp 54874231 --pheno pheno_rn2.txt --pheno-name Diagnosed --covar covar_rn2.txt --covar-variance-standardize --covar-name Age, Sex, Challenge, Vaccine, Dose --logistic --freq case-control --out Fcar

