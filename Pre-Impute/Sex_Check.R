#Update sex ####
#Need to use pre-imp data as we filtered autosomes

library(dplyr)
library(data.table)
library(stringr)
library(tidylog)

setwd("~/GWAS_22/new_gwas/T1T2/sex_check")
data <- "T1T2.fam"
covar <- "~/GWAS_22/new_gwas/meta/covar.txt"

#Check covar file ####
data <- fread(data) #this has IDs FID = IID hence how they got joined at somepoint
covar <- fread(covar) #1 is male 2 is female

#Make update sex file ####
#fam file has .cel extension in ID name
sex <- select(covar,GENOTYPING_ID,Sex)
colnames(sex) <- c("V2","SEX")

#add in fam column 
fam <- data[,c(1:2)]
fam$V2 <- as.character(fam$V2)
sex <- left_join(fam, sex, by = "V2") #matched 105 rows
  sex<-unique(sex)
  
#Reorer so is FID, IID, Sex (as requested by plink)
sex <- sex[,c(1,2,4)]
colnames(sex) <- c("FID", "IID", "Sex")

#Recode sex
#sex M=1, F=2
sex$Sex <- ifelse(sex$Sex == "M", 1, ifelse(sex$Sex == "F", 2, NA))

#duplicate samples, check properly later
sex <- sex[-c(19,45),]

sex <- drop_na(sex, Sex)

write.table(sex, file = "~/GWAS_22/new_gwas/T1T2/sex_check/update_sex.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#Try with plink2 first
#--update-sex expects a file with FIDs and IIDs in the first two columns, and sex information (1 or M = male, 2 or F = female, 0 = missing


#Then run sex-check with y counts

hh <- fread("rsv_redo.hh", header = F)
#Family ID, Within-family ID,Variant ID
bim <- rename(bim, SNP = V2)
hh <- rename(hh, SNP = V3)
chr_check <- left_join(hh, bim, by = "SNP")
table(chr_check$V1.y)
#Chr 23 = X chromosome where all the hetezygous snps are, isn't this what we would expect?


#More file cleaning
#Load SNP info
bim <-fread("T1T2.bim")

#Filter SNPs in autosomal chromosomes")
table(bim$V1)

xy <- bim %>% filter(V1 == "23" | V1 == "24" | V1 == "25")

#Remove duplicates 
duplicated_snps <- xy[duplicated(paste(xy$V1,xy$V4)),]
#remove snps where A > T, C > G
#PLINK can't tell what strand to merge/read these SNPs from
atcg_snps <-  xy[((xy$V5 == "A") & (xy$V6 == "T")) |
                   ((xy$V5 == "T") & (xy$V6 == "A")) |
                   ((xy$V5 == "C") & (xy$V6 == "G")) |
                   ((xy$V5 == "G") & (xy$V6 == "C"))]
#Exclusion list
exclude <- rbind(duplicated_snps, atcg_snps)
exclude <- exclude$V2
write.table(file="sex_check/exclude_snps_xy.txt", exclude, sep = "\t", quote=F, row.names = F, col.names = F)


#Results #####
library(ggplot2)

sexcheck <- fread("T1T2_redo.sexcheck")
ggplot(sexcheck, aes(x=F)) + geom_histogram(bins =25) + labs(title = "Plink Sex Check -hh", x = " \nX-Chr Homozygosity (F-statistic)") #600W, H=385
#below 0.2 = Female, above 0.8 = male (looks right from the data)

#Heterozygosity check #####################
#145 sample mis-sexed seems like a lot!
#136021 het. haploid genotypes present (see rsv_sex_clean.hh ); many commands treat these as missing.

hh <- fread("rsv_sex_clean.hh", header = F)
#Family ID, Within-family ID,Variant ID

bim <- rename(bim, SNP = V2)
hh <- rename(hh, SNP = V3)

chr_check <- left_join(hh, bim, by = "SNP")
table(chr_check$V1.y)
