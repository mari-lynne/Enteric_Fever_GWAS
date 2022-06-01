#Update sex ####
#Need to use pre-imp data as we filtered autosomes

library(dplyr)
library(data.table)
library(stringr)
library(tidylog)
library(dplyr)
library(tidyr)
library(RColorBrewer)

#VAST covar file, need new file with updated cat_IDs
covar <- "~/GWAS_22/new_gwas/meta/final_metadata/all_covar.txt"

#covar <- "~/GWAS_22/new_gwas/meta/covar.txt"
prefix <- "VAST_PATCH" #Plink prefix of data

setwd(paste("~/GWAS_22/new_gwas/",prefix,"/sex_check", sep=""))

#Check covar file ####
data <- fread(paste(prefix,".fam", sep =""))

#This has IDs FID = IID hence how they got joined at somepoint
covar <- fread(covar)
#1 is male 2 is female

#Make update sex file ####

sex <- select(covar,IID,Sex)
colnames(sex) <- c("V2","SEX")

#add in fam column 
fam <- data[,c(1:2)]
fam$V2 <- as.character(fam$V2)
sex <- left_join(fam, sex, by = "V2") #matched 105 rows
sex<-unique(sex)
  
#Reorer so is FID, IID, Sex (as requested by plink)
#sex <- sex[,c(1,2,4)]
colnames(sex) <- c("FID", "IID", "Sex")

#Recode sex
#sex M=1, F=2
sex$Sex <- ifelse(sex$Sex == "M", 1, ifelse(sex$Sex == "F", 2, NA))

#duplicate samples, check properly later
#sex <- sex[-c(19,45),]

sex <- drop_na(sex, Sex)

#write.table(sex, file = "~/GWAS_22/new_gwas/T1T2/sex_check/update_sex.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#Save update file
write.table(sex, file = "update_sex.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#Exlude xy dodgy snps in R

#More file cleaning
#Load SNP info
bim <-fread(paste(prefix,".bim", sep =""))

#Filter SNPs in autosomal chromosomes")
sum <- table(bim$V1) #weird xy category with 22415 snps
# 17599  1320   712
sum

#T1T2 P1Tyger
#xy <- bim %>% filter(V1 == "23" | V1 == "24" | V1 == "25")

#VAST
xy <- bim %>% filter(V1 == "X" | V1 == "XY" | V1 == "Y")

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
write.table(file="exclude_snps_xy.txt", exclude, sep = "\t", quote=F, row.names = F, col.names = F)

#Run sex_check.sh


#Plot ####
sexcheck <- fread(paste(prefix,"_redo.sexcheck", sep=""))

sexcheck$SNPSEX <- as.factor(sexcheck$SNPSEX)
levels(sexcheck$SNPSEX) <- c("M","F")

ggplot(sexcheck, aes(x=F, fill =SNPSEX)) + geom_histogram(bins=25, color="black") + labs(title = paste(prefix,"Sex Check"), x = " \nX-Chr Homozygosity (F-statistic)") +scale_fill_brewer(palette="Dark2")
#Males have only one copy of snp therefore homozygous
#Females may have two copies of an X SNP therefore are heterozygous



#Write fail list ####

fail <- filter(sexcheck, STATUS == "PROBLEM", PEDSEX != "0")
write.table(file = paste(prefix, "_failSex.txt", sep = ""),fail, sep = "\t", quote=F, row.names = F, col.names = F)

ggplot(fail, aes(x=F, fill =SNPSEX)) + geom_histogram(color="black") + labs(title = paste(prefix,"Sex Check"), x = " \nX-Chr Homozygosity (F-statistic)") +scale_fill_brewer(palette="Dark2")


#impute sex for missing phenotypes

impute <-  filter(sexcheck, PEDSEX == "0")
write.table(file = paste(prefix, "_impSex.txt", sep = ""),impute, sep = "\t", quote=F, row.names = F, col.names = F)








