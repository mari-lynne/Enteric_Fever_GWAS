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

#sort phenotype file (I believe it has all the genetic IDs so make file just with genetic ID and outcome, the other can have the covar information, rename to IID :))


#enteric.psam #split by _ FID, IID #is ordered by IID
#pheno #reorder by IID descending, then paste the psam FIDs to match

#check for duplicates

setwd("~/GWAS/all_enteric/clean_data/post-impute/Plink")


#cant modify order of psam file or number of samples
#so just have to remove duplicates from covar file for now
#which is in covar rn2

#so now just need to edit psam file to zero, try constant fIID #this just removed FID column
#remove first bit f string by_

#duplicates add a number to differentiate
#split IID first


library(stringr)
library(stringi)
library(dplyr)

#Update psam ####

test <-  enteric %>%
  dplyr::mutate(IID = str_extract(IID, "[^_-]+$"))

dups <- unique(test)
stri_duplicated(test$IID)

#make new column

test_2 <- test %>% mutate(dup_check = stri_duplicated(test$IID))
#highlights first duplicate entry

#remove duplicates
dups <- test_2 %>% filter(dup_check == TRUE)

#edit IID column

test_2$IID <- ifelse(test_2$dup_check == TRUE, 
                         paste(test_2$IID , "26",sep = ""),
                     paste(test_2$IID , "",sep = ""))

#rm dup check column to tidy
psam_update <-  unique(test_2)

#Update pheno ####
pheno <- fread("pheno.txt")

pheno2 <-  pheno %>%
  dplyr::mutate(IID = str_extract(GENOTYPING_ID, "[^_-]+$"))

#make new column

pheno2 <- pheno2 %>% mutate(dup_check = stri_duplicated(pheno2$IID))
#highlights first duplicate entry

#edit IID column

pheno2$IID <- ifelse(pheno2$dup_check == TRUE, 
                     paste(pheno2$IID , "26",sep = ""),
                     paste(pheno2$IID , "",sep = ""))

#rm dup check column to tidy
pheno_update <-  unique(pheno2)

#covar update ####

covar2 <-  covar %>%
  dplyr::mutate(IID = str_extract(IID, "[^_-]+$"))

#make new column

covar2 <- covar2 %>% mutate(dup_check = stri_duplicated(covar2$IID))
#highlights first duplicate entry


#edit IID column

covar2$IID <- ifelse(covar2$dup_check == TRUE, 
                     paste(covar2$IID , "26",sep = ""),
                     paste(covar2$IID , "",sep = ""))

#rm dup check column to tidy
covar_update <-  unique(covar2)

#test IDs
library(tidylog)

test <- left_join(psam_update, covar_update, by = "IID")

#9 unmatched rows covar >pheno
#psam > covar = 6 unmatched rows
#psam > covar = 10 unmatched rows

#anti  return all rows from x without a match in y
missing_covar <- anti_join(psam_update, covar_update, by = "IID")
missing_pheno <- anti_join(psam_update, pheno_update, by = "IID")
missing <- anti_join(pheno_update, covar_update, by = "IID")

#double_check recoding with paste 26
#highlights first instance, that is different in psam and pheno file
#psam file order is P1_Tyger, T1T2, VAST_PATCH
#the dupicates highligted will be in P1Tyger/T1

#add study code temporarily to pheno file?
#check pheno/covar order, maybe easiest to reorder by study then redo

covar$Study <- as.factor(covar$Study)

levels(covar$Study) #not order

#reorder by factor to match psam ####

#with(covar, covar[order])
#Reupdate covariate file

covar <- read.csv("~/GWAS/all_enteric/clean_data/post-impute/Plink/covar.txt", sep="")
View(covar)

u <- unique(covar$GENOTYPING_ID) #309 unique IDs

#add study col to pheno

#?dplyr::bind_cols

test <- bind_cols(pheno, covar) #match up by rows

test

category_order <-c("P1", "TYGER", "T1", "T2", "VAST", "PATCH")

covar <- covar %>% 
  arrange(factor(Study, levels = category_order))



#remake pheno and covar
?dplyr::select

pheno <- select(test, GENOTYPING_ID...1, Diagnosed)
covar <- test[,3:19]

?dplyr::rename_with

pheno <- rename(pheno, IID = GENOTYPING_ID...1)
covar <- rename(covar, IID = GENOTYPING_ID...3)

#redo ID updating pheno ####

pheno2 <-  pheno %>%
  dplyr::mutate(IID = str_extract(IID, "[^_-]+$"))

#make new column

pheno2 <- pheno2 %>% mutate(dup_check = stri_duplicated(pheno2$IID))
#highlights first duplicate entry

#edit IID column

pheno2$IID <- ifelse(pheno2$dup_check == TRUE, 
                     paste(pheno2$IID , "26",sep = ""),
                     paste(pheno2$IID , "",sep = ""))

pheno_update <-  unique(pheno2)

#redo covar ####

covar2 <-  covar %>%
  dplyr::mutate(IID = str_extract(IID, "[^_-]+$"))

#make new column

covar2 <- covar2 %>% mutate(dup_check = stri_duplicated(covar2$IID))
#highlights first duplicate entry
#edit IID column
covar2$IID <- ifelse(covar2$dup_check == TRUE, 
                     paste(covar2$IID , "26",sep = ""),
                     paste(covar2$IID , "",sep = ""))

covar_update <-  unique(covar2)

#test with merges again ####
#covar > pheno = 9 unmatched rows 
#psam > covar = 6 unmatched rows
#psam > covar = 10 unmatched rows

#anti  return all rows from x without a match in y
missing_covar <- anti_join(psam_update, covar_update, by = "IID")
missing_pheno <- anti_join(psam_update, pheno_update, by = "IID")
missing <- anti_join(pheno_update, covar_update, by = "IID")

#pheno update and covar update match
#However there are still 7 mismatches between the pheno/covar files and the psam files
#Is probably a 0 issue - suspicions were correct (due to regex)
#Add zeros in excel

write.table(pheno_update, file = "pheno_update.txt", row.names = F, quote = F, sep = "\t")
write.table(covar_update, file = "covar_update.txt", row.names = F, quote = F, sep = "\t")
write.table(missing_covar, file = "zeroIDs2update")

psam_update <- psam_update[,1:2]

#Recode covar and pheno variables ####

#Tried to run assoc model with updated files however still getting issues, now with categorisation
#Error: 'Dose' entry on line 128 of covar_update.txt is categorical, while earlier entries are not.
#(Case/control and quantitative phenotypes must all be numeric/'NA'.
#Categorical phenotypes cannot be 'NA'--use e.g. 'NONE' to represent missing #this is quite annoying, message plink ppl
#categorical values instead--or start with a number.


#diagnosed = 2, nTD = 1, NA = NONE

#Recode pheno file

pheno_update$Diagnosed <- ifelse(pheno_update$Diagnosed == "1",2,
                            ifelse(pheno_update$Diagnosed == "0",1,"NONE"))

#recode sex, female = 1, male = 2


write.table(psam_update, file = "enteric2.psam", row.names = F, quote = F, sep = "\t")


#Update original PSAM
#tab delimited and has #IID i think the hashtag was the missing token
#get rid of double quotes
