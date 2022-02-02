#Date 28/01/22
setwd("~/GWAS_22/Enteric_Fever_GWAS/Pre-Impute/ID_updating")

#make new branch git 
git checkout -b ID_updates
#update remote
git push --set-upstream origin ID_updates

#Background ####
#Previously when merging T1T2P1 with VAST_PATCH there were some ID issues in that the genetics(blinded) IDs of some of the VP samples were unfortunately the same as the T1T2 study.
#I tried to work around it by modifying the files to add _26 to the end of duplicate IDs in VAST/PATCH...
#But because they were in associated files etc. this was more difficult the you would think/might be confusing for future analysis

#Proposed solution update IIDs from the start to avoid confusion/changing midway through
#new VAST Patch IID consists of 'OLD_IID-FID' concatenated e.g  73-4

setwd("~/GWAS_22/Enteric_Fever_GWAS/Pre-Impute/ID_updating")

#Steps ####
#use update IDs flag in plink with a .txt file with the old and new IDs respectively, see link:
https://www.cog-genomics.org/plink/1.9/data#update_indiv

#Get IDs from fam files (in all_enteric/clean_data/raw_data/VAST_PATCH)
#Export as .csv
#Make update files as below 

plink2 --bfile --update-IDs --make-bed --out newfile

--update-ids expects input with the following four fields:
  
1) Old family ID
2) Old within-family ID
3) New family ID
4) New within-family ID

setwd("~/GWAS_22/new_gwas/VAST_PATCH") 
getwd()
cd ~/GWAS_22/new_gwas/VAST_PATCH
library(data.table)
library(dplyr)
library(stringr)
ids <-fread("new_vp_IDs.csv")

ids$IID_new <- str_c(ids$IID_new, "-", ids$FID)

write.table(ids, file = "new_vp_IDs2.txt", row.names = F, col.names = F, quote = F, sep = "\t")

#test in Plink
system("plink2 --bfile VAST_PATCH --update-ids new_vp_IDs2.txt --make-bed --out VAST_PATCH_2")

#check new fam file
check <- fread("VAST_PATCH_2.fam")

#success IDs updated with concatenated IDs
#save new genotyping IDs to metadata ID folder 

write.csv(check[,2], file = "new_VAST_IID.csv")

setwd("~/GWAS_22/new_gwas/meta")
IID <- check[,2]

enteric <- read.csv("enteric_IDs_Feb22.csv")
?str_extract
IID <- IID %>% mutate(Split_ID = str_extract(V2, "\\d*(?=-)"))

colnames(IID) <- c("cat_IID", "GENOTYPING_ID")
new <- left_join(enteric, IID, by = "GENOTYPING_ID")

new <- select(new, -Geno)
write.csv(new, file = "new_genoIDs.csv")
