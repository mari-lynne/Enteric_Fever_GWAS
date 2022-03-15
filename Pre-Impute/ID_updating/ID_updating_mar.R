#ID sorting
library(dplyr)
library(data.table)
library(stringi)
library(stringr)
library(tidylog)

#Aims ####

#1)Need to double check all IDs
#2)Need to update pheno/covar files to include cat_IDs from the fam files
#3)Need to update all files to include relevant PATCH rechallenge data 

# 1) ID checking:

setwd('~/GWAS_22/new_gwas/meta/')

#check sample IDs in fam file

fam <- fread("enteric.fam")
#need to remove first number and delimiter (this was added during the vcf>plink process)
fam <- fam %>% mutate(Split_ID = str_extract(V2, "(?<=_).*"))
#look behind ?<= underscore _ #Do for any characters . #As many times *

#Update fam file V2 col
fam <- fam %>% select(V1, Split_ID, V3, V4, V5, V6)
write.table(file = "~/GWAS_22/new_gwas/Post-imp/Merged/Assoc/enteric.fam", fam, row.names = F, col.names = F, quote=F, sep="\t")

#Reload fam file

fam <- fread("~/GWAS_22/new_gwas/Post-imp/Merged/Assoc/enteric.fam")

VP <- str_split_fixed(fam$V2, "-", n = 2)

?str_split

#Check pheno/covar files

covar <- fread("~/GWAS_22/new_gwas/meta/covar.txt")
pheno <- fread("~/GWAS_22/new_gwas/meta/pheno.txt")
#Note these include now added in data for geno samples 41, 103, and 74 
#Samples 121, 73 no metadata or info that I could find, poss withdrew

#2) ID matching

#Update IDs of pheno/covar files to use the IDs from fam files (cat_ID)
ID <- read.csv("~/GWAS_22/new_gwas/meta/new_genoID_march.csv")
ID <- distinct(ID) #468

#There are duplicate catID in ID files
#Only VAST_PATCH samples had an extra pos ID added - new_VPIDs2.txt
#Updated genoID_march file with correct cat IDs

#Pheno and Covar file IDs need to be updated with catID column 
#cat IDs are used in fam file
#Joining info is in the newIDs file
covar <- fread("~/GWAS_22/new_gwas/meta/covar.txt")
pheno <- fread("~/GWAS_22/new_gwas/meta/pheno.txt")

#need a study column in IDs then can join to covars
#without study we would run into dup problems
#cp columns over into pheno file (once ordered)

#look ahead as we are looking ahead FROM the pattern searched for (-)
ID <- ID %>% mutate(Study = str_extract(META_ID, "\\w*(?=-)"))

#just need genotyping ID CAT ID and study
join <- dplyr::select(ID, GENOTYPING_ID, cat_IID, Study)

ID$Study<- as.factor(ID$Study)
covar$Study <- as.factor(covar$Study)

covar2 <- left_join(covar, ID, by = c("Study", "GENOTYPING_ID"))

#Matched 315 rows 3 extra rows in covar bc samples were removed due to fail checks
#Tidy table

covar2 <- select(covar2, cat_IID, 2:17)

#Bind pheno to covar
pheno2 <- bind_cols(covar2, pheno)
pheno2 <- select(pheno2, cat_IID, Diagnosed)

#check
test <- left_join(pheno2, covar2, by = "cat_IID") 
test <- distinct(test) 



#319 rows remaining #extra patch sample with no cat ID

#3)Update to include PATCH re-challenge IDs ####

#Plan:
#link the people who were rechallenged in PATCH with their previous IDs from the studies they were originally challenged in.

#Linking info found
#Join genoIDs march by Participant ID, new column prev ID has the info
#Will then need to duplicate entries in PLINK or VCF files (including the geno info
#New dup geno ID will be prev_ID plus_2 to signify duplicated entry)
#Double check there aren't any duplicates between them
There are duplicates between T1 and P1, merge with genotyping ID table by study as well
#Update PLINK files
#Update genotyping ID to include these additional entries Prev_ID_2 for the rechallenge samples
#Extra pheno will be the PATCH diagnosis info and covars for these samples
#make sure covar has prev challenge added
#Tempted to make new var if the prev challenge was homo or heterologous
#so covar one is prev challenge = P or T
#prev_challenge type = homo or hetero

#Get prev_study genotyping IDs ####

#Need to go from prev_ID in patch file thru lab_IDs in ID file to get geno_IDs

#Load linking info 
link <- fread("~/GWAS_22/new_gwas/meta/patch_linking_IDs.csv")
link <- select(link, PARTICIPANT_ID, Prev_ID, Study)

#Use Prev_ID to match to LAB_ID in ID table
colnames(link) <- c("PARTICIPANT_ID", "LAB_ID", "Study")

link2 <- link %>% filter(!is.na(LAB_ID)) #rm 9 rows

ID <- ID %>% filter(!is.na(GENOTYPING_ID)) #316 rows so only matching to Genotyped samples

link2 <- dplyr::left_join(link, ID, by = c("Study", "LAB_ID"))
#Will have 2 diff participant IDs .x is their patch ID .y is the prev participant ID I think

#clean
link2 <- distinct(link2)
link2 <- link2 %>% filter(!is.na(GENOTYPING_ID)) #27 rows removed 34 extra samples :)

#rename cat_IID eventually to be cat_IID + PATCH

#Get PATCH pheno covar data ####
#Based on participant ID

meta <- fread("enteric_meta_march.csv")
meta <- filter(meta, Study == "PATCH")
#format as covar file can fill out/sort meta file later with all Patch from the diagnosis file 

link2 <- rename(link2, PARTICIPANT_ID =  PARTICIPANT_ID.x)
extra_meta <- left_join(link2, meta, by = "PARTICIPANT_ID")

#need to rename .x's
extra_meta <- extra_meta %>% rename_with(~str_remove(., '\\.x'))

#Rename cat_IID to include PATCH
extra_meta$cat_IID <- str_c(extra_meta$cat_IID, extra_meta$Study.y, sep = "_")

#53 matched rows 1 row only in extra that's unmatched
#Make RC covar and pheno file ####
colnames <- colnames(covar2)
colnames

patch_covar <- select(extra_meta, colnames)
patch_pheno <- select(extra_meta, cat_IID, Diagnosed)

write.table(patch_covar, file = "patch_rc_covar.txt", row.names = F, quote=F, sep="\t")
write.table(patch_covar, file = "patch_rc_pheno.txt", row.names = F, quote=F, sep="\t")

#PLINK -  make keep list ####

#based on the CAT IDs (use pre str_c version)
extra_meta2 <- left_join(link2, meta, by = "PARTICIPANT_ID")
old <- select(extra_meta2, cat_IID, META_ID)
new <- select(extra_meta, cat_IID, META_ID)

keep <- select(old, cat_IID)
update_name <- left_join(old, new, by = "META_ID")
update_name = select(update_name, - META_ID)

#write files
write.table(keep, file = "keep_IDs.txt", row.names = F, quote=F, sep="\t")
write.table(update_name, file = "update_name.txt", row.names = F, col.names = F, quote=F, sep="\t")

#filter for these samples in PLINK save as an extra dataset
setwd("~/GWAS_22/new_gwas/Post-imp/Merged/Assoc")

system("plink2 --bfile enteric --keep keep_IDs.txt --make-bed --out rechallenge")
#33 samples (1 must have been a dodge one)

#Update sample names 
system("plink2 --bfile rechallenge --update-ids update_name.txt --make-bed --out rechallenge_rn")

#Concatenation ####

system("plink --bfile enteric --bmerge rechallenge_rn --out all_enteric") 
#348 People :)
#6260 same-position warnings 9,970,728 varients

#Concat pheno and covar files 

all_pheno <- rbind(pheno2, patch_pheno)
#352 varients

all_pheno <- filter(all_pheno, !is.na(cat_IID))

all_covar <- rbind(covar2, patch_covar)
all_covar <- filter(all_covar, !is.na(cat_IID))

fam <- fread("all_enteric.fam")

#488 and 573 are missing a pheno, had looked previously and they are nowhere to be found :(
#Also update sex 

write.table(all_pheno, file = "all_pheno.txt", row.names = F, quote=F, sep="\t")
write.table(all_covar, file = "all_covar.txt", row.names = F, quote=F, sep="\t")

#Also update sex
#Have already removed old chromosomes 
#Do by old cat ID and covar

update_sex <- select(covar, GENOTYPING_ID, Sex)
('1' = male, '2' = female, '0' = unknown)

update_sex$Sex <- update_sex$Sex <- ifelse(update_sex$Sex == "M", 1, ifelse(update_sex$Sex == "F", 2, 0))

write.table(update_sex, file = "update_sex.txt", row.names = F)
