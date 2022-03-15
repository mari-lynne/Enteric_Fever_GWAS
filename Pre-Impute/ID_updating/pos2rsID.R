#Update rs_IDs ####

#Aims:
#Update Chr:position snps to rsID format

#Methods:
#From USCS database download dbsnp rs IDs and positions
#Reformat ton include rs ID and chrom:pos:ref:alt unfo
#Left join to enteric bim file based on pos Ids
#Replace pos Ids with rs IDs in bim file

setwd("~/GWAS_22/new_gwas/Post-imp/Merged/Assoc")

library(data.table)
library(dplyr)
library(stringr)
library(tidylog)

#Download and check files ####

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz 
#dbsnp 153 is latest build sort for common snps however there is no ftp link yet
#Online version will crash so use 151 (2017) for now

rs_full <- fread("snp151Common.txt")
bim_check <- fread("all_enteric_QC.bim")

head(bim_check)
chr1:710225:T:A dbsnp = rs185127847
chr1rs <- filter(rs_full, V2 == "chr1")
filter(rs_full, V5 == "rs185127847")
#matches :)


#Pre-cleaning :) ####

#Filter autosomes
rs_new <- rs_full %>% filter(V2 != c("chrX", "chrY", "chr23", "chr24"))

#Filter includes Topmed
rs_new <- rs_new %>% filter(str_detect(V21, 'TOPMED'))

#Filter inconsistent alleles
rs_new <- rs_new %>% filter(V19 != "InconsistentAlleles")


#update-name ####
#requires new ID old ID
#Update the bim file in R based of the dbsnp table

#Reformat dbsnp dataframe ####
#from rs full make the same ID column
#cat chr:pos:ref:alt
#cols chr:pos = V2:V4:
#col ref:alt = V10 then will need to sub / for :
#^=old ID, then keep the V5 col as new ID

rs_new <- rs_new %>% mutate(chrpos = str_c(V2, V4, sep = ":"))
rs_new$V10 <- rs_new$V10 %>% str_replace_all("\\/", ":")

#make new col
rs_new <- rs_new %>% mutate(bim_ID = str_c(chrpos, V10, sep = ":"))

#Filter to just have bim_ID and rs ID cols

rs_IDs <- rs_new %>% select(bim_ID, V5)
colnames(rs_IDs) <- c("bim_ID", "rs_ID")
head(rs_IDs)

bim_check <- rename(bim_check, bim_ID = V2)

update_names <- left_join(bim_check, rs_IDs, by = "bim_ID")
#rows only in x    5,865,028
#> rows only in y  (12,045,095)
#> matched rows      2,592,203
check_snps <- anti_join(bim_check, rs_IDs, by = "bim_ID")
head(check_snps)
#Not a duplicate issue, probably to do with multi-allelic snps
#also this was just a list of common snps so they could be missing from database
#when checking a few of the non-matched snps they were all rare
#leave for now may need to redo with all dbsnps later
#Could I instead try grepl or str_detect
#Match if grepl contailns first part of the string for multi-allelics

update_names <- left_join(bim_check, rs_IDs, by = "bim_ID")

#now need to update ID column with rs_IDs
#replace NAs in rs_ID with data from bim_ID
update_names <- update_names %>%
  mutate(V2 = coalesce(rs_ID, bim_ID))

new_bim <- select(update_names, V1, V2, V3, V4, V5, V6)

#Write bim file ####
write.table(new_bim, file = "all_enteric_QC2.bim", row.names = F, quote=F, col.names=F, sep="\t")

