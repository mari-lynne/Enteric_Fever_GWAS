#Feb 21 22 ###############################################################
#Files are updated IDs and saved as new bed files
#QC each before testing merge

#Pre-impute tool QC ####

#Do per study 
#VAST ####

setwd("~/GWAS_22/new_gwas/VAST_PATCH") 
library(data.table)
library(stringr)
library(dplyr)
library(tidylog)

#Update RS ID's from the GSA
#update Global Screening array plink file with Snp coordinates using update name function (first column =bp, second column = rsID, in 37 format)
remap <- fread("~/GWAS_22/new_gwas/GSA-24v3-0_rsids.txt")

#old snp name col 1, new snp name col2
system("./plink --bfile VAST_PATCH/VAST_PATCH_2 --update-name GSA-24v3-0_rsids.txt --make-bed --out VAST_rename")

#read data library(data.table)
bim<-fread("VAST_rename.bim")

#Filter SNPs in autosomal chromosomes")
chrom <- c(seq(1,23),"X")
bim <- bim[bim$V1%in%chrom,]
#722145 snps left

#keep only SNPS containing rs IDs #
bim$snp <- rownames(bim)
rs_bim <- bim %>% dplyr::filter(grepl("rs", V2))

#Remove duplicates ####
bim.dups<-rs_bim[duplicated(paste(rs_bim$V1,rs_bim$V4)),]
#subset all SNP IDs that are in col one and col 4 #9663 removed
duplicated_snps<-bim.dups$V2

#remove snps where A > T, C > G as PLINK can't tell what strand to merge/read these SNPs from
atcg_snps <- rs_bim$V2[((bim$V5 == "A") & (bim$V6 == "T")) |
                         ((bim$V5 == "T") & (bim$V6 == "A")) |
                         ((bim$V5 == "C") & (bim$V6 == "G")) |
                         ((bim$V5 == "G") & (bim$V6 == "C"))]

#make exclusion SNP list ##
exclude <- c(duplicated_snps, atcg_snps)

write.table(file="exclude_snps_rm.txt", exclude, sep = "\t", quote=F, row.names = F, col.names = F)

#SNPs to keep ####
write.table(rs_bim$V2, paste0("ValidSNPs_rm.txt"), sep = "\t", quote = F, col.names = F, row.names = F)

#Update PLINK files ####
system("./plink2 --bfile VAST_rename --extract ValidSNPs_rm.txt --make-bed --out VAST_rename_rm")
system("./plink2 --bfile VAST_rename_rm --exclude exclude_snps_rm.txt --make-bed --out VAST_rename_rm")

#Liftover ###

#Make BED file:
# Make BED file ####
system("awk '{print "chr" $1, $4 -1, $4, $2 }' VAST_rename_rm.bim > VAST_rename_rm_lift")
#do in terminal

#Liftover go online to https://genome.ucsc.edu/cgi-bin/hgLiftOver
#submit file at the bottom of the page
#Successfully converted 661430 records, Conversion failed on 250 records

# ectract mapped variants snp IDs
system("awk '{print $4}' VAST_38.bed > VAST_38_rm.snps")
# ectract updated positions
system("awk '{print $4, $3}' VAST_38.bed > VAST_38_rm.pos")

#update plink file with SNPs that we could extract, and update the snp coordinates using update map function
system("./plink2 --bfile VAST_rename_rm --rm-dup --make-bed --out VAST_rename_rm_dup")
#exclude duplicates again from rm-dup list
system("./plink2 --bfile VAST_rename_rm --exclude VAST_rename_rm_dup.rmdup.mismatch --make-bed --out VAST_rename_rm_dup")
system("./plink2 --bfile VAST_rename_rm_dup --extract VAST_38_rm.snps --update-map VAST_38_rm.pos --make-bed --out VAST_remap")

#Error: --update-map variant ID 'rs113994127' appears multiple times in dataset
#Fix with extra rm-dup step using PLINK2

#Test merge #### (these have all been cleaned and updated to 38)

#Merge
system("./plink --bfile t1t2/T1T2_remap --bmerge VAST_remap.bed VAST_remap.bim VAST_remap.fam --make-bed --out merge_rm")

#Error: Identical A1 and A2 alleles on line 2693 of VAST_remap.bim.
remap <- fread("VAST_remap.bim")
remap[2693,]
#remove all snps with missing V5 and V6 info
missing <- remap %>% filter(V5 == "." & V6 == ".")
#save V2 as exclusion snpID 
write.table(file="missing_snps_rm.txt", missing$V2, sep = "\t", quote=F, row.names = F, col.names = F)
system("./plink2 --bfile VAST_remap --exclude missing_snps_rm.txt --make-bed --out VAST_remap")


#remove all files with rm in bash shell
#rm *rm*
#move all files to new VAST folder
#mkdir vast_patch
#mv *VAST* vast_patch/

#T1T2 ####
#attempted to merge T1 P1 before hand but got error so do separately first
#system("./plink --bfile T1T2 --bmerge P1Tyger.bed P1Tyger.bim P1Tyger.fam --make-bed --out T1T2") centimorgan position errors