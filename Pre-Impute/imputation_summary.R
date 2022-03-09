# Imputation from the top ####

#New approach steps #####
#QC studies separately
#Liftover all files separately
#Run will Rayner pre-impute QC separately
#Merge output of Will QC vcf files
#Zip VCF files
#upload merged VCF files to TOPMED
#Post-imputation QC

#Feb 05 2022 ####
#QC Studies (test merge)
#Liftover
#Will Rayener pre-impute QC
#Merge all files 
#Prepare for TOPMED Imputation

#Original files:
library(dplyr)
library(data.table)

#Raw_data "~/GWAS_22/new_gwas/StudyName")

setwd("~/GWAS_22/new_gwas/T1T2")

#convert all to bed files ####
./plink --file P1Tyger --make-bed --out P1Tyger
plink2 --map T1T2 --ped T1T2 --make-bed --out T1T2
#./plink --file Rdurury --make-bed --out VAST_PATCH 

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

#read data library(data.table)
bim<-fread("T1T2.bim")

#Filter SNPs in autosomal chromosomes")
chrom <- c(seq(1,23),"X")
bim <- bim[bim$V1%in%chrom,]
#722005 snps left

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
system("./plink --bfile T1T2 --extract ValidSNPs_rm.txt --make-bed --out T1T2_rm")
system("./plink --bfile T1T2_rm --exclude exclude_snps_rm.txt --make-bed --out T1T2_rm")

#Liftover ###

#Make BED file:

# Make BED file ####
system("awk '{print "chr" $1, $4 -1, $4, $2 }' T1T2_rm.bim > T1T2_rm_lift")
#do in terminal

#Liftover go online to https://genome.ucsc.edu/cgi-bin/hgLiftOver
#submit file at the bottom of the page
#Successfully converted 661430 records, Conversion failed on 250 records

# ectract mapped variants snp IDs
system("awk '{print $4}' T1T2_38.bed > T1T2_38_rm.snps")
# ectract updated positions
system("awk '{print $4, $3}' T1T2_38.bed > T1T2_38_rm.pos")

#update plink file with SNPs that we could extract, and update the snp coordinates using update map function
system("./plink --bfile T1T2_rm --extract T1T2_38_rm.snps --update-map T1T2_38_rm.pos --make-bed --out T1T2_remap")

#Clean folders ####

#remove all files with rm
#rm *rm*
#move all files to new VAST folder
#mkdir t1t2
#mv *T1T2* t1t2/

#P1 Tyger updating ####

#read data library(data.table)
bim<-fread("P1Tyger.bim")

#Filter SNPs in autosomal chromosomes")
chrom <- c(seq(1,23),"X")
bim <- bim[bim$V1%in%chrom,]
#722005 snps left

#keep only SNPS containing rs IDs #
bim$snp <- rownames(bim)
rs_bim <- bim %>% dplyr::filter(grepl("rs", V2))

#SNP annotation an issue for global screening array, clean what I can then hopefully liftover should updata
#what are ilmnseq, ilmndup data?
#Try with file from biostarts, have contacted illumina 

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
system("./plink --bfile P1Tyger --extract ValidSNPs_rm.txt --make-bed --out P1Tyger_rm")
system("./plink --bfile P1Tyger_rm --exclude exclude_snps_rm.txt --make-bed --out P1Tyger_rm")

#Liftover ####

# Make BED file ###
system("awk '{print "chr" $1, $4 -1, $4, $2 }' P1Tyger_rm.bim > P1Tyger_rm_lift")
#do in terminal

#Liftover go online to https://genome.ucsc.edu/cgi-bin/hgLiftOver
#submit file at the bottom of the page
#Successfully converted 661430 records, Conversion failed on 250 records

# ectract mapped variants snp IDs
system("awk '{print $4}' P1Tyger_38.bed > P1Tyger_38_rm.snps")
# ectract updated positions
system("awk '{print $4, $3}' P1Tyger_38.bed > P1Tyger_38_rm.pos")

#update plink file with SNPs that we could extract, and update the snp coordinates using update map function
system("./plink --bfile P1Tyger_rm --extract P1Tyger_38_rm.snps --update-map P1Tyger_38_rm.pos --make-bed --out P1Tyger_remap")

#Clean files ####

#remove all files with rm
#rm *rm*
#move all files to new P1 folder
#mkdir p1tyger
#mv *P1Tyger* p1tyger/

#Set Plink PATH ###
export PATH=$PATH:/home/GWAS/all_enteric/clean_data/
export PATH=/home/GWAS/all_enteric/clean_data/:$PATH 

./plink --bfile P1Tyger_remap --geno 0.1 --mind 0.1 --make-bed --out P1Tyger_remap
./plink --bfile P1Tyger_remap --geno 0.05 --make-bed --out P1Tyger_remap
./plink --bfile T1T2_remap --geno 0.1 --mind 0.1 --make-bed --out T1T2_remap
./plink --bfile T1T2_remap --geno 0.05 --make-bed --out T1T2_remap
./plink --bfile VAST_remap --geno 0.1 --mind 0.1 --make-bed --out VAST_remap #2 people removed due to --geno
./plink --bfile VAST_remap --geno 0.05 --make-bed --out VAST_remap

#TOPMED Pre-impute checking ####
#make perl script later, apparently bash is not good for editing text files
#For now skip this section, go straight to Will Rayner Script

#Update Reference Genome ####

#convert ref genome to plink 
#update ref allele to fasta ref just to be on safe side
#then can --keep allele-ref-order

bcftools norm TOPMED38.vcf.gz -f GRCh38.fa --check-ref ws -Oz -o TOMED38_ref.vcf.gz

bcftools view -m2 -M2 -v snps all_enteric.vcf.gz -Oz -o all_enteric_bi.vcf.gz

./plink --vcf TOPMED38.vcf.gz --biallelic-only strict --keep-allele-order --out reference

#bcftools norm TOPMED38.vcf.gz -f GRCh38.fa --check-ref ws -Oz -o TOMED38_ref.vcf.gz

#plink --vcf TOPMED38.vcf.gz --a2-allele GRCh38.fa [1-based column index of ref alleles] [1-based column index of variant IDs] --recode vcf --real-ref-alleles

#read data library(data.table)
bim<-fread("reference.bim")

#Filter SNPs in autosomal chromosomes")
chrom <- c(seq(1,23),"X")
bim <- bim[bim$V1%in%chrom,]
#722005 snps left

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

plink --bfile reference --exclude reference.dups --make-bed --out reference

plink --bfile reference --bmerge VAST_remap --merge-mode 5 --out indiv_full

system("./plink --bfile /reference --bmerge refname --make-bed --out refname$filename_merged_rm")

#remove atcg from ref
#remove duplicates
#remove indels
#compare snp ID positions btween refd
#update position based on ref
#check chromosome position
#compare ref alt alleles
#update to ref alt alles of reference genome 
#compare frequency of a SNP, if the difference in frequency is greater that 20 percent, remove
#to do so make frequecy files of both snp and refence genome, then compare freq columns by snpID
#check for mismatches by merging, then using the flip snp output as below

system("./plink --bfile /reference/T1T_remap --bmerge refname --make-bed --out refname$filename_merged_rm")
#flip mismatched snps in datasets then remerge
system("./plink --bfile filename_clean --flip refname$filename_merged_rm.missnp --make-bed --out filename$flip_rm")

#for now just do will rayner script ##

#Will Rayner Pre-impute checking ####

#VAST
system("./plink --bfile vast_patch/VAST_remap --freq --out VAST_remap")

system("perl pre_impute.pl -b vast_patch/VAST_remap.bim -f VAST_remap.frq -r PASS.VariantsTOPMED38.tab.gz -h")

#T1T2
system("./plink --bfile t1t2/T1T2_remap --freq --out T1T2_remap")

system("perl pre_impute.pl -b t1t2/T1T2_remap.bim -f T1T2_remap.frq -r PASS.VariantsTOPMED38.tab.gz -h")


#P1Tyger 
system("./plink --bfile p1tyger/P1Tyger_remap --freq --out T1T2_remap")

system("perl pre_impute.pl -b p1tyger/P1Tyger_remap.bim -f P1Tyger_remap.frq -r PASS.VariantsTOPMED38.tab.gz -h")

#Update rs IDs ####

#update plink file with SNPs that we could extract, and update the snp coordinates using update map function # this was already done in VAST section
GSA-24v3-0_rsids.txt
system("./plink2 --bfile VAST_PATCH_rm --extract VAST_38_rm.snps --update-map GSA-24v3-0_rsids.txt --sort-vars --make-bed --out VAST_remap")



