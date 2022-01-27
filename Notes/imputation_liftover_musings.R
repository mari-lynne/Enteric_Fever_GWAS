#Directories

#set PLINK 1.9 to PATH
# #qcdir='~/qcdir' #rename these directories
#refdir='~/reference'
name='VAST_PATCH'
refname='Topmed38'

#Plan: #####

#Liftover enteric ####

#convert_enteric_files.bim to > BED for liftover using 
#
awk '{print "chr" $1, $4 -1, $4, $2 }' VAST_PATCH.bim > VAST_PATCH.tolift
#double check x and y chromosome codings, if as 23/24 need to recode for BEDsed 's/chr23/chrX/' | sed 's/chr24/chrY/' > \VAST_PATCH.tolift

cd map_ped 

#not working, try with map file instead 

#map format is 1  snp1   0  5000650 :chr SNP_ID GD Position)
1	1:103380393	124.1418	103380393

#BED format is
#1) chrom - The name of the chromosome (e.g. chr3, chrY)
#2)chromStart - The starting position of the feature in the chromosome.
#Chrom End 

awk '{print $1, $4-1, $4}' VAST_PATCH > VAST_PATCH.bed

#liftover study snps to 38 reference genome
online tool https://genome.ucsc.edu/cgi-bin/hgLiftOver 

#download lifted files

#update enteric files:
## ectract mapped variants
awk '{print $4}' $refdir/VAST_PATCH38 > $refdir/VAST_PATCH38.snps
# ectract updated positions
awk '{print $4, $3}' $refdir/VAST_PATCH38  > $refdir/VAST_PATCH38.pos

#remap SNP IDs and positions in PLINK 37 - 38 
plink --bfile $refdir/VAST_PATCH \
--extract $refdir/VAST_PATCH38.snps \
--update-map $refdir/VAST_PATCH38.pos \
--make-bed \
--out $refdir/VAST_PATCH38
mv $refdir/VAST_PATCH38.log $refdir/log

#copy and paste these files into GWAS directory for use
#cd to GWAS folder

#Rerun above steps for P1_Tyger and T1T2 data sets 

#Merge enteric #####

#Once data is updated merge enteric files

#remove duplicates in R #####
setwd('/home/mari/GWAS')
#install.packages("data.table")
library(data.table)
bim<-fread("VAST_PATCH38.bim")
bim.dups<-bim[duplicated(paste(bim$V1,bim$V4)),] 
duplicated_loci<-bim.dups$V2
write.table(file="duplicated_loci.txt",duplicated_loci,quote=FALSE,row.names = F,col.names = F)
system("./plink --bfile VAST_PATCH38 --exclude duplicated_loci.txt --make-bed --out VAST_PATCH38")

#T1T2_remove duplicates
bim<-fread("T1T2_38.bim")
bim.dups<-bim[duplicated(paste(bim$V1,bim$V4)),] 
duplicated_loci<-bim.dups$V2
write.table(file="duplicated_loci.txt",duplicated_loci,quote=FALSE,row.names = F,col.names = F)
system("./plink --bfile T1T2_38 --exclude duplicated_loci.txt --make-bed --out T1T2_38")


#Merge the VAST files T1T2 data first
# This will give us a misssnip file that we can later use to check for mismatches/flips
system("./plink --bfile VAST_PATCH38 --bmerge T1T2_38 --make-bed --out VAST_PATCH_T1T2")

#flip mismatched snps in og datasets then remerge
system("./plink --bfile VAST_cleaned --flip VAST_PATCH_T1T2-merge.missnp --make-bed --out VAST_flipped.hapmap-snps")

#remerge with the flipped snps
system("./plink --bfile VAST_flipped.hapmap-snps --bmerge VAST_PATCH_T1T2 --make-bed --out VAST_PATCH_T1T2merge")

#Merge this file with P1_Tyger file

