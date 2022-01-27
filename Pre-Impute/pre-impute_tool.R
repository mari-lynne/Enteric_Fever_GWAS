#VAST_PATCH Pre-processing ####
setwd('/home/mari/GWAS/all_enteric')
library(data.table)
library(stringr)
library(dplyr)

#set filename and directories
wddir = "/home/mari/GWAS"
qcdir = "/home/mari/GWAS/QC"
filename = "VAST_PATCH"
path2plink = "/home/mari/GWAS/plink"
log = "/home/mari/GWAS/log"

#read data library(data.table)
bim<-fread(filename)

#Filter SNPs in autosomal or X chromosome")
chrom <- c(seq(1,23),"X")
bim <- bim[bim$V1%in%chrom,]

#keep only SNPS containing rs IDs #
bim$snp <- rownames(bim)
rs_bim <- bim %>% dplyr::filter(grepl("rs", V2))
#these are snps to keep ^ rs_bim
ValidSNPs <- rs_bim$V2

write.table(file="ValidSNPs_rm.txt", ValidSNPs, sep = "\t", quote = FALSE, col.names = F, row.names = F)

# Filter in valid SNPs
#input your file name here, precode it filename =

system("./plink --bfile filename --extract ValidSNPs_rm.txt --make-bed --out gwas_rm")

#filter CNV varients #exclude for now
CNV_bim <- bim %>% dplyr::filter(grepl("CNV", V2))

#remove duplicates #
bim.dups<-bim[duplicated(paste(bim$V1,bim$V4)),]
duplicated_snps <-bim.dups$V2 

bim <- read.table("gwas_rm.bim", stringsAsFactors=F)
rm_snps <- bim$V2[((bim$V5 == "A") & (bim$V6 == "T")) |
                 ((bim$V5 == "T") & (bim$V6 == "A")) |
                 ((bim$V5 == "C") & (bim$V6 == "G")) |
                 ((bim$V5 == "G") & (bim$V6 == "C"))]

#make exlusion list ##
exclude <- cbind(duplicated_snps, rm_snps)

write.table(file="exclude_snps_rm.txt", exclude, quote=FALSE,row.names = F,col.names = F)

#extract duplicates
system("./plink --bfile gwas_rm --exclude exclude_snps_rm.txt --snps-only just-acgt --make-bed --out filename$clean")

#Run Liftover ####

#convert to USCS BED file
#Add in chr and selecting relevant cols (-1 position)
system("awk '{print "chr" $1, $4 -1, $4, $2 }' filename$clean.bim > lift_snps")

liftOver lift_snps hg38ToHg19.over.chain.gz lifted.bed unlifted.bed

#get lifted snp IDs
system("awk '{print $4}' lifted.bed > lifted.snps")
# ectract updated positions
system("awk '{print $4, $3}' lifted.bed > lifted.pos")

#update hapmap plink file snp coordinates using update map function
system("./plink --bfile filename$clean --extract lifted.snps --update-map lifted.pos --make-bed --out filename$clean")

#tidy files ####
rm *rm*
mv *.log* log/

#Compare with Reference genome ####  
  #convert reference genome vcf to plink  
refname = "TOPMED_38"
#have function in scrpt to detect multiple input types, convert to plink 1.9

./plink --vcf refname --a2-allele --keep-allele-order --out refname
  
#merge with ref genome to get a list of snps that need to be flipped

#merging ####
system("./plink --bfile filename_clean --bmerge refname --make-bed --out refname$filename_merged_rm")

#flip mismatched snps in datasets then remerge
system("./plink --bfile filename_clean --flip refname$filename_merged_rm.missnp --make-bed --out filename$flip_rm")

#generate frequency file of ref genome 
./plink --freq --bfile filename_clean --out
./plink --freq --bfile refname --out

#column 5 has MAF left join filename_clean, refname $2 (or SNP)
filename_clean_rm <- left_join(filename_clean, refname, by = "SNP") %>%
mutate(MAF_dif = MAF.x - MAF.y) %>% if(MAF_dif >0.2,
{
  filter(filename_clean)
} else{
  print("PASS")
})
  
#run pre-impute QC  #####
plink --bfile filename_clean_rm --geno 0.1 --mind 0.1 --out
plink --bfile filename_clean_rm --geno 0.05 --mind 0.05  --out 

#Prepare for upload to TOPMED ####

#converts plink files to vcf, splits by chromosome then zips
#run script
bash pre_impute_top.sh

#need to add

#upload data #####
  
#Alt, run pre-impute script ######
  
filename = #use clean file from above

./plink --freq --bfile filename --out

perl HRC-1000G-check-bim.pl -b filename.bim  -f filename.frq -r PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz -h sh Run-plink.sh

#the run-plink script needs to be modified to include ./ use sed
sed -i 's+plink+./plink+g' Run-plink.sh

bash Run-plink.sh #writes.vcf files
#warning underscore in sample_ID's, might need to modify

#compress files to gz using bcf tools
for i in {1..22}; do bcftools sort filename-updated-chr$i.vcf -Oz -o filename$i.vcf.gz; done
echo sort DONE

#upload to TOPMED server
  

