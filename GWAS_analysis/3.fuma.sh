#!/bin/bash

dir=~/GWAS_22/gwas_final/merge/typhoid/QC
fuma_dir=~/GWAS_22/gwas_final/merge/typhoid/assoc/fuma
meta_dir=~/GWAS_22/gwas_final/meta

name=typhoid2.IBD

cd $dir
mkdir log
log=$dir/log

# # 1. Print UCSC BED file ------------------------------------------------------------------
# awk '{print "chr" $1, $4 -1, $4, $2 }' $name.bim > ${name}_lift.bed

# echo BED file DONE

# # 2. Download (wget) chain files 38>37 
# wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -O hg38ToHg19.over.chain.gz
        
# echo Chain file download DONE

# # 3. liftBed: liftOver the bed file into build 37 
# liftOver ${name}_lift.bed hg38ToHg19.over.chain.gz ${name}_lift_out.bed ${name}_lift_out_unlifted.bed

# echo Liftover DONE

# # 4. Update data with liftover SNPs ------------------------------------------------------

# # extract mapped variants snp IDs
# awk '{print $4}' ${name}_lift_out.bed > ${name}_lift_map.snps
# # extract updated positions
# awk '{print $4, $3}' ${name}_lift_out.bed > ${name}_lift37coord.pos

# echo Make SNP var and coord files DONE

# #Update plink file with extracted snps and coordinates
# plink2 --bfile ${name} --extract ${name}_lift_map.snps --update-map ${name}_lift37coord.pos --sort-vars --make-pgen --out ${name}_remap37
# #Error: Fixed-width .bed/.pgen output doesn't support sorting yet.  Generate a regular sorted .pgen first, and then reformat it
# plink2 --pfile ${name}_remap37 --make-bed --out ${name}_remap37

# echo Update Plink files DONE

# # 5 Gchr37 Assoc Analysis ------------------------------------------------------

# plink2 --bfile ${name}_remap37 \
# --maf 0.1 \
# --pheno ${meta_dir}/pheno_typhoid.txt \
# --pheno-name Diagnosed \
# --covar ${meta_dir}/typhoid_pcacovar.txt \
# --covar-name sex, age, chall_vax, PC1, PC2, PC3, PC4, PC5 \
# --covar-variance-standardize \
# --glm firth \
# --out ${fuma_dir}/typhoid

# echo Association testing DONE

# Move log files

mv *.log $log

# 6) filter for ADD ------------------------------------------------------------

# First print header then filter
cd $fuma_dir

awk 'NR==1{print $1, $2, $4, $6, $12, $9, $10} NR!=1{if ($7 == "ADD") {print $1, $2, $4, $6, $12, $9, $10}}' \
typhoid.Diagnosed.glm.firth > typhoid_37_gwas.txt

# TODO update this as awk script - but for now rename header and tab seperate in R, then gzip file

