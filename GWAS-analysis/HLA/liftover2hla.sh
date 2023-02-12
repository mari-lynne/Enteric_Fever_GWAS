#!/bin/bash

main_dir=~/GWAS_22/gwas_final/merge/HLA-TAPAS-master
dir=${main_dir}/mydata
refdir=${dir}/reference
name=typhoid2.IBD_remap37
refname=1000G.bglv4.bgl.phased
mkdir log
log=$dir/log

# run as source ~/GWAS_22/Enteric_GWAS/liftover2hla.sh

cd ${dir}

# # 1. Print UCSC BED file ------------------------------------------------------------------
# awk '{print "chr" $1, $4 -1, $4, $2 }' ${name}.bim > ${name}_lift.bed

# echo "UCSC BED file DONE"

# # 2. Download (wget) chain files 38>37 
# wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz -O ${refdir}/hg19ToHg18.over.chain.gz

# echo "Chain file download DONE"

# # 3. liftBed: liftOver the bed file into build 36

# liftOver ${name}_lift.bed ${refdir}/hg19ToHg18.over.chain.gz ${name}_lift_out.bed ${name}_lift_out_unlifted.bed

# echo "Liftover DONE"

# # 4. Update data with liftover SNPs ------------------------------------------------------

# # extract mapped variants snp IDs
# awk '{print $4}' ${name}_lift_out.bed > ${name}_lift_map.snps
# # extract updated positions
# awk '{print $4, $3}' ${name}_lift_out.bed > ${name}_lift36coord.pos

# echo "Make snp variants and coord files DONE"

# #Update plink file with extracted snps and coordinates
# plink2 --bfile ${name} --extract ${name}_lift_map.snps --update-map ${name}_lift36coord.pos --sort-vars --make-pgen --out ${name}_36

# plink2 --pfile ${name}_36 --make-bed --out ${name}_36

echo "Update Plink files DONE"

# get HLA genes
plink2 --bfile ${name}_36 --from-kb 29299 --to-kb 33884 --chr 6 --make-bed --out hla_hg18 

mv *.log ${log}

