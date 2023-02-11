#!/bin/bash

main_dir=~/GWAS_22/gwas_final/merge/hla_tapas
dir=${main_dir}/mydata
refdir=${main_dir}/reference
name=typhoid2.IBD
refname=1000G.bglv4.bgl.phased


# run as source ~/GWAS_22/Enteric_GWAS/test.sh
cd ${main_dir}

# # Plink flip snps
# plink --bfile ${dir}/hla_chr6 --flip MyHLA-TAPAS/hla_chr6_impute.REF.bglv4.MERGED.FOUNDERS-merge.missnp --make-bed --out ${dir}/hla_6_flip

# plink --bfile ${dir}/hla_chr6 --flip MyHLA-TAPAS/hla_chr6_impute.REF.bglv4.MERGED.FOUNDERS-merge.missnp --make-bed --out ${dir}/hla_6_flip

# # Remove duplicates
# awk '{print $2, $4}' ${dir}/hla_6_flip.bim | cut -f -2 | uniq -d | cut -f -1 > multivar
# awk '{print $1}' multivar > multivar_snps
# # Do the same for duplicate name snps
# awk '{print $2}' ${dir}/hla_6_flip.bim | uniq -d >> multivar_snps

# plink2 --bfile ${dir}/hla_6_flip --rm-dup --make-bed --out ${dir}/hla_6_flip
# cat all_phase3_qc2.rmdup.mismatch multivar_snps | uniq > exclude_list
# plink2 --bfile ${dir}/hla_6_flip --exclude exclude_list --make-bed --out ${dir}/hla_6_flip_qc


# python HLA-TAPAS.py \
#     --target mydata/hla_6_flip_qc \
#     --reference reference/1000G.bglv4.bgl.phased \
#     --hped-Ggroup reference/1000G.EUR.Ggroup.hped \
#     --pheno mydata/pheno_typhoid.phe \
#     --hg 19 \
#     --out MyHLA-TAPAS/hla_6_impute \
#     --mem 4g \
#     --nthreads 4

 python3 -m SNP2HLA \
	--target mydata/typhoid2.IBD \
	--out mydata/hla_chr6_impute \
	--reference reference/1000G.bglv4 \
	--nthreads 2 \
	--mem 4g