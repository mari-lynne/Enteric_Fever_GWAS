#!/bin/bash

dir=~/GWAS_22/gwas_final/merge/typhoid/vast/igg/chains

# Step 1: set up vars
# Step 2: run assoc code - check output
# Step 3: Filter by significance (use awk)

cd ${dir}

# select gene region

 plink2 --bfile ~/GWAS_22/gwas_final/merge/typhoid/QC/typhoid2.IBD \
--maf 0.1 \
--chr 14 \
--from-bp 105576437 \
--to-bp 106889844 \
--make-bed \
--out ighv_snps

# edit header
#for i in *.txt; do awk 'if (NR==1) {sed 's/.*\\.//'}'else print > ${i}; done
# sed default greedy
# update fixed output in R beforehand

#mkdir ${dir}/all

for i in *.txt; do plink2 --bfile ~/GWAS_22/gwas_final/merge/typhoid/QC/typhoid2.IBD \
  --pheno ~/GWAS_22/gwas_final/meta/ig-pheno.txt \
  --pheno-name ighv_usage \
  --covar ${i} \
  --covar-name age, bmi, sex, chall_vax, PC1, PC2, PC3, PC4, PC5 \
  --covar-variance-standardize \
  --glm \
  --adjust \
  --out ${dir}/all/${i}; done

# Filter output

#$7 = Test
# $12 = P val

for i in *.glm.linear; do awk '($12<0.05 && $7 == "ADD") {print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12}' ${i} > sig.${i};done

cd ${dir}/all
# save snps P < 0.000005
for i in *.adjusted; do awk '($10<0.000005)' ${i} > sig.${i} ;done
#for i in *.adjusted; do awk '($4<0.05)' ${i} > no_mt.sig.${i} ;done