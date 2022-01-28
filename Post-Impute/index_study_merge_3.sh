!#/bin/bash

#make index file for merging imputed data sets
#renamed from all_chromosomes attempts to study

bcftools index VAST_PATCH_imp.vcf.gz
bcftools index P1T1_imp.vcf.gz

echo DONE index

bcftools merge P1T1_imp.vcf.gz VAST_PATCH_imp.vcf.gz -Oz -o all_enteric.vcf.gz

echo DONE merging 
#took approx 1.5hs

#convert to vcf using plink2

./plink2 --vcf all_enteric.vcf.gz --double-id --allow-extra-chr 0 --make-bed --out all_enteric

echo DONE PLINK convert

#multiallelic varients, should have removed them during merge stage (I think theyre kept in as default)
bcftools view -m2 -M2 -v snps all_enteric.vcf.gz -Oz -o all_enteric_bi.vcf.gz

echo DONE varient pruning 

#cd Plink
./plink2 --vcf all_enteric_bi.vcf.gz --double-id --out enteric

echo DONE unziping vcf

#Notes ##
#Problems with vcf header ID could fix at an earlier stage? for now updated using sed and chr_rename.sh script

