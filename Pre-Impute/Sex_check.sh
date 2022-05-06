!#/bin/bash

#Aims:

#Check sex of genotyped individuals in enteric studies

#Data
#~/GWAS_22/new_gwas/VAST_PATCH/VAST_PATCH2.bim
#~/GWAS_22/new_gwas/T1T2/T1T2.bim
#~/GWAS_22/new_gwas/P1Tyger/P1Tyger.bim"

#Make update sex file from covar file
#See Sex_Check.R

plink2 \
--bfile T1T2 \
--update-sex update_sex.txt \
--make-bed \
--out T1T2_sex