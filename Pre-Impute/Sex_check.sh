!#/bin/bash
set -uex
#set so doesn't close upon erroring out 

#Aims:

#Check sex of genotyped individuals in enteric studies

#Data
#~/GWAS_22/new_gwas/VAST_PATCH/VAST_PATCH2.bim #renamed to VAST_PATCH in sexdir
#~/GWAS_22/new_gwas/P1Tyger/P1Tyger.bim"

#See Sex_Check.R for making update sex files
#CD to sex dir in terminal then run script as source ~/GWAS_22/Enteric_Fever_GWAS/Pre-Impute/Sex_check.sh


data=VAST_PATCH
sexdir=~/GWAS_22/new_gwas/VAST_PATCH/sex_check
echo $sexdir

cd $sexdir

plink \
--bfile $sexdir/${data} \
--update-sex update_sex.txt \
--make-bed \
--out $sexdir/${data}_sex

#Error: Invalid centimorgan position on line 58721 of $sexdir/${data}.bim
#Need to use plink 1.9

plink \
--bfile $sexdir/${data}_sex \
--chr 23-24 \
--maf 0.01 \
--geno 0.05 \
--mind 0.05 \
--hwe 0.000001 \
--make-bed \
--out $sexdir/${data}_sex_chr

#exclude list in r
plink \
--bfile $sexdir/${data}_sex_chr \
--exclude exclude_snps_xy.txt \
--make-bed \
--out $sexdir/${data}_sex_qc

#Prune SNPs
plink \
--bfile $sexdir/${data}_sex_qc \
--set-hh-missing \
--indep-pairphase 20000 2000 0.5 \
--out $sexdir/${data}_sex_LD

#Exclude pruned SNPs
plink \
--bfile $sexdir/${data}_sex_qc \
--extract $sexdir/${data}_sex_LD.prune.in \
--make-bed \
--out $sexdir/${data}_sex_clean

plink \
--bfile $sexdir/${data}_sex_clean \
--set-hh-missing \
--make-bed \
--out $sexdir/${data}_sex_hh

#Add in ycount statistic
plink \
--bfile $sexdir/${data}_sex_hh \
--check-sex ycount  \
--out $sexdir/${data}_redo

mkdir log

mv *.log log