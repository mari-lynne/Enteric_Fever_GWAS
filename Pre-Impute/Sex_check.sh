!#/bin/bash

#Aims:

#Check sex of genotyped individuals in enteric studies

#Data
#~/GWAS_22/new_gwas/VAST_PATCH/VAST_PATCH2.bim
#~/GWAS_22/new_gwas/T1T2/T1T2.bim
#~/GWAS_22/new_gwas/P1Tyger/P1Tyger.bim"

#Make update sex file from covar file
#See Sex_Check.R

plink \
--bfile T1T2 \
--update-sex update_sex.txt \
--make-bed \
--out T1T2_sex

#Error: Invalid centimorgan position on line 58721 of T1T2.bim #works with old plink

plink \
--bfile T1T2_sex \
--chr 23-24 \
--maf 0.01 \
--geno 0.05 \
--mind 0.05 \
--hwe 0.000001 \
--make-bed \
--out T1T2_sex_chr

#exclude list in r
plink \
--bfile T1T2_sex_chr \
--exclude exclude_snps_xy.txt \
--make-bed \
--out T1T2_sex_qc

#Warning: 69485 het. haploid genotypes present (see T1T2_sex_qc.hh ); many
#Total genotyping rate is 0.997583.
#5804 variants and 318 people pass filters and QC.

plink \
--bfile T1T2_sex_qc \
--split-x b37 no-fail\
--make-bed \
--out T1T2_sex_chr

#43 chromosome codes changed, is that 43 variants?
#hh genotypes stayed the same so not an autosome region problem

plink \
--bfile T1T2_sex_qc \
--set-hh-missing \
--indep-pairphase 20000 2000 0.5 \
--out T1T2_sex_LD

#Warning: 66700 het. haploid genotypes present (see T1T2_sex_LD.hh )
#5804 variants and 318 people pass filters and QC.
#Pruned 2065 variants from chromosome 23, leaving 3696.
#Pruned 14 variants from chromosome 25, leaving 29.
#Pruning complete.  2079 of 5804 variants removed.


plink \
--bfile T1T2_sex_qc \
--extract T1T2_sex_LD.prune.in \
--make-bed \
--out T1T2_sex_clean
#Warning: 42368 het. haploid genotypes present (see T1T2_sex_clean.hh ); many
#3725 variants and 318 people pass filters and QC.

plink \
--bfile T1T2_sex_clean \
--set-hh-missing \
--make-bed \
--out T1T2_sex_hh

#65 mismatches with het-haplotypes removed
#145 mismatches without

plink \
--bfile T1T2_sex_hh \
--check-sex ycount  \
--out T1T2_redo
