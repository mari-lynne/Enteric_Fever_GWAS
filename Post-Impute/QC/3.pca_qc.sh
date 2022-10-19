#!/bin/bash

# Script for calculating ancestry and study PCAs with plotting
# source ~/GWAS_22/Enteric_GWAS/Post-Impute/QC/3.pca_qc.sh

# Variables and Directories -----------------------------------------------------------
name=typhoid2
ref=all_hg38

dir=~/GWAS_22/gwas_final/merge/typhoid/QC
qcdir=${dir}/ancestry
refdir=${qcdir}/ref
log=${qcdir}/log

mkdir ${dir}/pca
pcadir=${dir}/pca
cd ${qcdir}

# Update fam/psam file and extract pop data ---------------------------------------------

# Earlier conversion of .psam to .fam lost pheno/pop data of 1KG participants
# Eigenvec file has FID (0s) IID PC1-10
# Get ID and pop from original psam file 
# psam file needs an updated FID column

# Add header to fam file using sed
# Paste fam file into psam to add the FAM column
sed -i '1i FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE' ${refdir}/${ref}.fam 
paste ${refdir}/${ref}.fam ${refdir}/${ref}.psam \
| cut -f 1,7-12 > ${qcdir}/${ref}.psam
# Remove extra columns 
awk '{print $1, $2, $6, $7, $8}' ${qcdir}/${ref}.psam > ${qcdir}/1kG.ID2Pop.txt

# Merged study/ancestry PCA ---------------------------------------------------------------

 plink2 \
 --bfile ${ref}-${name}.merged.no_sib \
 --pca \
 --out ${ref}-${name}.LD 

# Study PCA ------------------------------------------------------------------------
 
# PCA on LD pruned, IBD and sib checked data
 
 plink2 \
 --bfile ${qcdir}/${name}.LD \
 --pca \
 --out ${qcdir}/${name}.LD 

mv *.eigenval ${pcadir}
mv *.eigenvec ${pcadir}
mv 1kG.ID2Pop.txt ${pcadir}
mv *.log ${log}

# Plot in R -------------------------------------------------------------------------
#Rscript ${qcdir}/pop_pca_plot.R