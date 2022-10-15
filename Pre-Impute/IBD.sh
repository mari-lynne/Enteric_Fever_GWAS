#!/bin/bash

# IBD script checking
#source ~/GWAS_22/Enteric_Fever_GWAS/Pre-Impute/IBD.sh

# Oct 22 - double check with all_enteric_qc3.cleaned (i.e pca file)

#dir=/home/mari/GWAS_22/new_gwas/Post-imp/merged
# home/mari/GWAS/all_enteric/clean_data/post-impute/Plink
dir=/home/mari/GWAS_22/new_gwas/just_typhoid
data=typhoid  #plink file prefix

cd /home/mari/GWAS_22/new_gwas/just_typhoid/QC

plink \
--bfile $dir/${data} \
--genome \
--out $data

#Remove one sample from each pair with pi-hat (% IBD) above threshold (0.1875 below):
awk '$10 >= 0.1875 {print $2, $4, $10}' $data.genome > pair_list.txt
awk '$10 >= 0.1875 {print $1, $2}' $data.genome | uniq > $data.outliers.txt 
wc -l $data.outliers.txt

echo outlier list - Done

#Use outlier list to filter samples in PLINK data 
plink2 \
--bfile $dir/${data} \
--remove $data.outliers.txt \
--make-bed \
--out $data.IBD


# Initial run 349 - 309 samples but that's because RC samples are in there
# Try without patch samples = 319 participants
# Double check VAST rpt samples too 
# 313 participants, 5 participants fail

#0 8585 8648 46 1 14 97
# IID1 IID2 PI_HAT
# 8585 8243 0.5248 #Cousins/sibling pair
# 8648 8414 0.9999 #sample done twice?
# 46 28 0.5065 #
# 14 31 0.4884

# So I removed legit siblings from file already, but, ran pca with repeat samples needs to be done separately