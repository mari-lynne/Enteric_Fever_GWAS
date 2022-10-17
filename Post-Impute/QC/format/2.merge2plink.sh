#!/bin/bash
set -ux

# 3.merge.sh

# Merge separate imputed study data
# Directory and Variables -----------------------------------------------
STUDY1=VAST
STUDY2=P1T1

DIR=~/GWAS_22/gwas_final
DIR1=${DIR}/${STUDY1}
DIR2=${DIR}/${STUDY2}

#mkdir merge
OUT_DIR=${DIR}/merge

# Index and merge
bcftools index ${DIR1}/${STUDY1}_concatR2.bcf \
bcftools index ${DIR2}/${STUDY2}_concatR2.bcf \
bcftools merge -m none ${DIR2}/${STUDY2}_concatR2.bcf ${DIR1}/${STUDY1}_concatR2.bcf -Ob -o ${OUT_DIR}/merged_data.bcf 
echo "Finished index and merge of ${STUDY1} and ${STUDY2}"

# -m none   ..  no new multiallelics, output multiple records instead
# Might have same ids, potensh filter for typhoid IDs first

# Plink conversion
plink2 --bcf ${OUT_DIR}/merged_data.bcf \
--max-alleles 2 \
--double-id \
--hwe 0.00001 \
--maf 0.05 \
--make-bed \
--out ${OUT_DIR}/merged_data

echo "Plink conversion DONE"

mkdir ${DIR}/log
mv *.log log