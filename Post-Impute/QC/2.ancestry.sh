#!/bin/bash
set -ux

# Notes:
# Script to impute ancestry of study data based on 1000g data in GChr38 form
# Run Script as: source ~/GWAS_22/Enteric_GWAS/Post-Impute/QC/2.ancestry.sh

# Steps:
# Download ref data
# Prune for linkage disequilibrium
# Clean both our ref and study data for AT/CG SNPs and duplicates/CHR mismatches
# Check SNP annotation, either update ref to CHR:POS:BP or study data to rsIDs
# Merge data sets (Plink1.9)

# Variables and Directories -----------------------------------------------------------

# name = study PLINK files
# ref = reference genome, in our case hg38
# qcdir = file to output and perform most of the script in
dir=~/GWAS_22/gwas_final/merge/QC

mkdir -p ${dir}/ancestry
qcdir=${dir}/ancestry

mkdir -p ${dir}/ref
refdir=${qcdir}/ref

mkdir -p ${refdir}/plink_log
log=${refdir}/plink_log

name=typhoid 
ref=all_hg38

# Download reference data ------------------------------------------------------------

cd ${refdir}

pgen=https://www.dropbox.com/s/e5n8yr4n7y91fyp/all_hg38.pgen.zst?dl=1
pvar=https://www.dropbox.com/s/cy46f1c8yutd1h4/all_hg38.pvar.zst?dl=1
sample=https://www.dropbox.com/s/3j9zg103fi8cjfs/hg38_corrected.psam?dl=1
wget $pgen
mv 'all_hg38.pgen.zst?dl=1' all_hg38.pgen.zst
plink2 --zst-decompress all_hg38.pgen.zst > all_hg38.pgen
wget $pvar
mv 'all_hg38.pvar.zst?dl=1' all_hg38.pvar.zst
wget $sample
mv 'hg38_corrected.psam?dl=1' all_hg38.psam

# Convert ref genome to plink  bed
plink2 \
--pfile ${refdir}/${ref} vzs --max-alleles 2 \
--allow-extra-chr \
--autosome \
--make-bed \
--out ${refdir}/${ref}

echo "1000G data download - DONE"

# 1) Tidy Study and Reference data ----------------------------------------------------
cd ${qcdir}

# 1a) Calculate Study and ref LD --------------------------------------------------------
plink2 \
--bfile ${dir}/${name}.IBD \
--indep-pairwise 50 5 0.2 \
--out ${qcdir}/${name}.LD

plink2 \
--bfile ${refdir}/${ref} \
--indep-pairwise 50 5 0.2 \
--out ${qcdir}/${ref}.LD

# 1b) Prune Study and Ref Data ---------------------------------------------------------
plink2 \
--bfile ${dir}/${name}.IBD \
--allow-extra-chr \
--extract ${qcdir}/${name}.LD.prune.in \
--make-bed \
--out ${qcdir}/${name}.LD

plink2 \
--bfile ${refdir}/${ref} \
--allow-extra-chr \
--extract ${qcdir}/${ref}.LD.prune.in \
--make-bed \
--out ${refdir}/${ref}

# 2) Remove AC-GT SNPs ----------------------------------------------------------------

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    ${qcdir}/${name}.LD.bim  > \
    ${qcdir}/${name}.acgt

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    ${refdir}/${ref}.bim  > \
    ${qcdir}/${ref}.acgt

# save ref_genome without ac/gt snps list to qc directory

echo AC-GT SNP list done

# 2a) Clean study data
plink2 \
--bfile ${qcdir}/${name}.LD \
--exclude ${qcdir}/${name}.acgt \
--make-bed \
--out ${qcdir}/${name}.no_acgt

# 2b) Clean ref data
plink2 \
--bfile ${refdir}/${ref} \
--allow-extra-chr \
--exclude ${qcdir}/${ref}.acgt \
--make-bed \
--out ${refdir}/${ref}.no_acgt


# 3) Remove duplicated snps and keep just atcg snps -----------------------------------

plink2 \
--bfile ${qcdir}/${name}.no_acgt \
--snps-only just-acgt \
--rm-dup exclude-all \
--make-bed \
--out ${qcdir}/${name}.cleaned

plink2 \
--bfile ${refdir}/${ref}.no_acgt \
--snps-only just-acgt \
--rm-dup exclude-all \
--make-bed \
--out ${qcdir}/${ref}.cleaned

rm -f *.no_acgt

# 4) Filter reference data for study SNPs ---------------------------------------------
# all files from here on in in qcdir
# 4a) Update ref chromsome annotation in bim file 
awk '{if ($1 != 0) print $2,"chr"$2}' ${ref}.cleaned.bim > updateFormat.txt

plink2 \
--bfile ${ref}.cleaned \
--update-name updateFormat.txt \
--make-bed \
--out ${ref}.cleaned

echo "${ref} chr annotation update- DONE"


# 4b) Filter reference data for study SNPs -----------------------------------------------

awk '{print $2}' ${name}.cleaned.bim > keep_list

plink2 \
--bfile ${ref}.cleaned \
--extract keep_list \
--make-bed \
--out ${ref}.forMerge

echo "Keep list - Done"


# 5a) Test merge _---------------------------------------------------------------------

# This gives us a list of snps which we can exclude
# You could also flip the snps and try remerging, but that's not necessary here

plink \
--bfile ${name}.cleaned \
--bmerge ${ref}.forMerge --merge-mode 6 \
--out ${ref}.merge_failures

grep "chr" ${ref}.merge_failures.log |\
awk 'BEGIN {OFS="\t"} {
if ($2 == "Multiple")
	print $7;
else
	print $3;
}'| sed -e "s/'//g" > exclude_list.txt 
#awk prints fields containing SNP IDs in log file
#Sed removes single quotes from output

echo "Exclude List - Done"

plink2 \
--bfile ${qcdir}/${ref}.forMerge \
--exclude exclude_list.txt \
--make-bed \
--out ${qcdir}/${ref}.cleanMerge

rm -f ${ref}.forMerge 

# 5b) Remerge -------------------------------------------------------------------------

plink \
--bfile ${name}.cleaned \
--bmerge ${ref}.cleanMerge \
--make-bed \
--out ${ref}.merged

plink \
--bfile ${ref}.merged \
--geno 0.01 \
--maf 0.01 \
--hwe 0.0001 \
--make-bed \
--out ${ref}.merged

# 6) Relatedness check of Ref data ----------------------------------------------------

# Fliter related participants in IKG data and use IBD
# https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3
# Download file - add wget code

nosib=https://www.dropbox.com/s/129gx0gl2v7ndg6/deg2_hg38.king.cutoff.out.id?dl=1

wget $nosib

plink2 --bfile ${ref}.merged \
--remove deg2_hg38.king.cutoff.out.id?dl=1 \
--make-bed \
--out ${ref}.merged.no_sib

plink \
--bfile ${ref}.merged.no_sib \
--genome \
--missing \
--out ${ref}.merged

# Remove one sample from each pair with pi-hat > 0.1875:
awk '$10 >= 0.1875 {print $1, $2}' ${qcdir}/${ref}.merged.genome |\
uniq > ${qcdir}/${ref}.merged.outliers.txt 

wc -l ${ref}.merged.outliers.txt

echo "outlier list - Done"

# Filter related samples in PLINK data 
plink2 \
--bfile ${qcdir}/${ref}.merged.no_sib \
--remove ${qcdir}/${ref}.merged.outliers.txt \
--make-bed \
--out ${qcdir}/${ref}.merged.IBD

# File cleaning -----------------------------------------------------------------------

# Move log files
mv *.log $log
rm *.nosex
rm *.zst

echo "File cleaning done"