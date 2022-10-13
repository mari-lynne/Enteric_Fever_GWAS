#!/bin/bash
#set -uex

# Notes:
# Script to impute ancestry of study data based on 1000g data in GChr38 form
# Run Script as source ~/GWAS_22/Enteric_Fever/GWAS/Post-Impute/Post-imputation_QC/ancestry_check_38.sh

# Steps:
# Download ref data
# Prune for linkage disequilibrium
# Clean both our ref and study data for AT/CG SNPs and duplicates/CHR mismatches
# Check SNP annotation, either update ref to CHR:POS:BP or study data to rsIDs
# Merge data sets (Plink1.9) 
# PCA of data colour code by ancestry (in R)

# Set up directories -------------------------------------------------------

name=all_enteric_QC3 #name of study PLINK files
refname=all_hg38

studydir=/home/mari/GWAS_22/new_gwas/QC

#Set up folder for reference genome download
#mkdir -p $studydir/ref
refdir=$studydir/ref

#Set up log folder
#mkdir -p $refdir/plink_log
log=$refdir/plink_log

#Set up folder for ancestry qc of study and reference data
#mkdir -p $studydir/ancestrty
qcdir=$studydir/ancestry
# #qcdir will contain the cleaned study and refernce data

# #Download refernce data -------------------------------------------------------

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

echo 1000g data download DONE

# Move files into new folders
mv *$refname* $refdir 
cp *$name* $qcdir

# Convert to Bed
plink2 \
--pfile $refdir/$refname vzs --max-alleles 2 \
--allow-extra-chr \
--autosome \
--make-bed \
--out $refdir/$refname
mv $refdir/$refname.log $log

# 1) Tidy Study and Reference data ---------------------------------------------------

cd $qcdir

# 1a) Prune Study Data --------------------------------------------------------------
plink2 \
--bfile $qcdir/$name \
--indep-pairwise 50 5 0.2 \
--out $qcdir/$name.LD

plink2 \
--bfile $qcdir/$name \
--allow-extra-chr \
--extract $qcdir/$name.LD.prune.in \
--make-bed \
--out $qcdir/$name.LD

# 2) Remove AC-GT SNPs -----------------------------------------------------
# (these snps are difficult to merge)

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $qcdir/$name.LD.bim  > \
    $qcdir/$name.acgt

cd $refdir
awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $refdir/$refname.bim  > \
    $qcdir/$refname.acgt #save ref_genome without ac/gt snps list to qc directory
cd $qcdir

echo AC-GT SNP list done

# plink2 \
--bfile $qcdir/$name.LD \
--exclude $qcdir/$name.acgt \
--make-bed \
--out $qcdir/$name.no_acgt  #clean study data
   
# plink2 \
--bfile $refdir/$refname \
--allow-extra-chr \
--exclude $qcdir/$refname.acgt \
--make-bed \
--out $qcdir/$refname.no_acgt 


# 2b) Remove duplicated snps and keep just atcg snps -------------------------------------------------
plink2 \
--bfile $qcdir/$name.no_acgt \
--snps-only just-acgt \
--rm-dup exclude-all \
--make-bed \
--out $qcdir/$name.cleaned

# plink2 \
--bfile $qcdir/$refname.no_acgt \
--snps-only just-acgt \
--rm-dup exclude-all \
--make-bed \
--out $qcdir/$refname.cleaned

cd $qcdir

# rm -f *$refname.no_acgt* 

# 4) Filter reference data for study SNPs -------------------------------------------------
# 4a) Update ref chromsome annotation in bim file 
awk '{if ($1 != 0) print $2,"chr"$2}' $refname.cleaned.bim > updateFormat.txt
echo make update name file DONE

plink2 \
--bfile $refname.cleaned \
--update-name updateFormat.txt \
--make-bed \
--out $refname.cleaned


# 4b) Filter reference data for study SNPs -----------------------------------------------
awk '{print $2}' $name.cleaned.bim > keep_list

echo keep list - Done

plink2 \
--bfile $refname.cleaned \
-extract keep_list \
--make-bed \
--out $refname.forMerge

echo Keep list - Done


# Test merge _----------------------------------------------------------------
# This gives us a list of snps which we can exclude
# You could also flip the snps and try remerging, but that's not necessary here

plink \
--bfile $name.cleaned \
--bmerge $refname.forMerge --merge-mode 6 \
--out $refname.merge_failures


grep "chr" $refname.merge_failures.log |\
awk 'BEGIN {OFS="\t"} {
if ($2 == "Multiple")
	print $7;
else
	print $3;
}'| sed -e "s/'//g" > exclude_list.txt 
#awk prints fields containing SNP IDs in log file
#Sed removes single quotes from output

echo Exclude List - Done

plink2 \
--bfile $qcdir/$refname.forMerge \
--exclude exclude_list.txt \
--make-bed \
--out $qcdir/$refname.cleanMerge

rm -f $refname.forMerge 
# Remerge -----------------------------------------------------------------------

plink \
--bfile $name.cleaned \
--bmerge $refname.cleanMerge \
--make-bed \
--out 1KG.merged

plink \
--bfile 1KG.merged \
--geno 0.01 \
--maf 0.01 \
--hwe 0.0001 \
--make-bed \
--out 1KG.merged

# Move log files
mv *.log $log
rm *.nosex

echo File cleaning done

# Relatedness check ----------------------------------------------------------------
# Fliter related participants in IKG data and study data
# calculate IBD using plink --genome flag and filter using pi-hat score
# Should rerun this pre-pcs further analysis

# Remove related samples 1KG
# https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3
# Download file - add wget code

nosib=https://www.dropbox.com/s/129gx0gl2v7ndg6/deg2_hg38.king.cutoff.out.id?dl=1

wget $nosib

plink2 --bfile 1KG.merged \
--remove deg2_hg38.king.cutoff.out.id?dl=1 \
--make-bed \
--out 1KG.merged.no_sib

plink \
--bfile 1KG.merged.no_sib \
--genome \
--missing \
--out 1KG.merged

# Remove one sample from each pair with pi-hat (% IBD) above threshold (0.1875 below):
awk '$10 >= 0.1875 {print $1, $2}' $qcdir/1KG.merged.genome |\
uniq > $qcdir/1KG.merged.outliers.txt 

wc -l 1KG.merged.outliers.txt

echo outlier list - Done

#Use outlier list to filter samples in PLINK data 
plink2 \
--bfile $qcdir/1KG.merged.no_sib \
--remove $qcdir/1KG.merged.outliers.txt \
--make-bed \
--out $qcdir.1KG.merged.IBD
# Works :)

# PCA ----------------------------------------------------------------------------------------------------------------

 plink2 \
 --bfile 1KG.merged.IBD\
 --pca \
 --out 1KG.merged

# Update fam file and extract pop data -------------------------------------------------------------------------------

# Earlier conversion of .psam to .fam lost pheno/pop data of 1KG participants
# Eigenvec file has FID (0s) IID PC1-10
# get ID and pop from original psam file 
# psam file needs an updated FID column

sed -i '1i FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE' $refdir/$refname.fam #Add header to fam file using sed
paste $refdir/$refname.fam $refdir/$refname.psam |\ #paste fam file into psam to add the FAM column
cut -f 1,7-12 > $qcdir/$refname.psam #remove extra columns 

awk '{print $1, $2, $6, $7, $8}' $qcdir/$refname.psam > $qcdir/1kG.ID2Pop.txt #file for R to use


# Plot in R -----------------------------------------------------------------------------------------------------------

Rscript $qcdir/pop_pca_plot.R