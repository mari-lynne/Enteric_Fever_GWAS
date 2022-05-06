### File tracker
GWAS study data saved in:
~/GWAS_22/new_gwas/VAST_PATCH
~/GWAS_22/new_gwas/P1Tyger
~/GWAS_22/new_gwas/T1T2

Scripts are saved in:
~/GWAS_22/Enteric_Fever_GWAS/Pre-Impute..Post-Impute..ID_updating

### Work flow:

1)Prepare Plink data
2)Prepare and QC data for imputation
3)Topmed Imputation
4)Post-imputation processing
5)Association testing - GWAS
6)SNP annotation and mapping

1) **Prepare Plink data**

Original Plink data saved as .map, .ped, converted to plink binary format for compatibility with plink 1.9/2, and faster processing times

1a) Updated original VAST_PATCH plink files with concatenated IDs, see ID_updating/New_Ids.R, saved as VAST_PATCH2.bed
1b) Updated P1Tyger files to bed
1c) Updated T1T2 files to bed. No QC steps run at this point

2) **Prepare data for imputation**

Pre-impute script: 
- QC studies separately
- Liftover all files separately 37 > 38
- Run Will Rayner pre-impute QC per study > see Run.plink.sh
- Merge output of Will QC vcf files
- Zip VCF files
- Upload merged VCF files to TOPMED

Notes:
- I Would like to tidy this script/generalise it.
- Also add liftover commandline code
- Intergrate sex check/new ancestry scripts


3) **Topmed Imputation**
- P1Tyger, T1T2 merged > P1T1. Upload separately from VAST_PATCH study
- SNPs imputed using Topmed server and regerence genome https://imputation.biodatacatalyst.nhlbi.nih.gov/#!
- Eagle V.4 phasing and R2 cut off of 0.3

4) **Post-imputation processing**

4a) Download, merge and format:
- Data from Topmed is in separate chr vcf files, therefrore need to unzip the data and concatanate chr's for easier processing and merging, see post-impute_download.sh
- There was issues with Chr annotation in VCF files, fixed using sed in chr_rename.sh
- Rezipped then concat/merge? (P1T1, VAST_PATCH) using bcftools merge

4b)Run QC for SNPs, and merge studies
- SNPs filtered in Post_Impute_QC.RMD
- P1T1, VAST_PATCH merged with bmerge in PLINK - 9 mill snps across the studies
- Included rechallange participants as extra data (see ID_updates_march.R script
- Output is all_enteric plink files

4c) Check population stratification of the merged data set
- Variables recoded in post-impute_QC.R script
- PCA on all the enteric participants using reported ethicity data
- PCAs saved with covar file to inclde in association analysis (first 5 selected)
- Rs Ids updated in pos2rsID.R 

5) **Association testing**

GWAS_testing.RMD
- Run association tests of enteric fever challenge data across all enteric fever studies - typhoid and paratyphoid combined
- Run association tests for paratyphoid and typhoid separately
- Visualise GWAS results using a Manhattan Plot
- Associate the Fc receptor region SNPs with enteric fever
-Also contains code for recoding pheno/covar files

6) **SNP annotation and functional mapping**

Aims: To functionally characterise the SNPs that are significant in our typhoid association results

6a) Initially attempted to format data for FUMA, requires liftover back to GCHr 37, however the tool online siginifigance thresholds are too stringent for our analysis, see fuma.R 

6b) SNP_annotation.Rmd:
- Formatted data for SNP nexus (using just typhoid assoc as the associations between organisims were not common)
- Downloaded SNP annotation results from publically available databases
- Summarised data by SNP

7) **Fc Receptor analysis**

- Work completed in the Fc_receptor repositry
- Includes cis-EQTL analysis and haplotype analysis



