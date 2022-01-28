Current Workflow:


- Have QCd data for VAST_PATCH, (P1,Tyger,T1,T2)- these were merged
- Rename P1TygerT1T2 to P1T1 for short
- Imputed snps for VAST_PATCH and P1T1 study cohorts separately
- Downloaded chr data, unzipped and concatanated chrs
- Rezipped then merged x2 .vcf.gz files using bcftools merge
- Merged data saved as all_enteric
- Convert to plink format

Improvements:

- Will try to fix merging pre-imputation and strand checking etc, can use these same scripts/bcftools
- Upload to topmed as one file
- Then download, unzip, concatanate
- Also check sample_IDs pre-imputation, I think the best solution would be to modify the FAM column as a study ID, that way when we make IDs using --double ID it should give unique IDs and save the massive hassle I went through
- I also started a pre-impute tool so need to go back to formatting that, could do similar for post (maybe post-Phd lol)
