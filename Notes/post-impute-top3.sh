#!/bin/bash'

##post-impute VAST-PATCH #####
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/428199/ffcceb5a822a8fd0cae918811b87a65d1eabb94c6ece0d993223e66124534880 | bash

echo download DONE

for f in *.zip; do unzip -P "A6lKZ5esXcD8lM" -d /home/mari/GWAS/all_enteric/VAST_PATCH/post-imp $f; done


echo unzip DONE

#uncompress .gz files to .vcf files

#for f in *.vcf.gz; do gunzip $f; done  

#echo decompress DONE

#musing ####

#merge vcf files 
#now in vcf folder
#can try cat in bash or vcf tools 

#for f in *.vcf; do bcftools concat $f -o VAST_PATCH_merge.vcf; done
#echo vcf merge DONE #only merges one file, basically doesnt work sequentitally/adds


#need to do i [1 in 22 sort of thing

#chr11.dose.vcf
#name=$(basename ${i} .vcf)

#for i in {1..22}*.vcf; do bcftools concat name=$(chr${i} .dose.vcf) -Oz -o VAST_PATCH_merge.vcf; done

#add zip output too


#convert to plink2 (keeps dosage data)

#./plink2 --vcf /home/mari/GWAS/TOPMED_impute/post-impute/VAST_PATCH_merge.vcf --out VAST_PATCH_full

#note merging of study data does not work with plink2 yet it seems
#so re-run pipe line for T1T2, P1Tyger and merge as vcf.gz files
#hopefully merge will work better after all strand align/lifted to the same ref genome 

#pre-process all study data merge when in plink 1.9 then do imputation on all
