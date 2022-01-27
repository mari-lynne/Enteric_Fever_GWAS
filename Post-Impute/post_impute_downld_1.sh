#!/bin/bash'

#Post-impute 
#Download and unzip files from topmed server
#files will be as separate chromosomes therefore concat with bcftools 

#Unzip

for f in *.zip; do unzip -P "vPFwBegJFr26kX" -d  ~/GWAS/all_enteric/P1T1merge/post-imp $f; done

echo Unzip DONE

#bcftools concat
#make files variable first #then concat

#redid concat using concat.sh script instead as this aligned it to the ref genome I think

#files=$(for i in $(seq 1 22);do echo ~/GWAS/all_enteric/P1T1merge/post-imp$i.dose.vcf.gz;done)
#bcftools concat $files -Oz -o allchromosomes.vcf.gz

#echo Concat DONE



