#!/bin/bash'

#Post-impute 

#Unzip

#for f in *.zip; do unzip -P "vPFwBegJFr26kX" -d  ~/GWAS/all_enteric/P1T1merge/post-imp $f; done


echo Unzip DONE

#bcftools concat
#make files variable first #then concat

files=$(for i in $(seq 1 22);do echo ~/GWAS/all_enteric/P1T1merge/post-imp$i.dose.vcf.gz;done)
bcftools concat $files -Oz -o allchromosomes.vcf.gz

echo Concat DONE



