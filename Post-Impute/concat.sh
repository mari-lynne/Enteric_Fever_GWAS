#!bin/bash'

bcftools concat chr1.dose.vcf.gz chr2.dose.vcf.gz chr3.dose.vcf.gz chr4.dose.vcf.gz chr5.dose.vcf.gz chr6.dose.vcf.gz chr7.dose.vcf.gz chr8.dose.vcf.gz chr9.dose.vcf.gz chr10.dose.vcf.gz chr11.dose.vcf.gz chr12.dose.vcf.gz chr13.dose.vcf.gz chr14.dose.vcf.gz chr15.dose.vcf.gz chr16.dose.vcf.gz chr17.dose.vcf.gz chr18.dose.vcf.gz chr19.dose.vcf.gz chr20.dose.vcf.gz chr21.dose.vcf.gz chr22.dose.vcf.gz -Ou | 
bcftools view -Ou -i 'R2>0.3' |
bcftools norm -Ou -m -any |
bcftools norm -Ou -f GRCh38.fa |
bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' -o new_allchromosomes.converted.R2_0.3.vcf.gz

echo Concat Done

#Done for both P1T1 and VAST/PATCH before merging
#Try redoing with for loop like before
#files=$(for i in $(seq 1 22);do echo ~/GWAS/all_enteric/P1T1merge/post-imp$i.dose.vcf.gz;done)
#bcftools concat $files -Oz -o allchromosomes.vcf.gz
