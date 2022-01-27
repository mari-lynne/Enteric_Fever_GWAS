#!/bin/bash

#Rrun-plink script needs to be modified to include ./ use sed
#sed -i 's+plink+./plink+g' Run-plink.sh

#--chr-output chrM --recode vcf

#sed '/\alleles\a --chr-output chrM' Run-plink.sh

#sed -i --real-ref-alleles --recode vcf

#bash Run-plink.sh #writes.vcf files
#warning underscore in sample_ID's, might need to modify
#add --const-fid

#convert vcf to plink files using PLINK2 (binary format)

#for f in {1..22}; do ./plink2 --vcf VAST_38_remap-updated-chr$f.vcf --out VAST_remap38-updated-chr$f; done
#echo VCF to PLINK DONE

#compress files to gz using bcf tools
#for i in {1..22}; do bcftools sort VAST_38_remap-updated2-chr$i.vcf -Oz -o VAST_38_remap-updated2-chr$i.vcf.gz; done
#echo sort DONE

#bcftools sort VAST_38_remap-updated2-chr18.vcf -Oz -o VAST_38_remap-updated2-chr18.vcf.gz


