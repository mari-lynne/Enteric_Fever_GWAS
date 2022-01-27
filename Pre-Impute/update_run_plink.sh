!#/bin/bash

#Run-plink script needs to be modified to include ./ use sed
sed -i 's+plink+./plink+g' Run-plink.sh

#+ sign is the same as using back slashes but we changed it as we want to modify a backlash in this instance

sed 's/alleles/& --chr-output chrM/g' Run-plink.sh > Run-plink2.sh

echo updating plink commands Done

bash Run-plink2.sh

echo vcf conversion Done

#compress files to gz using bcf tools
for i in {1..22}; do bcftools sort VAST_remap-updated-chr$i.vcf -Oz -o VAST_remap-updated-chr-chr$i.vcf.gz; done

echo sort and vcf compression DONE

mkdir vcf_gz_toimpute
mv *.gz* vcf_gz_toimpute

mkdir log
mv *.log* log/

rm *nosex*

echo files organising DONE
