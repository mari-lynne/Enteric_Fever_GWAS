#!/bin/bash
set -ux

# plink_data=vast
# plink_dir=~/GWAS_22/gwas_final/merge/typhoid/vast
# outdir=~/GWAS_22/gwas_final/merge/typhoid/vast/iga

plink_data=typhoid2.IBD
plink_dir=~/GWAS_22/gwas_final/merge/typhoid/QC
outdir=~/GWAS_22/gwas_final/eQTL/vast_tyg/Nov
data=D7_f2r

# same as outdir in Rscript :: EQTL = ~/GWAS_22/gwas_final/eQTL # Or adnob/antibody measure
cd ${outdir}
log=${outdir}/log

plink \
--bfile ${plink_dir}/${plink_data} \
--extract ${data}_snp_keep.txt \
--allow-no-sex \
--r2 \
--out ${data}