#!bin/bash

plink --bfile hla_chr6 --recode --out hla_chr6



$ python HLA-TAPAS.py \
    --target mydata/hla_chr6 \
    --reference reference/1000G.bglv4 \
    --hped-Ggroup reference/1000G.EUR.Ggroup.hped \
    --pheno mydata/pheno_typhoid.phe \
    --hg 19 \
    --out MyHLA-TAPAS/hla_chr6 \
    --mem 4g \
    --nthreads 4
