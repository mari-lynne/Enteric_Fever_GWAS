#Aims merge all our enteric data sets
#save all files in the same folder

#Need our reference genome converted to BED for Liftover:
#download tomed files

#convert_enteric_files.bim to 

#convert either vcf > Plink > BED using 
#
awk '{print "chr" $1, $4 -1, $4, $2 }' VAST_PATCH.bim | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > VAST_PATCH.tolift

#used online and got a .bed file with lifted co-ordinates NCBI build 36 to CGRCh37
  # output reads Successfully converted 1439670 records: Conversion failed on 946 records. 
  #download file and rename as HapMapIII_CGRCh37 #save in refdir folder
  
  #get mapped positions
  #awk is a command that you can use to deal with text or strings and prints them/reformats $ means the fourth field or column
  
  #try without removeing a t c-g snps first
  
#setup directories
  
#or use plink QC script for each study against the ref genome
  perl HRC-1000G-check-bim.pl -b VAST_PATCH_QC3_rmdup.bim  -f VAST_PATCH_freq.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h sh Run-plink.sh
  
  #qcdir='~/qcdir' #rename these directories
  #refdir='~/reference'
  name='VAST_PATCH'
  refname='Topmed38'
  
  mkdir -r $qcdir/plink_log
  
  #Plan: #####
  
  #convert_enteric_files.bim to > BED for liftover using 
  #
  awk '{print "chr" $1, $4 -1, $4, $2 }' VAST_PATCH.bim | \
  sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > \
  VAST_PATCH.tolift
  
  #liftover study snps to 38 reference genome
  online tool https://genome.ucsc.edu/cgi-bin/hgLiftOver 
  
  #update enteric files:
  ## ectract mapped variants
  awk '{print $4}' $refdir/VAST_PATCH38 > $refdir/VAST_PATCH38.snps
  # ectract updated positions
  awk '{print $4, $3}' $refdir/VAST_PATCH38  > $refdir/VAST_PATCH38.pos
  
  #remap SNP IDs and positions in PLINK 37 - 38 
  plink --bfile $refdir/VAST_PATCH \
  --extract $refdir/VAST_PATCH38.snps \
  --update-map $refdir/VAST_PATCH38.pos \
  --make-bed \
  --out $refdir/VAST_PATCH38
  mv $refdir/VAST_PATCH38.log $refdir/log
  
  #copy and paste these files into GWAS directory for use
  #cd to GWAS folder
  
  #Rerun above steps for P1_Tyger and T1T2 data sets 
  
  
  
  