#Date 28/01/22
setwd("~/GWAS_22/Enteric_Fever_GWAS/Pre-Impute/ID_updating")

#make new branch git 
git checkout -b ID_updates
#update remote
git push --set-upstream origin ID_updates

#Background ####
#Previously when merging T1T2P1 with VAST_PATCH there were some ID issues in that the genetics(blinded) IDs of some of the VP samples were unfortunately the same as the T1T2 study.
#I tried to work around it by modifying the files to add _26 to the end of duplicate IDs in VAST/PATCH...
#But because they were in associated files etc. this was more difficult the you would think/might be confusing for future analysis

#Proposed solution is to add the study code in the FID column of the FAM file (I do still need to check that this will stay separate in vcf etc or might have to re-separate)

#AIMS ####

#Update .fam FID files from VAST_PATCH to include study code
#PATCH == 3 VAST samples == 4 T1T2 == 1 Tyger == 2

setwd("~/GWAS_22/Enteric_Fever_GWAS/Pre-Impute/ID_updating")

#Steps ####
#So modifying fam files directly etc is not advised as was my hesitance previously and took a lot! of triple checking
#Plink actually has an updateids flag which I'll use :) #I will make the new file in excel - sorry R

https://www.cog-genomics.org/plink/1.9/data#update_indiv

#Get IDs from fam files (in all_enteric/clean_data/raw_data/VAST_PATCH)
#Export as .csv
#Make update files as below 

./Plink --bfile --update-IDs --out newfile

--update-ids expects input with the following four fields:
  
1) Old family ID
2) Old within-family ID
3) New family ID
4) New within-family ID

#SET PLINK TO PATH ####
#Save Plink in folder in PATH directory so can be called dircetly

echo $PATH
#add line to shell start up file
export PATH=$PATH:/place/with/the/file

#my version
#edit the .bashrc file in vim and add PATH file line in Vim

export PATH="$PATH:/home/mari/Plink"
source ~/.bashrc


