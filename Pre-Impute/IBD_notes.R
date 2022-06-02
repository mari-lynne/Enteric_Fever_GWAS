#IBD notes ####

#IBD ####
#Remove with pi-hat scores

pairs <- fread("pair_list.txt")

#We want to remove the sibling pair as the PC will be swayed by the presence of siblings rather than diverging by actual ancestry/population structure differences.
#Therefore also need to remove the RPT samples

#Need to remove these samples before doing ancestry plots and calculating PCAs for model
#Also do need to separate the RC samples from the data set if calculating PCA

#So ideally do above, filter sibs from main all_enteric_QC2 data
#Then also split based on covar typhoid, and paratyphoid sets
#The paratyphoid set will have the patch rc participants but bc they got typhoid first, there should not be a duplicate

#names in meta/final_metadata/T_covar

#Double check if you can see any immediate labelling issues
#8648, 8414
#5019	Tyger_D00_8414	8414	TYGER-5019	Tyger_D00_8414
#5020	Tyger_D00_8822	8648	TYGER-5020	Tyger_D00_8822

#The samples are one after another in covar file, they are not in fam file
#Different ages tho so can't be twins

#Check sex
#8822 comes up as a problem sex so remove that one of the pair :)
#8585 also problem sex 

#So to remove:
#Rpt samples 1-rpt, 97 rpt-26

#Genuine relatives, probs siblings
#One from each pair, either
#46 28 0.5065
#14 31 0.4884

#Keep based on most interesting pheno
#remove vi-ps participants maybe until I do a separate gwas for vaccine efficacy

#Using cat and cut made failqc_ID_list.txt

#filter IDs

#Need to get FIDs for samples, they are all 0
#337 people remaining > all_enteric_QC3