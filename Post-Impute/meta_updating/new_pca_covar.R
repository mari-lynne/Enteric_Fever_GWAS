#New pheno files #from ancestrt/QC
#wd "/home/mari/GWAS_22/new_gwas/QC/ancestry"

#Need to be split by typhoid paratyphoid
#Include newly calculated PCs (enteric.pop.pca.R
eigenval <- fread("1KG.merged.eigenval")
pve <- data.frame(PC = 1:10, pve = (eigenval/sum(eigenval)*100))
pve$PC <- as.factor(pve$PC)


#PCA no IKG ####
pca2 <- fread("all_enteric_QC3.cleaned.eigenvec")
eigenval2 <- fread("all_enteric_QC3.cleaned.eigenval")

pheno <- fread("~/GWAS_22/new_gwas/meta/final_metadata/all_pheno_num.txt")

merge2 <- full_join(pca, pheno, by = "IID")
merge2$Diagnosed<- as.factor(merge2$Diagnosed)
ggplot(data=merge2,aes(PC1, PC2, colour = Diagnosed)) + geom_point()
ggplot(data=merge2,aes(PC2, PC3, colour = Diagnosed)) + geom_point()

merge3 <- full_join(PCA, pheno, by = "IID")) #337

#Get covars
#T_covar <- fread("~/GWAS_22/new_gwas/meta/final_metadata/T_covar.txt")#328
covar <- fread("~/GWAS_22/new_gwas/meta/final_metadata/all_covar_june.txt")

library(tidylog)
T_PCA <- left_join(T_covar, merge2, by = "IID")#matched = 319, 9 rows only in x, 18 rows only in y

PCA_covar <- left_join(covar, merge2, by = "IID") #matched rows 331, rows only in x = 15, rows only in y = 6 #326 matched rows
PCA_covar <- left_join(merge2,covar, by = "IID")
missing <- anti_join(merge2,covar, by = "IID") #These were failures mostly nothing to find

table(PCA_covar$Vaccine) #Ty21a recoded as 2 instead of 1 by mistake #corrected

PCA_covar = PCA_covar[,-c(8:13)] #remove extra PCs

#filter for just typhoid P = 3, T = 1, TN = 2
T_covar2 <- filter(T_PCA, Challenge != "3") #257 participants
#filter for just paratyphoid P = 3, T = 1, TN = 2
P_covar2 <- filter(T_PCA, Challenge == "3") #71 participants

#Save
write.table(T_covar2, file = "~/GWAS_22/new_gwas/meta/final_metadata/T_covar_Jun.txt", row.name = F, col.names = T, sep = "\t", quote = F)

write.table(P_covar2, file = "~/GWAS_22/new_gwas/meta/final_metadata/P_covar_Jun.txt", row.name = F, col.names = T, sep = "\t", quote = F)

save.image(file = "~/GWAS_22/new_gwas/meta/final_metadata/pheno_covars_Jun.RData")
