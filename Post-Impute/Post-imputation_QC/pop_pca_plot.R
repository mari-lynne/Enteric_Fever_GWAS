###Date: 2022/06/02
###Purpose: Plot user data projected on 1KG PCs, greying out 1KG individuals

setwd("/home/mari/GWAS_22/new_gwas/QC/ancestry")

##Load packages
library(ggplot2)
library(colorspace)
library(readr)
library(data.table)
library(dplyr)

##Initialise command line arguments
args <- commandArgs(TRUE)

##Load data
root <- args[1]

PCA <- read_delim("1KG.merged.eigenvec", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
##Rename columns
colnames(PCA) <- c("FID","IID", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

#merge with Population data
pop <- read_delim("1kG.ID2Pop.txt", 
                         delim = " ", escape_double = FALSE, 
                         trim_ws = TRUE)

merge <- full_join(PCA, pop, by = "IID")


#colour NAs as grey = our study data
#overlay that with the pop data
##Define colour palette
KG_Palette<-heat_hcl(length(unique(merge$Population)), h = c(300, 75), c. = c(35, 95), l = c(15, 90), power = c(0.8, 1.2), fixup = TRUE, gamma = NULL, alpha = 0.8)

#Large pop
KG_Palette2<-heat_hcl(length(unique(merge$SuperPop)), h = c(300, 75), c. = c(35, 95), l = c(15, 90), power = c(0.8, 1.2), fixup = TRUE, gamma = NULL, alpha = 0.8)

#Arrange so NAs are at the bottom of DF
#Rename NAs

merge$SuperPop <- merge$SuperPop %>% replace_na("Study")
merge$SuperPop<- as.factor(merge$SuperPop)
levels(merge$SuperPop)

merge$Population <- merge$Population %>% replace_na("ZStudy")
merge$Population<- as.factor(merge$Population)
levels(merge$Population)

merge <-merge[order(merge$SuperPop),]

pdf("KG.pop_strat_PCA.pdf")
ggplot(merge,aes(PC1,PC2, colour=SuperPop)) +geom_point() + scale_colour_manual(values = KG_Palette2) +theme_bw()
dev.off()


pdf("KG.subpop_strat_PCA.pdf")
ggplot(merge,aes(PC1,PC2, colour=Population)) +geom_point() + scale_colour_manual(values = KG_Palette) +theme_bw()
dev.off()

a=ggplot(merge,aes(PC1,PC2, colour=SuperPop)) +geom_point() + scale_colour_manual(values = KG_Palette2)
b=ggplot(merge,aes(PC2,PC3, colour=SuperPop)) +geom_point() + scale_colour_manual(values = KG_Palette2)
c=ggplot(merge,aes(PC3,PC4, colour=SuperPop)) +geom_point() + scale_colour_manual(values = KG_Palette2)
d=ggplot(merge,aes(PC4,PC5, colour=SuperPop)) +geom_point() + scale_colour_manual(values = KG_Palette2)


relabel <- c("ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","STU","TSI","YRI","Study")

e=ggplot(merge,aes(PC1,PC2, colour=Population)) +geom_point() + scale_colour_manual(values = KG_Palette,labels = relabel)+theme(legend.key.size = unit(0.5, "cm"), legend.position = "none")

f=ggplot(merge,aes(PC4,PC5, colour=Population)) +geom_point() + scale_colour_manual(values = KG_Palette, labels = relabel) +theme(legend.key.size = unit(0.5, "cm"), legend.position = "right")

#Add pc variance ####

eigenval <- fread("1KG.merged.eigenval")
pve <- data.frame(PC = 1:10, pve = (eigenval/sum(eigenval)*100))
pve$PC <- as.factor(pve$PC)

#Plot PVE ####
PC_Palette <- heat_hcl(10, h = c(320, 50), c. = c(35, 95), l = c(15, 90), power = c(0.8, 1.2), fixup = TRUE, gamma = NULL, alpha = 1)


PC_Palette <- heat_hcl(10, h = c(30, -160), c = c(80, NA, 45), l = c(38, 79), power = c(0.85, 1.0), fixup = TRUE, gamma = NULL, alpha = 1)

g <- ggplot(pve, aes(PC, V1, fill = PC)) + geom_bar(stat = "identity") + ylab("Proportion of Variance (%)\n") + scale_fill_manual(values = PC_Palette) + theme(legend.position = "none")

h <-ggplot(pve, aes(PC, V1, colour = PC)) +geom_line(group="PC", colour ="#454b51")+ geom_point(size =2.1)  + ylab("\nProportion of Variance (%)\n") + scale_colour_manual(values = PC_Palette) + theme(legend.position = "none")

layout <- c(
  area(1,1, 1,2),
  area(1,3, 1,4),
  area(1,5),
  area(2,1, 2,3)
)

#plot(layout)

#Save plots to pdf
pdf("enteric_pop_pca.pdf")
(a|b)/(c|d) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'a', title = 'Population structure of enteric CHIM samples')
e + f + guide_area() + h + plot_annotation(tag_levels = list(c('e',"f", "g")), title = 'Sub-Population structure of enteric CHIM samples') + plot_layout(design = layout, guides = 'collect')
dev.off()

save.image(file ="enteric_pop_pca.RData")

#PCA no IKG ####
pca2 <- fread("all_enteric_QC3.cleaned.eigenvec")
eigenval2 <- fread("all_enteric_QC3.cleaned.eigenval")

ggplot(data=pca2,aes(PC1, PC2)) + geom_point()
ggplot(data=pca2,aes(PC2, PC3)) + geom_point()
pve2 <- data.frame(PC = 1:10, pve = (eigenval2/sum(eigenval2)*100))
pve2$PC <- as.factor(pve2$PC)

ggplot(pve2, aes(PC, V1, colour = PC)) +geom_line(group="PC", colour ="#454b51")+ geom_point(size =2.1)  + ylab("\nProportion of Variance (%)\n") + scale_colour_manual(values = PC_Palette) + theme(legend.position = "none")

#Use this in assoc model

#colour pca by outcome  ####

pheno <- fread("~/GWAS_22/new_gwas/meta/final_metadata/all_pheno_num.txt")

merge2 <- full_join(pca, pheno, by = "IID")
merge2$Diagnosed<- as.factor(merge2$Diagnosed)
ggplot(data=merge2,aes(PC1, PC2, colour = Diagnosed)) + geom_point()
ggplot(data=merge2,aes(PC2, PC3, colour = Diagnosed)) + geom_point()

merge3 <- full_join(PCA, pheno, by = "IID")
merge3$Diagnosed<- as.factor(merge3$Diagnosed)

merge3 %>% filter(Diagnosed != "NA") %>% ggplot(aes(PC1, PC2, colour = Diagnosed)) + geom_point(alpha =0.9) + scale_colour_manual(values=c("seagreen3", "sandybrown"), labels = c("nTD", "TD"))
merge3 %>% filter(Diagnosed != "NA") %>% ggplot(aes(PC4, PC5, colour = Diagnosed)) + geom_point(alpha =0.9) + scale_colour_manual(values=c("seagreen3", "sandybrown"), labels = c("nTD", "TD"))

x=merge2 %>% filter(Diagnosed != "NA") %>% ggplot(aes(PC1, PC2, colour = Diagnosed)) + geom_point(alpha =0.9) + scale_colour_manual(values=c("seagreen3", "sandybrown"), labels = c("nTD", "TD"))
y=merge2 %>% filter(Diagnosed != "NA") %>% ggplot(aes(PC2, PC3, colour = Diagnosed)) + geom_point(alpha =0.9) + scale_colour_manual(values=c("seagreen3", "sandybrown"), labels = c("nTD", "TD"))
z=merge2 %>% filter(Diagnosed != "NA") %>% ggplot(aes(PC3, PC4, colour = Diagnosed)) + geom_point(alpha =0.9) + scale_colour_manual(values=c("seagreen3", "sandybrown"), labels = c("nTD", "TD"))
za=merge2 %>% filter(Diagnosed != "NA") %>% ggplot(aes(PC4, PC5, colour = Diagnosed)) + geom_point(alpha=0.9) + scale_colour_manual(values=c("seagreen3", "sandybrown"), labels = c("nTD", "TD"))
                                                                                         
(x|y)/(z|za)+plot_layout(guides = "collect") +plot_annotation(tag_levels = "a", title = "PCA by Diagnosis")                                         

#Check covars