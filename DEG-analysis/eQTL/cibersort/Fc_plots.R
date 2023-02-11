setwd(R)
#load packages ####
library(dplyr)
library(data.table)
library(tidyr)
library(tidylog)
library(limma)
library(stringi)
library(janitor)
library(stringr)
library(splines)
library(patchwork)

#Plots ####

#1) Volcano Plots
#2) FcR expression/splines plots
#3) FcR box plots 
#4) Genotype Plots
#5) Cis-EQTL plots

#2) FcR Splines ####

FC1B <- Fc_exprs %>% ggplot(aes(x=time, y=fcgr1b, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 14.3) +
  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
  geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 5)) +xlab("") + ggtitle(expression(paste("FC",gamma,"R1",beta))) + ylab("")

FC2A<- Fc_exprs %>% ggplot(aes(x=time, y=fcgr2a, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 14.3) + ylim(10,13.3) +
  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
  geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 3)) +
  xlab("") + ggtitle(expression(paste("FC",gamma,"R2",alpha))) + ylab("")

FC3A <- Fc_exprs %>% ggplot(aes(x=time, y=fcgr3a, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 14.3) +
  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
  geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 4)) + 
  xlab("Time (days)") + ggtitle(expression(paste("FC",gamma,"R3",alpha))) + ylab("Expression\n")
#Not much interesting with 3a


                             
FCAR <- Fc_exprs %>% ggplot(aes(x=time, y=fcar, colour = diagnosis)) +
                               geom_point() + 
                               xlim(0, 14.3) +
                               scale_colour_manual(values=c("seagreen3", "sandybrown")) +
                               geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 4)) + xlab("") + ggtitle(expression(paste("FC",alpha,"R"))) + ylab("Expression\n")


                             
FCRLA <- Fc_exprs %>% ggplot(aes(x=time, y=fcrla, colour = diagnosis)) +
                               geom_point() + 
                               xlim(0, 14.3) +
                               scale_colour_manual(values=c("seagreen3", "sandybrown")) +
                               geom_smooth( method = "lm",aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 3)) + 
                               xlab("Time (days)") + ggtitle(expression(paste("FCRL",alpha))) + ylab("")
                             
 FCRLB <- Fc_exprs %>% ggplot(aes(x=time, y=fcrlb, colour = diagnosis)) +
                  geom_point() + 
                  xlim(0, 14.3) +
                  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
                  geom_smooth(method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 3)) + xlab("Time (days)") + ggtitle(expression(paste("FCRL",beta))) + ylab("")
  



List <- list(FCAR, FC1B, FC2A, FC3A, FCRLA, FCRLB) 

#Plot
Plot <- wrap_plots(List,ncol = 3,nrow = 2)

Plot + plot_annotation(tag_levels = 'A')  +
  plot_layout(guides = 'collect')

#Size W1242, H700, TIFF

#Fc Volcano Plots ####

#Requires setting up separate contrast matrices for each comparison
#I'm also subsetting the data as per comparison beforehand but not sure if that's necessary
#remove nas for any analysis variables
#Do at this stage so samples from both the pheno and exprs rows with na's are removed
pheno_exprs <- pheno_exprs %>% drop_na(time) %>% filter(time >=0)

#Contrast filter (change to suit contrasts)
filtered <- filter(pheno_exprs, time == "0"| time == "0.5" | time == "1"| time == "7"| time == "8"|timepoint3 == "TD")
#filter(pheno_exprs, time == "0"| time == "0.5")
pData2 <- filtered[,1:14]
exprs2 <- filtered[,15:10461]
#Reformat expression data
exprs2 <- sapply(exprs2, as.numeric)
exprs2 <- t(exprs2) #185 cols, 10447 rows

#The number of rows in phenoData must match the number of columns in assayData
#Row names of phenoData must match column names of the matrix / matrices in assayData.
#Make Time_ID row names in pheno data
rownames(pData2) <- pData2$time_id
colnames(exprs2) <- pData2$time_id

#2) Design Matrices ####
table(pData2$time) 
#Set up vars for design matrix ####
#make time as numeric if modelling across time or factor for direct comps
#The design matrix tells the linear model if a particular condition is on, therefore to calculate/model the expression for that value
#So if comparing TD vs nTD the on off switch will be for that (pre-select timepoint)
#If comparing expression between time points then the on off switch is needed for that too
pData2$time <- as.factor(pData2$time)
pData2$timepoint3 <- as.factor(pData2$timepoint3)
pData2$study_arm <- as.factor(pData2$study_arm)
#make new vaccine factor (will give more df)
pData2$vax <- ifelse(pData2$study_arm == "T1",0,
                     ifelse(pData2$study_arm == "Placebo",0,
                            ifelse(pData2$study_arm == "M01ZH09", 1,
                                   ifelse(pData2$study_arm == "Ty21a", 1,NA))))
pData2$vax <- as.factor(pData2$vax)

design_full <- model.matrix(~0 + time + diagnosis + age + sex + vax, data = pData2)
design_TD <- model.matrix(~0 + timepoint3 + diagnosis + age + sex + vax, data = pData2)
#Tidy names
colnames(design)<-gsub("vax","",colnames(design_full))
colnames(design)<-gsub("timepoint3","",colnames(design_full))
colnames(design)<-gsub("vax","",colnames(design_TD))
colnames(design)<-gsub("timepoint3","",colnames(design_TD))


#Blocking var ####
#Then we estimate the correlation between measurements made on the same subject:
corfit <- duplicateCorrelation(exprs2,design=design_full,block=pData2$part_number)
corfit_TD <- duplicateCorrelation(exprs2,design=design_full,block=pData2$part_number)
corfit$consensus #0.347
#Then this inter-subject correlation is input into the linear model fit:
fit <- lmFit(exprs2, design=design_full, block=pData2$part_number, correlation=corfit$consensus)

#Contrast matrices - Day 0 vs Day 1 ####
#Now we can make any comparisons between the experimental conditions in the usual way, example:
cm <- makeContrasts(time1-time0, levels = colnames(design_full))

ebayesfit <- eBayes(contrasts.fit(lmFit(exprs2,design=design_full,block=pData2$part_number,correlation=corfit$consensus), cm))
results <- topTable((ebayesfit),num = Inf)

#Volcano Plot ####
ggplot(data=results, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point() + theme_minimal()

deg <- results
#Genes of interest to highlight
gene_list <- c("fcgr3a","fcar","fcgr1b","fcrla","fcrlb","fcgr2a")
#grep labels on row.names
deg$gene_name<-row.names(deg) 
#Make data table with absolute FC values of genes of interest, in this case FCR genes
data= subset(deg, rownames(deg) %in% gene_list)

#For main data mutate a new variable, reg, if FC and P values are above/below a certain threshold
deg <- deg %>%
  mutate(reg = case_when(
    deg$logFC >= 0 & deg$adj.P.Val <= 0.05 ~ "UP",
    deg$logFC <= 0 & deg$adj.P.Val <= 0.05 ~ "DOWN",
    abs(deg$logFC) <= 0 & deg$adj.P.Val >= 0.05 ~ "no_change",
    abs(deg$logFC) <= 0 & deg$adj.P.Val <= 0.05 ~ "no_change",
    abs(deg$logFC) > 0 & deg$adj.P.Val >0.05 ~ "no_change"
  )) %>%
  mutate(reg = factor(reg, levels = c("UP", "no_change","DOWN")))

#Plot volcano plot 
deg %>% ggplot(aes(x=logFC,y=-log10(P.Value)))+ geom_point(aes(color=reg))
               
#Volcano Plot Labelled ####
#1 Day
D0_D1 <- deg %>% ggplot(aes(x=logFC, y=-log10(adj.P.Val),label=gene_name))+ geom_point(aes(color=adj.P.Val))+ scale_color_gradientn(colours = c("darkred","#a5342d", "orange", "yellow"), values=c(0,0.01, 0.05, 0.1,1)) + theme_light() + labs(title = "Baseline - 24h Post-Challenge")+
  geom_label_repel(data=data,size=4,direction="both",nudge_y =3,nudge_x =-0.1,angle= 60,vjust= 0,segment.size= 0.5,segment.color="black",fill="grey") + xlim(-1,1.5) + geom_hline(yintercept=1.3, linetype="dashed", color = "#3b3a39", size = 0.5) + theme(legend.position = "none")

D0_D1


#Day 0 vs TD ####
design_TD <- model.matrix(~0 + timepoint3 + diagnosis + age + sex + vax, data = pData2)
#Tidy names
colnames(design_TD)<-gsub("vax","",colnames(design_TD))
colnames(design_TD)<-gsub("timepoint3","",colnames(design_TD))
colnames(design_TD)<-gsub("\\+","plus",colnames(design_TD))
colnames(design_TD)<-gsub("\\-","minus",colnames(design_TD))
#Do I have to remove the other columns I'm not using? I will so they're not included in model
design_TD <- design_TD[,c(1,9,13:15)]
#Estimate cor
corfit_TD <- duplicateCorrelation(exprs2,design=design_TD,block=pData2$part_number)
corfit$consensus #0.221
#Then this inter-subject correlation is input into the linear model fit:
fit <- lmFit(exprs2, design=design_TD, block=pData2$part_number, correlation=corfit_TD$consensus)
#Contrast matrices ####
#Now we can make any comparisons between the experimental conditions in the usual way, example:
cm_TD <- makeContrasts(TD-D0, levels = colnames(design_TD))
ebayesfit <- eBayes(contrasts.fit(lmFit(exprs2,design=design_TD,block=pData2$part_number,correlation=corfit$consensus), cm_TD))
results <- topTable((ebayesfit),num = Inf)

#Volcano set up ####
deg <- results
#Genes of interest to highlight
gene_list <- c("fcgr3a","fcar","fcgr1b","fcrla","fcrlb","fcgr2a")
#grep labels on row.names
deg$gene_name<-row.names(deg) 
#Make data table with absolute FC values of genes of interest, in this case FCR genes
data= subset(deg, rownames(deg) %in% gene_list)
#For main data mutate a new variable, reg, if FC and P values are above/below a certain threshold
deg <- deg %>%
  mutate(reg = case_when(
    deg$logFC >= 0 & deg$adj.P.Val <= 0.05 ~ "UP",
    deg$logFC <= 0 & deg$adj.P.Val <= 0.05 ~ "DOWN",
    abs(deg$logFC) <= 0 & deg$adj.P.Val >= 0.05 ~ "no_change",
    abs(deg$logFC) <= 0 & deg$adj.P.Val <= 0.05 ~ "no_change",
    abs(deg$logFC) > 0 & deg$adj.P.Val >0.05 ~ "no_change"
  )) %>%
  mutate(reg = factor(reg, levels = c("UP", "no_change","DOWN")))
#Plot volcano plot 
deg %>% ggplot(aes(x=logFC,y=-log10(P.Value)))+ geom_point(aes(color=reg))
#Volcano Plot Labelled ####
#Day 0 -TD
TD <- deg %>% ggplot(aes(x=logFC, y=-log10(adj.P.Val),label=gene_name))+ geom_point(aes(color=adj.P.Val))+ scale_color_gradientn(colours = c("darkred","#a5342d", "orange", "yellow"), values=c(0,0.01, 0.05, 0.1,1)) + theme_light() + labs(title = "Baseline - Day of Diagnosis")+
  geom_label_repel(data=data,size=4,direction="both",nudge_y =4,nudge_x =-0.15,angle= 60,vjust= 0,segment.size= 0.5,segment.color="black",fill="grey") + geom_hline(yintercept=1.3, linetype="dashed", color = "#3b3a39", size = 0.5) +theme(legend.position = "none")
TD

top <-topTable(ebayesfit, num = 100)
write.csv(top, file = "top_table_T1T2_D0_TD.csv")


#nTD vs TD D1 #####
#Filter so just day 1 time point then make contrasts between diagnosis
#In design filter
design_D1<- model.matrix(~0 + diagnosis +timepoint3 + age + sex + vax, data = pData2)
#Tidy names
colnames(design_D1)<-gsub("timepoint3","",colnames(design_D1))
colnames(design_D1)<-gsub("diagnosis","",colnames(design_D1))
colnames(design_D1)<-gsub("\\+","plus",colnames(design_D1))
colnames(design_D1)<-gsub("\\-","minus",colnames(design_D1))
#Do I have to remove the other columns I'm not using? I will so they're not included in model
design_D1 <- design_D1[,c(1,2,4,14:16)]
#Estimate cor
corfit_D1 <- duplicateCorrelation(exprs2,design=design_D1,block=pData2$part_number)
corfit$consensus #0.221
#Then this inter-subject correlation is input into the linear model fit:
fit <- lmFit(exprs2, design=design_D1, block=pData2$part_number, correlation=corfit_D1$consensus)
#Contrast matrices ####
#Now we can make any comparisons between the experimental conditions in the usual way, example:
cm_D1 <- makeContrasts(TD-nTD, levels = colnames(design_D1))
ebayesfit <- eBayes(contrasts.fit(lmFit(exprs2,design=design_D1,block=pData2$part_number,correlation=corfit$consensus), cm_D1))
results <- topTable((ebayesfit),num = Inf)

#Volcano set up ####
deg <- results
#Genes of interest to highlight
gene_list <- c("fcgr3a","fcar","fcgr1b","fcrla","fcrlb","fcgr2a")
#grep labels on row.names
deg$gene_name<-row.names(deg) 
#Make data table with absolute FC values of genes of interest, in this case FCR genes
data= subset(deg, rownames(deg) %in% gene_list)
#For main data mutate a new variable, reg, if FC and P values are above/below a certain threshold
deg <- deg %>%
  mutate(reg = case_when(
    deg$logFC >= 0 & deg$adj.P.Val <= 0.05 ~ "UP",
    deg$logFC <= 0 & deg$adj.P.Val <= 0.05 ~ "DOWN",
    abs(deg$logFC) <= 0 & deg$adj.P.Val >= 0.05 ~ "no_change",
    abs(deg$logFC) <= 0 & deg$adj.P.Val <= 0.05 ~ "no_change",
    abs(deg$logFC) > 0 & deg$adj.P.Val >0.05 ~ "no_change"
  )) %>%
  mutate(reg = factor(reg, levels = c("UP", "no_change","DOWN")))
#Plot volcano plot 
deg %>% ggplot(aes(x=logFC,y=-log10(P.Value)))+ geom_point(aes(color=reg))
#Volcano Plot Labelled ####
#Day 0 -D1
D1 <- deg %>% ggplot(aes(x=logFC, y=-log10(adj.P.Val),label=gene_name))+ geom_point(aes(color=adj.P.Val))+ scale_color_gradientn(colours = c("darkred","#a5342d", "orange", "yellow"), values=c(0,0.01, 0.05, 0.1,1)) + theme_light() + labs(title = "nTD - TD (1 Day post-challenge)")+
  geom_label_repel(data=data,size=4,direction="y",nudge_y =4,nudge_x =-0.15,angle= 60,vjust= 0,segment.size= 0.5,segment.color="black",fill="grey") + geom_hline(yintercept=1.3, linetype="dashed", color = "#3b3a39", size = 0.5) + theme(legend.position = "none")
D1

#ntD vs TD @ Day 7 ####

#Weird error time 0 has disappeared from matrix even though it's in the pData2
#Error in eval(ej, envir = levelsenv) : object 'D0' not found

table(pData2$time)


design_D7<- model.matrix(~0 +diagnosis + age + sex + vax +time, data = pData2)
#make new vaccine factor (will give more df)
D0 <- ifelse(pData2$timepoint3 == "D0",1,0)

design_D7 <- cbind(design_D7, D0)
#Tidy names
colnames(design_D7)<-gsub("diagnosis","",colnames(design_D7))

#Do I have to remove the other columns I'm not using? I will so they're not included in model
design_D7 <- design_D7[,c(1:5,10,20)]
#Estimate cor
corfit_D7 <- duplicateCorrelation(exprs2,design=design_D7,block=pData2$part_number)
corfit$consensus #0.221
#Then this inter-subject correlation is input into the linear model fit:
fit <- lmFit(exprs2, design=design_D7, block=pData2$part_number, correlation=corfit_D7$consensus)
#Contrast matrices ####
#Now we can make any comparisons between the experimental conditions in the usual way, example:
cm_D7 <- makeContrasts(time7-D0, levels = colnames(design_D7))

ebayesfit <- eBayes(contrasts.fit(lmFit(exprs2,design=design_D7,block=pData2$part_number,correlation=corfit$consensus), cm_D7))
results <- topTable((ebayesfit),num = Inf)

#Volcano set up ####
deg <- results
#Genes of interest to highlight
gene_list <- c("fcgr3a","fcar","fcgr1b","fcrla","fcrlb","fcgr2a")
#grep labels on row.names
deg$gene_name<-row.names(deg) 
#Make data table with absolute FC values of genes of interest, in this case FCR genes
data= subset(deg, rownames(deg) %in% gene_list)
#For main data mutate a new variable, reg, if FC and P values are above/below a certain threshold
deg <- deg %>%
  mutate(reg = case_when(
    deg$logFC >= 0 & deg$adj.P.Val <= 0.05 ~ "UP",
    deg$logFC <= 0 & deg$adj.P.Val <= 0.05 ~ "DOWN",
    abs(deg$logFC) <= 0 & deg$adj.P.Val >= 0.05 ~ "no_change",
    abs(deg$logFC) <= 0 & deg$adj.P.Val <= 0.05 ~ "no_change",
    abs(deg$logFC) > 0 & deg$adj.P.Val >0.05 ~ "no_change"
  )) %>%
  mutate(reg = factor(reg, levels = c("UP", "no_change","DOWN")))
#Plot volcano plot 
deg %>% ggplot(aes(x=logFC,y=-log10(P.Value)))+ geom_point(aes(color=reg))
#Volcano Plot Labelled ####
#Day 0 -D7
D7 <- deg %>% ggplot(aes(x=logFC, y=-log10(adj.P.Val),label=gene_name))+ geom_point(aes(color=adj.P.Val))+ scale_color_gradientn(colours = c("darkred","#a5342d", "orange", "yellow"), values=c(0,0.001, 0.05, 0.1,1)) + theme_light() + labs(title = "nTD - TD (7 Days post-challenge)")+ geom_label_repel(data=data,size=4,direction="y",nudge_y =4,nudge_x =-0.15,angle= 60,vjust= 0,segment.size= 0.5,segment.color="black",fill="grey") + geom_hline(yintercept=1.3, linetype="dashed", color = "#3b3a39", size = 0.5)
D7

#Group Plot ####
List <- list(D0_D1, TD, D1, D7)
Plot <- wrap_plots(List,ncol = 2,nrow = 2)

Plot + plot_annotation(tag_levels = 'A')  +
  plot_layout(guides = 'collect')

#Could also do baseline - Day 7 time
               