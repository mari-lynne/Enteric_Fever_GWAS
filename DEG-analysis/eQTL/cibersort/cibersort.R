# Cibersort ####

# Filter time expression data -------------------------------------------------
setwd("~/GWAS_22/gwas_final/eQTL")

# RNAseq data 
# Load micro-array experiment data
load("T1_T2_with_rsn.norm.R")
data = T1_T2_autosomes
rm(T1_T2_autosomes)

exprs <- as.data.frame(exprs(data))
fData <- clean_names(data@featureData@data)
pData <- clean_names(setDT(data@phenoData@data, keep.rownames = TRUE)[])

## Set up new vars -------------------------------------------------------------
# Make new typhoid diagnosis variable, na's = 0
pData <- pData %>%
  mutate(diagnosis =
           ifelse(is.na(pData$day_of_typhoid_diagnosis_from_challenge),
                  0, 1))
# Rename vars
names(pData)[names(pData) == "days_since_challenge"] <- "time"
names(pData)[names(pData) == "timepoint3"] <- "visit"

# Time point
table(pData$time, pData$study_arm)
pData$visit <- str_replace_all(pData$visit, "D0\\+12", "D0.12h")
# Make 12_24h time point
pData <- pData %>% mutate(time_point = ifelse(visit %in% c("D0.12h","D1"), "12-24h", visit))

### Make new baseline time point ----------------------
# Make day of challenge V28 for vaccine samples day 0
pData <- pData %>% 
  mutate(time_point = ifelse(time_point == "V0", "V28", time_point))
# Make new baseline time point using V0 for vaccinees and D0 for challenge
pData <- pData %>% 
  mutate(time_point = ifelse(time_point %in% c("V0","D0"), "Baseline", time_point)) 
table(pData$time_point, pData$study_arm)

### New day 7 time point -------------------------------
pData <- pData %>% 
  mutate(time_point = ifelse(time %in% c("6.5","7","8","9") & str_detect(visit, "TD"), "D7", time_point))
pData <- pData %>% 
  mutate(time_point = ifelse(time_point == "D9" & study_arm == "M01ZH09", "D7", time_point))
pData <- pData %>% 
  mutate(time_point = ifelse(time_point == "D8" & study_arm == "Ty21a", "D7", time_point))

# Add exprs colnames to pData for matching tables by
pData$exprs_cols <- colnames(exprs)
# Dups
dup <- pData[!stri_duplicated(pData$part_number) & time_point == "D7",]

# Filter data with dup$part_number and a D7 time point

# 1) Filter DGE data for time point and genotyped samples ----------------------

pData <- pData %>% filter(time_point %in% c("Baseline", "12-24h", "D7"))
exprs <- exprs[,colnames(exprs) %in% pData$exprs_cols]

# Prep mixture file ------------------------------------------------------------
# Make expression gene table
fData$chromosome_name <- as.numeric(fData$chromosome_name)
fData <- fData %>%
  filter(chromosome_name %in% c(1:22) & !is.na(ilmn_gene))
# Filter exprs
exprs <- filter(exprs, row.names(exprs) %in% fData$ensembl_gene_id)

# Update gene names to symbol
exprs$ensembl_gene_id <- rownames(exprs)
exprs <- left_join(fData, exprs, by = "ensembl_gene_id")
# clean gene names
rownames(exprs) <- make_clean_names(exprs$ilmn_gene, case = "all_caps")
exprs <- exprs[,-c(1:30)]

# If no fData
# library(biomaRt)
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl)
# # Only assoc gene loc data with genes that are in our expression data
# filtered_genes <- genes[genes$ensembl_gene_id %in% rownames(exprs),]
# filtered_genes <- filtered_genes[(filtered_genes$chromosome_name %in% c(1:22)),] 
# # Fill in NAs
# filtered_genes[filtered_genes == ""]<- NA
# # filter na genes
# filtered_genes <- filter(filtered_genes, !is.na(hgnc_symbol))


  # Filter time points -----------------------------
pData$time_id <- str_c(pData$part_number, pData$time_point, pData$study_arm, sep = "_")
colnames(exprs) <- pData$time_id 

# Combat (adjust for batch effects) -----------------
exprs <- as.matrix(exprs)
colnames(combat)
mod = model.matrix(~as.factor(time_point), data=pData)

combat <- as.data.frame(sva::ComBat(dat=exprs,
                          batch=pData$array_experiment,
                          mod=mod))

combat <- tibble::rownames_to_column(combat, "Gene")

D0 <- combat[ ,grep("Baseline|Gene", colnames(combat))]
D1 <- combat[ ,grep("12-24h|Gene", colnames(combat))]
D7 <- combat[ ,grep("D7|Gene", colnames(combat))]

# Make tables
# In notepad add tab to header and write Gene as colname 
write.table(D0, file = "T1T2/D0_cibersort.txt", sep = "\t", quote = FALSE, row.names = F)
write.table(D1, file = "T1T2/D1_cibersort.txt", sep = "\t", quote = FALSE, row.names = F)
write.table(D7, file = "T1T2/D7_cibersort.txt", sep = "\t", quote = FALSE, row.names = F)

check <- fread("T1T2/T1T2_D0_cibersort.txt")

# Results ---------------------------------------------------------------------

# Cell -specific --------------------------------------------------------------
neut <- fread("~/GWAS_22/gwas_final/eQTL/T1T2/cibersort/D7/CIBERSORTxHiRes_Job18_Neutrophils_Window40.txt")
neut1 <- fread("T1T2/cibersort/D1/CIBERSORTxHiRes_Job19_Neutrophils_Window40.txt")
neut0 <- fread("T1T2/cibersort/D0/CIBERSORTxHiRes_Job20_Neutrophils_Window40.txt")


# Join tables
neut <- melt(neut, id.vars = c("GeneSymbol"), measure.vars = c(2:ncol(neut)))
colnames(neut) <- c("Gene","time_id", "Expression")
neut$cell_type <- rep("Neutrophil", nrow(neut))
neut$time_point <- rep("D7", nrow(neut))

neut1 <- melt(neut1, id.vars = c("GeneSymbol"), measure.vars = c(2:ncol(neut1)))
colnames(neut1) <- c("Gene","time_id", "Expression")
neut1$cell_type <- rep("Neutrophil", nrow(neut1))
neut1$time_point <- rep("D1", nrow(neut1))

neut0 <- melt(neut0, id.vars = c("GeneSymbol"), measure.vars = c(2:ncol(neut0)))
colnames(neut0) <- c("Gene","time_id", "Expression")
neut0$cell_type <- rep("Neutrophil", nrow(neut0))
neut0$time_point <- rep("Baseline", nrow(neut0))

neut <- bind_rows(neut0, neut1, neut)

# Plot 
neut1 %>% group_by(Gene) %>% summarise(avg = mean(Expression)) %>%
  ggplot(aes(x = Gene, y = avg, fill = Gene)) +
  geom_bar(width = 1, stat = "identity") +
  labs(title = "D7 - Neutrophil", x = "\nFc-Receptor", y = "Gene Expression (average)\n") +
  scale_fill_viridis_d() + theme_pubr() +
  theme(legend.position = "none")

neut  %>% group_by(Gene, time_point) %>% summarise(avg = mean(Expression)) %>%
  ggplot(aes(x = Gene, y = avg, fill = Gene)) +
  geom_bar(width = 1, stat = "identity") +
  labs(title = "Neutrophil FcR-Expression", x = "\nFc-Receptor", y = "Gene Expression (average)\n") +
  scale_fill_viridis_d() +
  theme(legend.position = "none") + facet_wrap(~time_point)

# Calculate fold change 
neut$id <- str_extract(neut$time_id, "\\d*(?=_)")
# test <- neut %>% group_by(id, time_point) %>% filter(nrow(id) == 3)
#group_by(id) %>%
#  mutate(FC =
           (Expression[time_point != "Baseline"] - Expression[time_point == "Baseline"])
         / (Expression[time_point == "Baseline"]))


# #neut  %>% filter(Gene == "FCAR") %>% filter(str_detect(time_id, "TD")) %>% group_by(time_point) %>% summarise(avg = mean(Expression)) %>%
#   ggplot(aes(x = time_point, y = avg, fill = time_point)) +
#   geom_bar(width = 1, stat = "identity") +
#   labs(title = "Neutrophil FcR-Expression", x = "\nFc-Receptor", y = "Gene Expression (average)\n") +
#   scale_fill_viridis_d() +
#   theme(legend.position = "none")



?facet_grid


# Bulk ------------------------------------------------------------------------
ciber_0 <- fread("~/RNA/cibersort/CIBERSORTx_Day0/CIBERSORTxGEP_Job13_GEPs.txt")
ciber_1 <- fread("~/RNA/cibersort/CIBERSORTxGEP_Job9_GEPs.txt")
ciber_1 <- fread("~/RNA/cibersort/CIBERSORTx_Day1/CIBERSORTxGEP_Job11_GEPs.txt")
ciber_7 <- fread("~/RNA/cibersort/CIBERSORTx_Day7/CIBERSORTxGEP_Job10_GEPs.txt")

gene_list <- c("FCGR3A","FCAR","FCGR2A","FCGR1B", "FCGR2C")

'%!in%' <- function(x,y)!('%in%'(x,y))

## Baseline -------------------------------------------------------------------

# extra subsets ciber_0 <- fread("cibersort/CIBERSORTx_Day0/CIBERSORTxGEP_Job12_GEPs.txt")
ciber_0 = subset(ciber_0, GeneSymbol %in% gene_list)

# Reform df for ggplot
ciber_0 <- melt(ciber_0, id.vars = c("GeneSymbol"), measure.vars = c(2:11))
colnames(ciber_0) <- c("FcR","Cell.Type", "Expression")

bp_D0 <- ciber_0 %>%
  filter(Cell.Type %!in% c("Eosinophils","Mast cells")) %>%
  ggplot(aes(x = FcR, y = Expression, fill = Cell.Type)) +
  geom_bar(width = 1, stat = "identity") +
  labs(title = "Baseline", x = "\nFcR") +
  ylim(0, 60000) +
  scale_fill_viridis_d() + theme_pubr() +
  theme(legend.position = "none")

bp_D0

## 12h and 24hs -----------------------------------------------------------

# Shows the relative expression of FcRs on those cell types
ciber_1 = subset(ciber_1, GeneSymbol %in% gene_list)

#Bar Chart
ciber_1 <- melt(ciber_1, id.vars = c("GeneSymbol"), measure.vars = c(2:11))
colnames(ciber_1) <- c("FcR","Cell.Type", "Expression")

bp_1 <- ciber_1 %>%
  filter(Cell.Type %!in% c("Eosinophils","Mast cells")) %>%
  ggplot(aes(x = FcR, y = Expression, fill = Cell.Type)) +
  geom_bar(width = 1, stat = "identity") +
  labs(title = "Day 1", x = "\nFcR") +
  ylim(0, 60000)  +
  scale_fill_viridis_d() + theme_pubr() +
  theme(legend.position = "right")

bp_1

## Results Day 7 -----------------------------------------------------------

# Could do a day7/8 combined
ciber_7 = subset(ciber_7, GeneSymbol %in% gene_list)

#Bar Chart
ciber_7 <- melt(ciber_7, id.vars = c("GeneSymbol"), measure.vars = c(2:11))
colnames(ciber_7) <- c("FcR","Cell.Type", "Expression")

bp_7 <- ciber_7 %>%
  filter(Cell.Type %!in% c("Eosinophils","Mast cells")) %>%
  ggplot(aes(x = FcR, y = Expression, fill = Cell.Type)) +
  geom_bar(width = 1, stat = "identity") +
  labs(title = "Day 7", x = "\nFcR") +
  ylim(0, 80000) +
  scale_fill_viridis_d() + theme_pubr() +
  theme(legend.position = "right")


bp_7

### Group plots -----------------------------------

(bp_D0 |bp_1) +
  plot_annotation(tag_levels = 'a', title = "FcR Expression") +
  plot_layout(guides = 'collect')



#1600, 750


# Original Results ####
og_ciber <- fread("CiberSortInput.txt")
rownames(og_ciber) <- og_ciber$Gene
og_ciber <-as.matrix(og_ciber)
exprs_1 <- og_ciber[,grep("V1_24", colnames(og_ciber))] #Just 12h data

# No eosinophil/mast cells
bp_1 <- ciber_1 %>% filter(Cell.Type != "Mast cells",Cell.Type != "Eosinophils") %>%
  ggplot(aes(x=FcR, y=Expression, fill=Cell.Type))+
  geom_bar(width = 1, stat = "identity") + labs(title = "Day 1\n", x="\nFcR", y="") + ylim(0,60000) + theme(legend.position = "none", axis.text.x = element_text(size= 11,face="bold")) +scale_fill_manual(values=wes)

bp_1

bp_D0<- ciber_0 %>% filter(Cell.Type != "Mast cells",Cell.Type != "Eosinophils") %>%
  ggplot(aes(x=FcR, y=Expression,fill=Cell.Type))+
  geom_bar(width = 1, stat = "identity") + labs(title = "Baseline\n", x="\nFcR") + ylim(0,60000) + theme(legend.position = "none", axis.text.x = element_text(size= 11,face="bold")) +scale_fill_manual(values=wes)

#+scale_fill_brewer(palette = "Accent")
wes <- c("#46AC8C", "#0B775E","#C6CDF7","#7294d4","#FD6467","#EBCC2A","#E58601", "#B40F20")

bp_D0

bp_7 <-  ciber2 %>% filter(Cell.Type != "Mast cells",Cell.Type != "Eosinophils") %>%
  ggplot(aes(x=FcR, y=Expression, fill=Cell.Type))+
  geom_bar(width = 1, stat = "identity") + labs(title = "Day 7\n",x="\nFcR", y="") + ylim(0,60000) + theme(axis.text.x = element_text(size= 11,face="bold"))+scale_fill_manual(values=wes)

bp_7


list <- list(bp_D0, bp_1, bp_7)
Plot <- wrap_plots(list,ncol = 3,nrow = 1)

Plot + plot_annotation(tag_levels = 'A')  +
  plot_layout(guides = 'collect')
