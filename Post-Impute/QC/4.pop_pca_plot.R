# Date: 2022/06/02
# Purpose: Plot user data projected on 1KG PCs, greying out 1KG individuals

# Updates 12/10/22
# Reran PCA of study data after pruning
# Have separated into just typhoid cohort ( n=253)

# Load data --------------------------------------------------------------------

setwd("/home/mari/GWAS_22/gwas_final/merge/enteric/QC/pca")

# Typhoid: typhoid2
# Enteric: enteric.IBD

# Load packages

library(dplyr)
library(tidylog)
library(patchwork)
library(ggplot2)
library(colorspace)
library(readr)
library(data.table)
library(forcats)
library(viridis)


# Initialise command line arguments
args <- commandArgs(TRUE)

# Load data
root <- args[1]

PCA <- read_delim(
  "all_hg38-typhoid2.LD.eigenvec",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Rename columns
colnames(PCA) <-
  c("FID","IID", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

# Merge with Population data
pop <- read_delim(
  "1kG.ID2Pop.txt",
  delim = " ",
  escape_double = FALSE,
  trim_ws = TRUE
)

pop <- pop %>% rename(IID = `#IID`)

merge <- left_join(PCA, pop, by = "IID")

# Plot ancestry set up ---------------------------------------------------------

# Rename NAs in Pop data as Study
merge$SuperPop <-
  merge$SuperPop %>%
  replace_na("Study") 

merge$Population <-
  merge$Population %>%
  replace_na("Study")

# Recode to factors
merge$SuperPop <-
  as.factor(merge$SuperPop)
levels(merge$SuperPop)

merge$Population <-
  as.factor(merge$Population)
levels(merge$Population)

# Order labels for plotting
merge$SuperPop <- fct_relevel(merge$SuperPop, "Study", after = Inf)
merge$Population <- fct_relevel(merge$Population, "Study", after = Inf)

merge <- merge[order(merge$SuperPop), ]
merge <- merge[order(merge$Population), ]

# Population  PCA Plots --------------------------------------------------------

# Define colour palette for populations
# Super pop

KG_Palette_Super <-
  heat_hcl(
    length(unique(merge$SuperPop)),
    h = c(300, 75),
    c. = c(35, 95),
    l = c(15, 90),
    power = c(0.8, 1.2),
    fixup = TRUE,
    gamma = NULL,
    alpha = 0.85
  )

KG_Palette <-
  heat_hcl(
    length(unique(merge$Population)),
    h = c(300, 75),
    c. = c(35, 95),
    l = c(15, 90),
    power = c(0.8, 1.2),
    fixup = TRUE,
    gamma = NULL,
    alpha = 0.8
  )

a <- 
  ggplot(merge,
       aes(PC1, PC2, colour = SuperPop)) +
  geom_point() +
  scale_colour_manual(values = (KG_Palette_Super)) +
  theme_bw() +
  ggtitle("Typhoid GWAS Population Stratification (1KG)")
a

# Extra PCs Superpop

b <-
  ggplot(merge, aes(PC2, PC3, colour = SuperPop)) +
  geom_point() +
  scale_colour_manual(values = KG_Palette_Super) +
  theme_bw() 

b
c <- 
  ggplot(merge, aes(PC3, PC4, colour = SuperPop)) +
  geom_point() + scale_colour_manual(values = KG_Palette_Super) +
  theme_bw()

d <- 
  ggplot(merge, aes(PC4, PC5, colour = SuperPop)) +
  geom_point() + scale_colour_manual(values = KG_Palette_Super) +
  theme_bw()

e <- 
  ggplot(merge, aes(PC5, PC6, colour = SuperPop)) +
  geom_point() + scale_colour_manual(values = KG_Palette_Super) +
  theme_bw()

f <- 
  ggplot(merge, aes(PC6, PC7, colour = SuperPop)) +
  geom_point() + scale_colour_manual(values = KG_Palette_Super) +
  theme_bw()

g <- 
  ggplot(merge, aes(PC7, PC8, colour = SuperPop)) +
  geom_point() + scale_colour_manual(values = KG_Palette_Super) +
  theme_bw()

h <- 
  ggplot(merge, aes(PC8, PC9, colour = SuperPop)) +
  geom_point() + scale_colour_manual(values = KG_Palette_Super) +
  theme_bw()


i <- ggplot(merge, aes(PC9, PC10, colour = SuperPop)) +
  geom_point() + scale_colour_manual(values = KG_Palette_Super) +
  theme_bw() 

i

pdf("all_ancestry_PCAs.pdf")
a
b
c
d
e
f
g
h
i
dev.off()


# Sub population PCA

pdf("KG.subpop_strat_PCA.pdf")
y <- 
  ggplot(merge,
         aes(PC1, PC2, colour = Population)) +
  geom_point() +
  scale_colour_manual(values = KG_Palette) +
  theme_bw()
y
dev.off()



# Add pc variance --------------------------------------------------------------

eigenval <- fread("all_hg38-typhoid2.LD.eigenval")
pve <- data.frame(PC = 1:10, pve = (eigenval/sum(eigenval)*100))
pve$PC <- as.factor(pve$PC)

# Plot PVE (ancestry) ------------------------------------------------------------------

PC_Palette <-
  heat_hcl(
    10,
    h = c(30, -160),
    c = c(80, NA, 45),
    l = c(38, 79),
    power = c(0.85, 1.0),
    fixup = TRUE,
    gamma = NULL,
    alpha = 1
  )

ggplot(pve, aes(PC, V1, fill = PC)) + geom_bar(stat = "identity") + ylab("Proportion of Variance (%)\n") + scale_fill_manual(values = PC_Palette) + theme(legend.position = "none")

ggplot(pve, aes(PC, V1, colour = PC)) +geom_line(group="PC", colour ="#454b51")+ geom_point(size =2.1)  + ylab("\nProportion of Variance (%)\n") + scale_colour_manual(values = PC_Palette) + theme(legend.position = "none") 

ggsave("pve_ancestry.png")

# PCA - Just study data  -------------------------------------------------------

# Load data
pca <- fread("enteric.LD.eigenvec")
eigenval <- fread("enteric.LD.eigenval")

# Basic PCA
ggplot(data = pca, aes(PC1, PC2)) + geom_point()

# PVE
pve <- data.frame(PC = 1:10, pve = (eigenval/ sum(eigenval) * 100))
pve$PC <- as.factor(pve$PC)
cumsum(pve$V1)

# Plot proportion of variance
ggplot(pve, aes(PC, V1, colour = PC)) +
  geom_line(group = "PC", colour = "#454b51") +
  geom_point(size = 2.1)  +
  ylab("\nProportion of Variance (%)\n") +
  scale_colour_manual(values = PC_Palette) +
  theme(legend.position = "none")
ggsave("pve_study.png")
# Use these PCs in GLM model

# PCA by outcome  -------------------------------------------------------
covar <- fread("~/GWAS_22/gwas_final/meta/covar_enteric.txt")
pca_covar <- left_join(pca, covar, by = "IID") # matched 250, 301
pca_covar$Diagnosed <- as.factor(pca_covar$Diagnosed)
pca_covar <- pca_covar[,c(1:18,26:27)]

write.table(pca_covar,
            file = "~/GWAS_22/gwas_final/meta/enteric_pcacovar.txt",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")

# Replot PCAS
pca_covar %>%
  filter(Diagnosed != "NA") %>%
  ggplot(aes(PC1, PC2, colour = Diagnosed)) +
  geom_point(alpha = 0.9) +
  scale_colour_manual(values = c("seagreen3", "sandybrown"),
                      labels = c("nTD", "TD")
                      )

pca_covar %>%
  filter(Diagnosed != "NA") %>%
  ggplot(aes(PC3, PC4, colour = Diagnosed)) +
  geom_point(alpha = 0.9) +
  scale_colour_manual(values = c("seagreen3", "sandybrown"),
                      labels = c("nTD", "TD")
                      )

# PCA by sex -------------------------------------------------------------------

pca_covar$sex <- as.factor(pca_covar$sex)
# M= 1, F = 2

a <-
  pca_covar %>%
  filter(sex != "NA") %>%
  ggplot(aes(PC1, PC2, colour = sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()
a

b <- 
  pca_covar %>%
  filter(sex != "NA") %>%
  ggplot(aes(PC2, PC3, colour = sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

c <- 
  pca_covar %>%
  filter(sex != "NA") %>%
  ggplot(aes(PC3, PC4, colour = sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

d <- 
  pca_covar %>%
  filter(sex != "NA") %>%
  ggplot(aes(PC4, PC5, colour = sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

e <- 
  pca_covar %>%
  filter(sex != "NA") %>%
  ggplot(aes(PC5, PC6, colour = sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

f <- 
  pca_covar %>%
  filter(sex != "NA") %>%
  ggplot(aes(PC6, PC7, colour = sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

g <- 
  pca_covar %>%
  filter(sex != "NA") %>%
  ggplot(aes(PC7, PC8, colour = sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

h <- 
  pca_covar %>%
  filter(sex != "NA") %>%
  ggplot(aes(PC8, PC9, colour = sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

i <- pca_covar %>%
  filter(sex != "NA") %>%
  ggplot(aes(PC9, PC10, colour = sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

(a|d) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", title = "PCA of imputed gwas data") 
ggsave("PCenteric-12-45.png")


pdf("all-enteric-PCAs.pdf")
a
b
c
d
e
f
g
h
i
dev.off()