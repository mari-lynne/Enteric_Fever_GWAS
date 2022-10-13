# Date: 2022/06/02
# Purpose: Plot user data projected on 1KG PCs, greying out 1KG individuals

# Load data --------------------------------------------------------------------

setwd("/home/mari/GWAS_22/new_gwas/QC/ancestry")

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
  "1KG.merged.eigenvec",
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

merge <- full_join(PCA, pop, by = "IID")

# Plot ancestry set up ---------------------------------------------------------

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

x <- 
  ggplot(merge,
       aes(PC1, PC2, colour = SuperPop)) +
  geom_point() +
  scale_colour_manual(values = (KG_Palette_Super)) +
  theme_bw()

dev.off()

# Sub population PCA

pdf("KG.subpop_strat_PCA.pdf")

y <- 
  ggplot(merge,
       aes(PC1, PC2, colour = Population)) +
  geom_point() +
  scale_colour_manual(values = KG_Palette) +
  theme_bw()

dev.off()

pdf("KG_ancestry_PCAs.pdf")
(x/c) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", title = "1KG ancestry PCA") 
dev.off()

ggsave("KG-ancestry-PCAs.png")


# Extra PCs

b <-
  ggplot(merge, aes(PC2, PC3, colour = SuperPop)) +
  geom_point() +
  scale_colour_manual(values = KG_Palette_Super)
c <- 
  ggplot(merge, aes(PC3, PC4, colour = SuperPop)) +
  geom_point() + scale_colour_manual(values = KG_Palette_Super) +
  theme_bw() + theme(legend.position = "none")


d <- 
  ggplot(merge, aes(PC4, PC5, colour = SuperPop)) +
  geom_point() + scale_colour_manual(values = KG_Palette_Super) + theme_bw()

e <- 
  ggplot(merge, aes(PC9, PC10, colour = SuperPop)) +
  geom_point() + scale_colour_manual(values = KG_Palette_Super)

e

ggplot(merge, aes(PC8, PC10, colour = Population)) +
  geom_point() + scale_colour_manual(values = KG_Palette)


# Add pc variance --------------------------------------------------------------

eigenval <- fread("1KG.merged.eigenval")
pve <- data.frame(PC = 1:10, pve = (eigenval/sum(eigenval)*100))
pve$PC <- as.factor(pve$PC)

# Plot PVE ---------------------------------------------------------------------
PC_Palette <-
  heat_hcl(
    10,
    h = c(320, 50),
    c. = c(35, 95),
    l = c(15, 90),
    power = c(0.8, 1.2),
    fixup = TRUE,
    gamma = NULL,
    alpha = 1
  )


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

g <- ggplot(pve, aes(PC, V1, fill = PC)) + geom_bar(stat = "identity") + ylab("Proportion of Variance (%)\n") + scale_fill_manual(values = PC_Palette) + theme(legend.position = "none")

h <-ggplot(pve, aes(PC, V1, colour = PC)) +geom_line(group="PC", colour ="#454b51")+ geom_point(size =2.1)  + ylab("\nProportion of Variance (%)\n") + scale_colour_manual(values = PC_Palette) + theme(legend.position = "none")

layout <- c(
  area(1,1, 1,2),
  area(1,3, 1,4),
  area(1,5),
  area(2,1, 2,3)
)
# plot(layout)

# Save plots to pdf ------------------------------------------------------------
pdf("enteric_pop_pca.pdf")
(a | b) / (c |d) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a', title = 'Population structure of enteric CHIM samples')

e + f + guide_area() + h +
  plot_annotation(tag_levels = list(c('e', "f", "g")), title = 'Sub-Population structure of enteric CHIM samples') + plot_layout(design = layout, guides = 'collect')
dev.off()

save.image(file = "enteric_pop_pca.RData")

# PCA no IKG --------------------------------------------------------------------

# Load data
pca <- fread("all_enteric_QC3.cleaned.eigenvec")
eigenval <- fread("all_enteric_QC3.cleaned.eigenval")

# Basic PCA
ggplot(data = pca, aes(PC1, PC2)) + geom_point()

# PVE
pve <- data.frame(PC = 1:10, pve = (eigenval/ sum(eigenval) * 100))
pve$PC <- as.factor(pve$PC)
cumsum(pve$V1)


library(ggfortify)
pca_residues = prcomp(msd[, -c(1:3,18)], scale. = TRUE) # calcs pca
plot(pca_residues$sdev^2/sum(pca_residues$sdev^2), type = "b") # scree plot
autoplot(pca_residues, data = msd, colour = "Timepoint", shape="Vaccine", frame = TRUE, frame.type = 'norm')




# Set up colour palette
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


# Plot proportion of variance
ggplot(pve, aes(PC, V1, colour = PC)) +
  geom_line(group = "PC", colour = "#454b51") +
  geom_point(size = 2.1)  +
  ylab("\nProportion of Variance (%)\n") +
  scale_colour_manual(values = PC_Palette) +
  theme(legend.position = "none")
# Use these PCs in GLM model

# Colour PCA by outcome  -------------------------------------------------------
pca <- fread("all_enteric_QC3.cleaned.eigenvec")
eigenval <- fread("all_enteric_QC3.cleaned.eigenval")
pheno <- fread("~/GWAS_22/new_gwas/meta/final_metadata/all_pheno_num.txt")
covar <- fread("~/GWAS_22/new_gwas/meta/final_metadata/all_covar_num.txt")

pca_covar <- full_join(pca, pheno, by = "IID") # matched 329 rows
pca_covar$Diagnosed <- as.factor(pca_covar$Diagnosed)
covar <- left_join(pca_covar, covar, by = "IID")


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

# Covar PCAs -------------------------------------------------------------------
covar$Sex <- as.factor(covar$Sex)
# M= 1, F = 2

a <-
  covar %>%
  filter(Sex != "NA") %>%
  ggplot(aes(PC1, PC2, colour = Sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

b <- 
  covar %>%
  filter(Sex != "NA") %>%
  ggplot(aes(PC3, PC4, colour = Sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

                                                                                         
(a|b) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", title = "PCA of imputed gwas data") 


ten <-   covar %>%
  filter(Sex != "NA") %>%
  ggplot(aes(PC9, PC10, colour = Sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

ten


c <- ggplot(pve, aes(PC, V1, colour = PC)) +
  geom_line(group = "PC", colour = "#454b51") +
  geom_point(size = 2.1)  +
  ylab("\nProportion of Variance (%)\n") +
  scale_colour_manual(values = PC_Palette) +
  ylim(0,20) +
  theme_bw() +
  theme(legend.position = "none")

c


# PCs proportion of variance only goes down to 9%?
# Either need to include more PCs - double check sib removal
# Include PCs 1-10 to account for structure