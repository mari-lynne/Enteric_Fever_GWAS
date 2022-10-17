# Date: 2022/06/02
# Purpose: Plot user data projected on 1KG PCs, greying out 1KG individuals

# Updates 12/10/22
# Reran PCA of study data after pruning
# Have separated into just typhoid cohort ( n=253)

# Load data --------------------------------------------------------------------

setwd("/home/mari/GWAS_22/gwas_final/merge/QC/pca")

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
  "all_hg38.LD.eigenvec",
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

c <- 
  ggplot(merge, aes(PC3, PC4, colour = SuperPop)) +
  geom_point() + scale_colour_manual(values = KG_Palette_Super) +
  theme_bw()

c

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

i <- ggplot(merge, aes(PC9, PC10, colour = Population)) +
  geom_point() + scale_colour_manual(values = KG_Palette) +
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

eigenval <- fread("all_hg38.LD.eigenval")
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
pca <- fread("typhoid.LD.eigenvec")
eigenval <- fread("typhoid.LD.eigenval")

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
pheno <- fread("~/GWAS_22/new_gwas/meta/final_metadata/all_pheno_num.txt")
covar <- fread("~/GWAS_22/new_gwas/meta/final_metadata/all_covar_num.txt")

# covar checks
covar_check <- covar %>% filter(Challenge != 3, Re_Challenge != 5)

# Should have 259 participants - I think use these to filter out enteric_QC3 

pca_covar <- left_join(pca, pheno, by = "IID") # matched 242
pca_covar$Diagnosed <- as.factor(pca_covar$Diagnosed)
covar <- left_join(pca_covar, covar, by = "IID")

# Make new covar file
# Make new categorical factors
covar$chall_vax <- str_c(rep("cv", length(covar)), covar$Challenge, covar$Vaccine, sep = "_")
covar$chall_vax <- as.factor(covar$chall_vax)

covar <- covar[,c(1:20,30)]

write.table(covar,
            file = "/home/mari/GWAS_22/new_gwas/just_typhoid/meta/typhoid_covar.txt",
            quote = F,
            row.names = F,
            sep = "\t")

pheno <- covar[,c(1,2,13)]
write.table(pheno,
            file = "/home/mari/GWAS_22/new_gwas/just_typhoid/meta/typhoid_pheno.txt",
            quote = F,
            row.names = F,
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
  ggplot(aes(PC2, PC3, colour = Sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

c <- 
  covar %>%
  filter(Sex != "NA") %>%
  ggplot(aes(PC3, PC4, colour = Sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

d <- 
  covar %>%
  filter(Sex != "NA") %>%
  ggplot(aes(PC4, PC5, colour = Sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

e <- 
  covar %>%
  filter(Sex != "NA") %>%
  ggplot(aes(PC5, PC6, colour = Sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

f <- 
  covar %>%
  filter(Sex != "NA") %>%
  ggplot(aes(PC6, PC7, colour = Sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

g <- 
  covar %>%
  filter(Sex != "NA") %>%
  ggplot(aes(PC7, PC8, colour = Sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

h <- 
  covar %>%
  filter(Sex != "NA") %>%
  ggplot(aes(PC8, PC9, colour = Sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

i <- covar %>%
  filter(Sex != "NA") %>%
  ggplot(aes(PC9, PC10, colour = Sex)) +
  geom_point(alpha = 0.95, size = 2) +
  scale_colour_manual(values = c("#a65ca6", "#fcd325"),
                      labels = c("Male", "Female")) +
  theme_bw()

(a|d) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", title = "PCA of imputed gwas data") 
ggsave("PCtyphoid-12-45.png")


pdf("all-typhoid-PCAs.pdf")
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




# Loop attempts ####
pc_names <- colnames(covar[,c(3:12)])
pc_plot <- covar[,c(3:12)]


for(i in (pc_names)) {
  print(colnames(covar[1, (i+1)]))
}


for(i in 1:ncol(pc_plot)) {
  for (j in 2:ncol(pc_plot)) {
    if (i == (j - 1))
      print(colnames(pc_plot[1, ..j]))
  }
}


for i1 in PC1, J = 2-10
for i2 in PC2, J = 2-10

for(i in 9:ncol(pc_plot)) {
  for (j in 10:ncol(pc_plot)) {
    if ((i == (j - 1)) && (j == (i +1)))
      plot <- covar %>%
        filter(Sex != "NA") %>%
        ggplot(aes_string(
          x = colnames(pc_plot[, ..i]),
          y = (colnames(pc_plot[, ..j])),
          colour = "Sex"
        )) +
        geom_point(alpha = 0.95,
                   size = 2) +
        scale_colour_manual(values =
                              c("#a65ca6", "#fcd325"),
                            labels =
                              c("Male", "Female")) +
        theme_bw()
    print(plot)
    Sys.sleep(0.5)
  }
}


area = list()  # because the actual function doesn't work
for (i in 1:ncol(pc_plot)) {
  for (j in 1:ncol(pc_plot)) {
    if (i == (j - 1)) {
      M[i, i] = 0
      next
    }
    selection = df[, c(i, j)]
    #area=integrate(f2, 1, 200, subdivisions = 500)
    area$value = mean(colSums(selection)) # something random to check
    M[i, j] = area$value
    M[j, i] = area$value
  }
}


for (i in 1:ncol(pc_plot)) {
  for (j in 2:ncol(pc_plot)) {
    if (i == (j - 1))
      covar %>%
      select(Sex, i, j) %>%
      print(colnames(i))
  }
}



ggplot(aes_string(
  x = covar[1, ..i],
  y = covar[1, ..j],
  colour = "Sex"
))  +
  geom_point(alpha = 0.95)
#print(pc_plot[1, c(i, j)])

pc_plot[1, c(1, 2)]

for(i in 1:ncol(pc_plot)) {
  for (j in 1:ncol(pc_plot)) {
    if ((i == (j - 1)) && (j == (i +1)))
      plot <- covar %>% select(i,j) %>%
        filter(Sex != "NA") %>%
        ggplot(aes_string(
          x = colnames(pc_plot[1, ..i]),
          y = colnames(pc_plot[1, ..j]),
          colour = "Sex"
        )) +
        geom_point(alpha = 0.95,
                   size = 2) +
        scale_colour_manual(values =
                              c("#a65ca6", "#fcd325"),
                            labels =
                              c("Male", "Female")) +
        theme_bw()
    print(plot)
    Sys.sleep(0.5)
  }
}








for(i in 1:ncol(pc_plot)) {
  print(colnames(pc_plot[1, ..i]))
}


# This goes through all of j for i, then pc2 v pc3,4,5
for(i in 1:ncol(pc_plot)) {
  for(j in 2:ncol(pc_plot)) {
  print(colnames(pc_plot[, ..i]))
    print(colnames(pc_plot[, ..j]))
  }
}

while()


# could make j a sequence to plot starting from pc2

while i == j-1


for(i in (pc_names)) {
  plot <- covar %>%
    filter(Sex != "NA") %>%
    ggplot(aes_string(x=i, y=(i + 1), colour = "Sex")) +
    geom_point(alpha = 0.95,
               size = 2) +
    scale_colour_manual(values =
                          c("#a65ca6", "#fcd325"),
                        labels =
                          c("Male", "Female")) +
    theme_bw()
  print(plot)
  Sys.sleep(0.5)
}
