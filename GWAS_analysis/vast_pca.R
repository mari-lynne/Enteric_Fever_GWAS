# Load data
setwd("~/GWAS_22/gwas_final/merge/typhoid/vast")

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

pca <- fread("vast.LD.eigenvec")
eigenval <- fread("vast.LD.eigenval")
pca_covar <- fread("covar_vast.txt")

# Basic PCA
ggplot(data = pca, aes(PC1, PC2)) + geom_point()

# PVE
pve <- data.frame(PC = 1:10, pve = (eigenval/ sum(eigenval) * 100))
pve$PC <- as.factor(pve$PC)
cumsum(pve$V1)

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
ggsave("pve_study.png")
# Use these PCs in GLM model

# PCA by outcome  -------------------------------------------------------
pca_covar$Diagnosed <- as.factor(pca_covar$Diagnosed)

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
  ggplot(aes(PC4, PC5, colour = Diagnosed)) +
  geom_point(alpha = 0.9) +
  scale_colour_manual(values = c("seagreen3", "sandybrown"),
                      labels = c("nTD", "TD")
  )

# PCA by sex -------------------------------------------------------------------

#TODO - Turn this into a function/feed in data from bash as before
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
ggsave("PCvast-12-45.png")


pdf("vast-PCAs.pdf")
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