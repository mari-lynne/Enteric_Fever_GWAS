# Cohorts
setwd("~/GWAS_22/gwas_final/eQTL/")

study <- c("/vast_tyg")
time_point <- c("/d7")
time <- c("d7")
var <- c("_TYG_")
pathway <- c("~/GWAS_22/gwas_final/InnateDB_genes.csv")

id_data <- c("~/GWAS_22/gwas_final/meta/geno_ids.csv")
plink <- c("~/GWAS_22/gwas_final/merge/typhoid/QC/typhoid2.IBD.fam")
assoc_data <- c("~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus2/typhoid_tophits.txt")
covar_data <- c("~/GWAS_22/gwas_final/meta/typhoid_pcacovar.txt")
out_dir <- paste0(getwd(), study, "/Nov/combat/redo")
'%!in%' <- function(x,y)!('%in%'(x,y))
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
data <- readRDS(file = "combat_final_d7.RDS")
samples <- data$samples
#samples$time_point <- as.factor(samples$time_point)levels(samples$time_point)

samples <- samples %>%
  mutate(time_point = ifelse(time_point3 == "TD", "TD", time_point)) %>%
  mutate(study = ifelse(batch2 == "1", "TyGER", "VAST"))

novax <- filter(samples, chall_vax == "T")
vax <- filter(samples, chall_vax != "T")

novax <- novax %>%
  mutate(time_point = ifelse(time_point %in% c("TD", "D7"), "D7-TD", time_point))

# all data summary
# use word table save as an image!! png

# RNA-seq deg cohorts no vax ---------------------------------------------------
novax %>% filter(time_point != "D0.12h") %>% tabyl(time_point, study, diagnosis) %>% adorn_totals("col") %>%
  kable(booktabs = TRUE, format = "latex")

# VAST cohort
load(file = "~/RNA/diff_expr_results/VAST_RNAseq.RData")

vast <- droplevels(data$samples)
table(vast$study_arm, vast$time_point3)
vast <- filter(vast, diagnosis %in% c("TD", "nTD"))
vast$time_point3 <- as.character(vast$time_point3)
vast <- vast %>%
  mutate(time_point_ = ifelse(time_point3 %in% c("V0"), "Baseline", time_point3))

vast %>% filter(time_point_ %in% c("Baseline", "V1", "V7")) %>% tabyl(time_point_, study_arm) %>% adorn_totals("col")%>% kable(booktabs = TRUE, format = "latex")
