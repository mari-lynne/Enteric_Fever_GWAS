#  new vast tyger -----------------------------------------------------------
  load(file = "VAST_Tyger_oct.RData")
  
#1) Update metadata vast to include extra vars from pheno (renamed to match tyger)
# 2) Add pheno/meta info into vast samples  (remove lab and p_id from meta to avoid dups)
# 3) Filter vast samples for matching colnames in tyger samples, and vice versa
# 4 ) bind rows

#1) 
  # Intersecting metacolumns
 vast_pheno <- clean_names(read.csv(file = "~/GWAS_22/gwas_final/meta/VAST_meta.csv"))
  vast_pheno <- vast_pheno %>% dplyr::rename(blood_quant_cfu_ml = bacterial_cfu_m_l,
                                             time_to_diagnosis_abx = time_to_diag_days,
                                             any_temperature_38 = temp_at_diagnosis_for_vast_this_is_fever_over_38)
  
  check <- intersect(names(tyger$metadata), names(vast_pheno))
  check
  vast_pheno <-  vast_pheno[ ,colnames(vast_pheno) %in% check]
  vast$meta_data <- left_join(vast$meta_data, vast_pheno, by = c("lab_id", "participant_id"))
  # refilter for matching meta cols
  check <- intersect(names(tyger$metadata), names(vast$meta_data))
  check
  vast$meta_data <-   vast$meta_data[ ,colnames(vast$meta_data) %in% check]
  
  # Filter tyger metadata 
  tyger$metadata <-  tyger$metadata[ ,colnames(tyger$metadata) %in% check]
  
  # 2) add to respective sample df (remove extra ids)
  extra <- intersect(names(vast$meta_data), names(vast$samples))
  vast$meta_data <- select(vast$meta_data, -c(extra))
  vast$samples <- bind_cols(vast$samples, vast$meta_data)
  names(vast$samples)
  
  # 2a) Filter tyger metatdata, add to samples
  extra <- intersect(names(tyger$metadata), names(tyger$samples))
  extra
  tyger$metadata <- select(tyger$metadata, -c(extra))
  tyger$samples <- bind_cols(tyger$samples, tyger$metadata)
  names(tyger$metadata)
  
  check <- intersect(names(tyger$samples), names(vast$samples))
  check
  tyger$samples <- tyger$samples[,colnames(tyger$samples) %in% check]
  dim(tyger$samples)
  dim(vast$samples)
  
  # New DEG list
  counts <- cbind(tyger$counts,vast$counts)
  samples <- rbind(tyger$samples,vast$samples)
  genes <- vast$genes
  dge <- DGEList(counts = counts, samples = samples, genes = genes)  
  
  # Filter time points for batch correction ------------------------------------
  # Filter contaminated
  # should filter for contaminated samples if I'm going to just use the uncontaminated 
  samples <- filter(samples, contaminated_sample == "No")
  
  # remove baseline samples in original pool who have lab ids in check
  #check <- filter(samples, sequence_pool == "resequence" & time == "0")
  #check <- filter(samples, (time == "0" & sequence_pool == "original" & lab_id %in% check$lab_id))
  #samples <- filter(samples, rownames(samples) %!in% rownames(check))
  
  # Clean diagnosis var
  samples <- filter(samples, diagnosis != "UNKNOWN - Not Challenged")
  
  # Clean diagnosis var
  table(samples$diagnosis)
  samples <- samples %>%
    mutate(diagnosis = ifelse(diagnosis == "1", "TD", diagnosis),
           diagnosis = ifelse(diagnosis == "0", "nTD", diagnosis))
  
  # Update counts 
  counts <- counts[ , colnames(counts) %in%  rownames(samples)]
  
  # Make new batch sequence pool var
  samples$batch_pool <- str_c(samples$batch, samples$sequence_pool, sep = "_")
  table(samples$batch_pool)
  
  # Harmonise time points -------------------------------------------------------
  
  # recode time_point3
  table(samples$time_point)
  # D0.12h <- D0+12
  # D0 <- D00
  samples <- samples %>%
    mutate(time_point3 = ifelse(time_point3 == "D0+12h", "D0.12h", time_point3)) %>%
    mutate(time_point3 = ifelse(time_point3 == "D00", "D0", time_point3))
  
  # Make day of challenge V28 for vaccine samples day 0
  samples <- samples %>%
    mutate(time_point = ifelse(time_point3 == "D0" & (study_arm %in% c("ViTCV", "ViPS")), "V28", time_point3))
  
  # Make new baseline time point, V0 for vaccinees and D0 for challenge
  samples <- samples %>%
    mutate(time_point = ifelse(time_point %in% c("V0","D0"), "Baseline", time_point)) 
  
  # Day 7 time point, around D7 for TD (time point3 still contains TD)
  samples <- samples %>% 
    mutate(time_point = ifelse(time %in% c("6","7","8") & (time_point3 == "TD"), "D7", time_point))
  
  table(samples$time_point3, samples$diagnosis)
  # filter vaccine timepoints and D01-14 which are study specific
  samples <- filter(samples, time_point %in% c("Baseline", "D0.12h") | time_point3 == "TD")
  
  # Just D7
  # samples <- filter(samples, time_point == "D7")
  tab <- table(samples$time_point3, samples$batch_pool, samples$diagnosis)
  addmargins(tab)
  
  samples <- samples %>%
    mutate(time_point3 = ifelse(time_point3 %in% c("V0","D0"), "Baseline", time_point3)) 
  
  # Merge TN and control samples
  samples <- samples %>% 
    mutate(chall_vax = ifelse(str_detect(study_arm, "TN|CTRL"), "T", study_arm))
  # vast were vaccinated, therefore should add in MenC as option for vax
  
  counts <- counts[ , colnames(counts) %in%  rownames(samples)]
  
  covar <- model.matrix(~chall_vax, data = samples)
  
  combat <- sva::ComBat_seq(counts=counts,
                            batch=samples$batch_pool,
                            group=samples$time_point3,
                            covar_mod = covar)
  
  
  
  # the other option is to just correct between batches 1,2,3. Then add in
  # batch pool again as a covar in any dge/eqtl analysis
  # or within study batch correction, then interstudy
  # or just exclude contaminated baseline samples
  # Tyger batch correction -----------------------------------------------------
  
  table(tyger$samples$time_point3, tyger$samples$batch)
  table(tyger$samples$time_point3, tyger$samples$contaminated_sample)
  tyger$samples <- filter(tyger$samples, contaminated_sample == "No")
  table(tyger$samples$time_point3, tyger$samples$batch)
  
  table(vast$samples$time_point3, vast$samples$batch)
  
  # theres no tyger d7 that could have also made things weird
  tyger$samples <- filter(tyger$samples, contaminated_sample == "No")
  tyger$samples <- filter(tyger$samples, time_point3 %in% c("D0.12h", "D00", "D01", "TD"))
  # have to correct for same time points in vast #
  tyger$counts <- tyger$counts[ , colnames(tyger$counts) %in%  rownames(tyger$samples)]
  
  tyger$samples <- droplevels(tyger$samples)
  covar <- model.matrix(~study_arm, data = tyger$samples)
  table(tyger$samples$study_arm, tyger$samples$batch)
  
  tyger_c <- sva::ComBat_seq(counts=tyger$counts,
                            batch=tyger$samples$batch,
                            group=tyger$samples$time_point3,
                            covar_mod = covar)
  
  
  tyger_bc <- DGEList(counts = tyger_c, samples = tyger$samples, genes = tyger$genes)
  tyger_bc$raw  <- tyger$counts
  
  
  # vast-bc ----------------------------------------------------------------
  table(vast$samples$time_point3, vast$samples$study_arm)
  vast$samples <- droplevels(vast$samples)
  table(vast$samples$time_point3, vast$samples$study_arm)
  
  # Make day of challenge V28 for vaccine vast$samples day 0
  samples <- vast$samples
  samples$time_point3 <- as.character(samples$time_point3)
  samples <- samples %>%
    mutate(time_point = ifelse(time_point3 == "D0" & (study_arm %in% c("ViTCV", "ViPS")), "V28", time_point3))
  # Make new baseline time point, V0 for vaccinees and D0 for challenge
  samples <- samples %>%
    mutate(time_point = ifelse(time_point %in% c("V0","D0"), "Baseline", time_point))
  
  table(samples$time_point, samples$study_arm)
  
  samples <- samples %>%
    mutate(time_point = ifelse(time_point == "D0+12h", "D0.12h", time_point))
  
  # Filter vast
  table(samples$time_point, samples$sequence_pool)
  table(samples$time_point, samples$diagnosis)
  samples <- filter(samples, diagnosis != "UNKNOWN - Not Challenged")
  samples <- filter(samples, time_point %in% c("D0.12h", "Baseline", "D7", "TD"))
  samples <- droplevels(samples)
  
counts <- vast$counts
  counts <- counts[ , colnames(counts) %in%  rownames(samples)]
  
  
  covar <- model.matrix(~study_arm, data = samples)
  table(samples$time_point, samples$sequence_pool)
  
  vast_c <- sva::ComBat_seq(counts=counts,
                             batch=samples$sequence_pool,
                             group=samples$time_point,
                             covar_mod = covar)
  

  
  # Save data -----------------------------------------------------------
  
  vast_bc <- DGEList(counts = vast_c, samples = samples, genes = vast$genes)
  vast_bc$raw  <- counts
   
  # recode tyger vars to match
  tyger_bc$samples <- tyger_bc$samples %>%
    mutate(time_point = ifelse(time_point3 == "D00", "Baseline", time_point3))
  
  # New DEG list
  counts <- cbind(tyger_bc$counts,vast_bc$counts)
  samples <- rbind(tyger_bc$samples,vast_bc$samples)
  genes <- vast_bc$genes
  dge <- DGEList(counts = counts, samples = samples, genes = genes)
  
  idx <- (dge$samples$time_point != "D7")
  dge <- dge[,idx]
  dim(dge)
  
  # remove d7 from merged
  dge$samples <- dge$samples %>% 
    mutate(chall_vax = ifelse(str_detect(study_arm, "TN|CTRL"), "T", study_arm))
  
  save.image(file = "vt_bc.nov8.RData")
 
  # final correction -----------------------------------------------------------
  table(dge$samples$batch, dge$samples$time_point)
  # filter contaminated samples, recode to b1
  dge$samples <- filter(dge$samples, contaminated_sample == "No")

  dge$samples <- dge$samples %>% 
    mutate(batch = ifelse(batch == 2, 1, batch))
  
  idx <- (dge$samples$time_point != "D01")
  dge <- dge[,idx]
  
  table(dge$samples$time_point, dge$samples$chall_vax)
  
  
  covar <- model.matrix(~chall_vax, data = dge$samples)
  
  dge_c <- sva::ComBat_seq(counts=dge$counts,
                            batch=dge$samples$batch,
                            group=dge$samples$time_point,
                            covar_mod = covar) 
   
  # TPM
  
  combat <- DGEList(counts = dge_c, samples = dge$samples, genes = vast$genes)
  combat$unstransformed  <- dge$raw
  combat$combat_counts <- dge_c
  
  
  tpm3 <- function(counts,len) {
    x <- counts/len
    return(t(t(x)*1e6/colSums(x)))
  }
  
  
  tail(vast$genes)
  tail(tyger$genes)
  
  #tpm1 <- tpm3(sub$counts, tyger$genes$transcript_length)
  #tpm3 <- tpm3(sub$counts, vast$genes$transcript_length)
  #tpm_bound <- cbind(tpm1, tpm3)
  
  
  tpm <- tpm3(combat$counts, combat$genes$transcript_length)

  combat$counts <- tpm
  combat$samples <- combat$samples %>%
    mutate(diagnosis = ifelse(diagnosis == "1", "TD", diagnosis),
           diagnosis = ifelse(diagnosis == "0", "nTD", diagnosis))
  
  saveRDS(combat, file = "combat_final.RDS")
  
  
  # Fever assoiations (hopefully with higher plat expression)
  qtls <- c("F2R","MPL","ITGB3","MMRN1","GP9","GP1BA", "CD9", "THBS1", "ANGPT1","PTGS1")
  
  idx <- (vast_bc$samples$time %in% c(6:9))
  sub <- vast_bc[,idx]
  dim(sub)
  
  # get genes and rsids
  idx <- sub$genes$gene_name %in% qtls
  sub <- sub[idx,]
  dim(sub$genes)
  
  exprs <- as.data.frame(t(sub$counts))
  names(exprs) <- sub$genes$gene_name
  pheno <- as.data.frame(sub$samples)
  pheno_exprs <- bind_cols(pheno, exprs)
  
  data <- pheno_exprs %>% filter(!is.na(diagnosis))
  
  
  # td ntd comparisons
  compare <- list(c("nTD", "TD"))
  
#  v_d7 <-#
    
    data %>% filter(batch == "3") %>%
  ggplot(aes(x = diagnosis, y = F2R, fill = diagnosis)) +
    geom_violin()  +
    geom_jitter(width = 0.1,
                size = 1.7,
                alpha = 0.5) +
    scale_colour_manual(values = c("#061c23")) + theme_bw() +
  labs(y = "Counts", title = "Day 7 - F2R expression") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 11.5),
      axis.text.y = element_text(size = 11.5)
    ) + scale_fill_viridis(discrete = TRUE) + stat_compare_means(
      comparisons = compare,
      method = "wilcox.test",
      label = "p.format",
      size = 5, 
      p.adjust.method = "none")
  
  
  #  plat_cors ----------------------------------------------------------------
    # corrupt file vast_pheno <- clean_names(read.csv(file = "~/GWAS_22/gwas_final/meta/VAST_meta.csv"))
    
  vast_pheno <- vast_pheno %>% dplyr::rename(blood_quant_cfu_ml = bacterial_cfu_m_l,
                                               time_to_diagnosis_abx = time_to_diag_days,
                                               any_temperature_38 = temp_at_diagnosis_for_vast_this_is_fever_over_38)
    
    
    plate <- left_join(vast_pheno, pheno_exprs, by = "lab_id")
    
    
    plate %>% filter(batch == "3") %>%
      ggplot(aes(x = diagnosis, y = platelets_x10_9_l, fill = diagnosis)) +
      geom_violin()  +
      geom_jitter(width = 0.1,
                  size = 1.7,
                  alpha = 0.5) +
      scale_colour_manual(values = c("#061c23")) + theme_bw() +
      labs(y = "Platelets x109/l", title = "Platelet Count", x = "Diagnosis") +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 11.5),
        axis.text.y = element_text(size = 11.5)
      ) + scale_fill_viridis(discrete = TRUE) + stat_compare_means(
        comparisons = compare,
        method = "wilcox.test",
        label = "p.format",
        size = 5, 
        p.adjust.method = "none")
    
    
    #(f2r_d0 |f2r_h12) + plot_annotation(tag_levels = "a", title = "F2R gene expression")
    
    # Plate cors redo ------------------------------------------------------------------
    
  # Clinical measurements
    library(readODS)
    plate <- read.ods(file = "~/GWAS_22/gwas_final/meta/vast_tyger_blood.ods", sheet = 1)
    names(plate) <- plate[1,]
    plate <- clean_names(plate)
    plate <- plate[-1,]
    cols <- plate[, c(6:ncol(plate))] %>% colnames
    plate[, cols] <- lapply(plate[, cols], as.numeric)
  
  # Pheno gwas and genotype info
    setwd("~/GWAS_22/gwas_final/eQTL/vast_tyg")
    #keep <- filter(pca_covar, study == "VAST" | study == "TYGER") %>% select(`#FID`, "IID")
    #write.table(keep, file = "vast_tyg_all_keep.txt", row.names = F, quote = F)
    
    gwas <- fread("~/GWAS_22/gwas_final/merge/typhoid/assoc/nexus2/typhoid_tophits.txt")
    geno <- fread("~/GWAS_22/gwas_final/eQTL/vast_tyg/chr5_8.traw") # chr5_loci.traw
  # Prep geno 
    id_data <- clean_names(fread("~/GWAS_22/gwas_final/meta/geno_ids.csv"))
    id_data$geno_id <- str_sub(id_data$cat_iid, start= -4)
    
    snp <- geno[,-c(1,3:6)]
    snp <-  t(snp)
    colnames(snp) <- snp[1, ]
    snp <- snp[-1, ]
    snp <- as.data.frame(snp)
    snp <- tibble::rownames_to_column(snp, var = "geno_id")
    colnames(snp) <- str_replace_all(colnames(snp), "\\:", "_")
    snp$geno_id <- str_sub(snp$geno_id, start= -4)
    snp <- left_join(snp, id_data, by = "geno_id")
  # Add pheno
    plate$lab_id <- as.numeric(plate$id)
    snp_pheno <- left_join(snp, plate, by = "lab_id")
    colnames(snp_pheno) <- str_replace_all(colnames(snp_pheno), "chr", "Chr")
    cols <- snp_pheno %>% dplyr::select(starts_with("Chr")) %>% colnames()
    snp_pheno[,cols] <- lapply(snp_pheno[,cols], factor)
  
    # chr5:77685638:A/G -rs snp # Chr5_77682126_C_T -eqtl snp # rs62364881pub
    #  8  133074620  rs2739156  T  A,C also
    
    # SNPs of interest to plot # , "Chr_8_133074620_T_A"
    snp_names <-
      c("Chr5_77682126_C_T", "Chr5_77685638_A_G", "Chr5_77775879_G_T", "Chr5_77795103_G_A")
    
    y = snp_pheno$platelets_x10_9_l
    for (i in snp_names) {
      plt <-
        ggplot(snp_pheno, aes_string(x = i, y = y)) +
        geom_boxplot(outlier.shape = NA, aes_string(fill = i)) +
        scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
        geom_jitter(width = 0.1,
                    size = 1.7,
                    alpha = 0.5) +
        scale_colour_manual(values = c("#061c23"))  +
        labs(y = "Platelet count", title = "TD") +
        theme_bw()  +
        theme(
          legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.title = element_text(size = 11.5),
          axis.text.y = element_text(size = 11.5)
        )
      print(plt)
      Sys.sleep(0.5)
    }
    
    
    
    #### eqtl OR high gene -----------
    compare <- list(c("0","1"), c("0", "2"))
    snp_pheno %>% filter(!is.na(Chr5_77682126_C_T), vaccine %!in% c("Vi-PS", "Vi-TCV")) %>%
    ggplot(aes(x = Chr5_77682126_C_T, y = platelets_x10_9_l, fill = Chr5_77682126_C_T)) +
      geom_boxplot(outlier.shape = NA) +
      scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
      geom_jitter(width = 0.1,
                  size = 1.7,
                  alpha = 0.5) +
      scale_colour_manual(values = c("#061c23"))  +
      labs(y = "Platelet count", title = "Baseline") +
      theme_bw()  +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 11.5),
        axis.text.y = element_text(size = 11.5)
      )
    
    
    ### rs gene ----
    snp_pheno %>% filter(!is.na(Chr5_77685638_A_G), study == "TYGER") %>%
      ggplot(aes(x = Chr5_77685638_A_G, y = platelets_x10_9_l, fill = Chr5_77685638_A_G)) +
      geom_boxplot(outlier.shape = NA) +
      scale_fill_manual(values = c("#A92D2F", "#CE7726", "#9B9759")) +
      geom_jitter(width = 0.1,
                  size = 1.7,
                  alpha = 0.5) +
      scale_colour_manual(values = c("#061c23")) +
      labs(y = "Platelet count\n", title = "Baseline - Tyger") +
      theme_bw()  +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 11.5),
        axis.text.y = element_text(size = 11.5)
      )
    
    
    # only trend seen in tyger - vast have vaccinated participants at baseline in control group though
    
    ### Diagnosis --------
   
     plate %>% filter(!is.na(outcome)) %>%
      ggplot(aes(x = outcome, y = platelets_x10_9_l, fill = outcome)) +
      geom_violin()  +
      geom_jitter(width = 0.1,
                  size = 1.7,
                  alpha = 0.5) +
      scale_colour_manual(values = c("#061c23")) + theme_bw() +
      labs(y = "Platelets x109/l", title = "Platelet Count", x = "Diagnosis") +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 11.5),
        axis.text.y = element_text(size = 11.5)
      ) + scale_fill_viridis(discrete = TRUE) + stat_compare_means(
        comparisons = compare,
        method = "wilcox.test",
        label = "p.format",
        size = 5, 
        p.adjust.method = "none")
    
    
   ### Day7 -----------------
    
    # Gene expression plots ----------------
    
    
    
   
   
    
    
   
   # Set up labels
    
    # Label P-val and OR
    gwas$SNP <- str_replace_all(gwas$SNP , "chr", "Chr")
    gwas$SNP <- str_replace_all(gwas$SNP, ":", "_")
    # Get summary data
    lab <- gwas[SNP == "Chr5_77682126_C_T", OR]
    lab <- signif(lab, digits = 2)
    lab_p <- gwas[SNP == "Chr5_77682126_C_T", P]
    lab_p <- signif(lab_p, digits = 2)
    data <- snp_pheno %>% filter(!is.na(Chr5_77682126_C_T))
  
  # chall_vax == "T"
    
  
  thbs1_d0 <-
    data %>% filter(chall_vax == "T") %>%
    ggplot(aes(x = diagnosis, y = THBS1, fill = diagnosis)) +
    geom_violin()  +
    geom_jitter(width = 0.1,
                size = 1.7,
                alpha = 0.5) +
    scale_colour_manual(values = c("#061c23")) + theme_bw() +
    labs(y = "Gene Expression (TPM)", title = "THBS1") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 11.5),
      axis.text.y = element_text(size = 11.5)
    ) + scale_fill_viridis(discrete = TRUE) 
  
  thbs1_d0
  
  ang_d0 <- data %>% filter(chall_vax == "T") %>%
    ggplot(aes(x = diagnosis, y = ANGPT1, fill = diagnosis)) +
    geom_violin()  +
    geom_jitter(width = 0.1,
                size = 1.7,
                alpha = 0.5) +
    scale_colour_manual(values = c("#061c23")) + theme_bw() +
    labs(y = "Gene Expression (TPM)", title = "ANGPT1") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 11.5),
      axis.text.y = element_text(size = 11.5)
    ) + scale_fill_viridis(discrete = TRUE) 
  
  ang_d0
  
vt_12 <- (f2r_d0|thbs1_d0|ang_d0) + plot_annotation(title = "12h", tag_levels = "a")
  

# # +
#   # stat_compare_means(
#       comparisons = compare,
#       method = "wilcox.test",
#       label = "p.format",
#       size = 5, 
#       p.adjust.method = "none", label.y = 55
#     ) +
#   scale_fill_manual(values = c("seagreen3", "sandybrown"))
  
  
  compare <- list(c("nTD", "TD"))
  data %>% filter(chall_vax == "T") %>%
    ggplot(aes(x = diagnosis, y = THBS1, fill = diagnosis)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = c("#A92D2F", "#CE7726")) +
    geom_jitter(width = 0.1,
                size = 1.7,
                alpha = 0.5) +
    scale_colour_manual(values = c("#061c23"))  +
    stat_compare_means(
      comparisons = compare,
      method = "wilcox.test",
      label = "p.format",
      size = 5, 
      p.adjust.method = "none"
    ) +
    labs(y = "Gene Expression (TPM)", title = "THBS1 - 12h") +
    theme_bw()  +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 11.5),
      axis.text.y = element_text(size = 11.5)
    ) 
  
  
  # symptom cors
  data %>% filter(batch == 3) %>%
    ggplot(aes(x = F2R, y = time_to_diagnosis_abx, colour = batch))+
    geom_point() + geom_smooth(method = lm) + stat_cor(method = "spearman")
  
  data %>%
    ggplot(aes(x = F2R, y = duration_bacteraemia_hrs, colour = batch))+
    geom_point() + geom_smooth(method = lm)  +
    xscale("log10") + yscale("log10", .format = TRUE) + stat_cor(method = "spearman")
  

  