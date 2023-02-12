#  new vast tyger -----------------------------------------------------------
setwd("~/RNA/oct")
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
  dim(dge)
  
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
  
  # Tyger-vast study nums checking nov25 -------
  
  str(vast_bc$samples)
  samples <- vast_bc$samples
  samples <- filter(samples, diagnosis %in% c("nTD", "TD")) %>% droplevels()
  D0 <- samples %>% filter(time_point == "Baseline")
  
 samples %>% tabyl(study_arm, time_point) %>% adorn_totals("row")
  
  #samples <- filter(samples, sequence_pool == "original")
  #contamids <- samples[duplicated(samples$lab_id), ]
  
  TD <- filter(samples, diagnosis == "1")
  nTD <- filter(samples, diagnosis == "0")
  table(TD$time_point3)
  
  
  # Start - Tyger batch correction ---------------------------------------------
  
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
  
  # remove d7 from merged - update 10/nov - keep d7 in
  # idx <- (dge$samples$time_point != "D7")
  # dge <- dge[, idx]
  # dim(dge)
  
  
  # Update D7 var to include TD
  table(dge$samples$time_point)
  
  dge$samples<- dge$samples %>% 
    mutate(time_point = ifelse(time %in% c(5:9) & (time_point3 == "TD"), "D7", time_point))
  
  dge$samples <- dge$samples %>% 
    mutate(chall_vax = ifelse(str_detect(study_arm, "TN|CTRL"), "T", study_arm))
  
  save.image(file = "vt_bc_d7.RData")
  load(file = "vt_bc_d7.RData")
 
  # final correction -----------------------------------------------------------
  table(dge$samples$batch, dge$samples$time_point)
  # filter contaminated samples, recode to b1
  dge$samples <- filter(dge$samples, contaminated_sample == "No")

  dge$samples <- dge$samples %>% 
    mutate(batch2 = ifelse(batch == 2, 1, batch))
  
  idx <- (dge$samples$time_point %!in% c("D01","TD"))
  dge <- dge[,idx]
  table(dge$samples$time_point, dge$samples$batch2)
  
  covar <- model.matrix(~chall_vax, data = dge$samples)
  dge_c <- sva::ComBat_seq(counts=dge$counts,
                            batch=dge$samples$batch,
                            group=dge$samples$time_point,
                            covar_mod = covar) 
   
  # TPM
  combat <- DGEList(counts = dge_c, samples = dge$samples, genes = vast$genes)
  combat$combat1_counts  <- dge$counts
  combat$combat2_counts <- dge_c
  
  
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
  
  #Add gene info
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','start_position','end_position'),mart = ensembl) %>% rename(gene_id = ensembl_gene_id)
combat$genes <- left_join(combat$genes, genes, by = c("gene_id"))
  
saveRDS(combat, file = "combat_final_d7.RDS")

# Summary data -------------

tyg_sum <- table(tyger$metadata$time_point3, tyger$metadata$seq_run)
tyg_sum

length(unique(tyger$metadata$lab_id)) # = 40

tyg_sum2 <- table(tyger$metadata$time_point3, tyger$metadata$strain)
tyg_sum2

vast_sum <- table(vast$meta_data$time_point3, vast$meta_data$vaccine_y)
vast_sum

length(unique(vast$meta_data$lab_id)) #103
ncol(comb_t)


