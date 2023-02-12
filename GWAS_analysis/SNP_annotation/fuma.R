#Fuma analysis

setwd('~/GWAS_22/new_gwas/Post-imp/Merged/Assoc/annotation')

#sort col names for FUMA
#SNP | snpid | markername | rsID: rsID
#CHR | chromosome | chrom: chromosome
#BP | pos | position: genomic position (hg19)
#A1 | effect_allele | allele1 | alleleB: affected allele (ALT)
#A2 | non_effect_allele | allele2 | alleleA: another allele (REF)
#P | pvalue | p-value | p_value | pval: P-value (Mandatory)
#OR: Odds Ratio
#Beta | be: Beta #log_SE
#SE: Sta

colnames(T_gwas) <- c("CHR", "BP", "SNP", "REF", "A2", "A1","FIRTH?","TEST","OBS_CT","OR", "Beta","Z_STAT","P","ERRCODE")

fuma <- select(T_gwas, SNP, CHR, BP, A1, A2, P, OR, Beta)

write.table(fuma, file = "fuma.txt", sep = "\t", quote = F, col.names = T, row.names = F)

library(data.table)
list <- fread("fuma.txt")

#Also format for SNP nexus


#chr4 100000 100001 (1based)

lift <- select(fuma, CHR, BP) %>% mutate(BP_1 = BP+1)
chr <- rep("chr", 4942013)
lift2 <- cbind(chr, lift)
lift2$chrN <- paste(lift2$chr, lift2$CHR, sep = "")
BED <- select(lift2, chrN, BP, BP_1)

write.table(BED, file = "~/GWAS_22/new_gwas/Post-imp/Merged/Assoc/annotation/bed.txt", sep = "\t", quote = F, col.names = F, row.names = F)

system("awk '{print "chr" $1, $4 -1, $4, $2 }' P1Tyger_rm.bim > P1Tyger_rm_lift")

# ectract mapped variants snp IDs
system("awk '{print $4}' P1Tyger_38.bed > P1Tyger_38_rm.snps")
# ectract updated positions
system("awk '{print $4, $3}' P1Tyger_38.bed > P1Tyger_38_rm.pos")

#update plink file with SNPs that we could extract, and update the snp coordinates using update map function
system("./plink --bfile P1Tyger_rm --extract P1Tyger_38_rm.snps --update-map P1Tyger_38_rm.pos --make-bed --out P1Tyger_remap")


#fuma_sig <- select(sig, SNP, CHR, BP, A1, A2, P, OR, Beta)
#write.table(fuma_sig, file = "fuma_sig.txt", sep = "\t", quote = F, col.names = T, row.names = F)