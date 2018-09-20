
#p_values_fp = '../../2018WGS/Data/pooled_GWAS/PresAbs_based_pseudoCount/All_chr_positive_pvals_PresAbs.csv'
#p_values_fp = '../../2018WGS/Data/pooled_GWAS/PresAbs_based_rmMissingGT/Presabs_corrected_pvals.csv'
#p_values_fp = '../../2018WGS/Data/pooled_GWAS/Count_based_rmMissingGT/All_chr_positive_pvals_countbased.csv'

#p_values_fp = '../../2018WGS/Data/pooled_GWAS/HLI_logitRegression_test/all_p_vals_count-based.csv'

# p_values_fp = '../../2018WGS/Data/MismatchEncoded/DR_genotype_encoding//p_values/all_p_values_agvhd24_DR_genotype_encoding.csv'

p_values_fp = 'Documents/NGSProject/2018WGS/Data/GWASH/unimputed/cleaned_BMT_cases/DR_genotypes/logitRegression_p_values_agvhd34_DR_genotype_encoding.csv'

original_pvals = read.csv(p_values_fp, header = F, col.names = c('SNP', 'p_value'), stringsAsFactors = F)

positive_pvals = original_pvals[original_pvals$p_value >= 0, ]
#positive_pvals$p_value = positive_pvals$p.values
map_file <- 'Documents/NGSProject/2018WGS/Data/GWASH/unimputed/cleaned_BMT_cases/cleaned_BMT_cases.bim'
SNPmap <- read.table(map_file, col.names = c('chromosome', 'SNP', 'DIST', 'POS', 'ALT1', 'ALT2'), stringsAsFactors = F)

merged_p_val <- merge(positive_pvals, SNPmap, by = 'SNP', sort = F)

original_pvals = merged_p_val

original_pvals$chromosome <- sapply(1:dim(original_pvals)[1], function(x) paste0('chr', original_pvals$chromosome[x]))

original_pvals$log_p = -log10(original_pvals$p_value)

qqnorm(original_pvals$log_p, main = "Normal Q-Q plot", 
       xlab = expression(Expected ~ ('-'~'log'[10] ~ italic(p))),
       ylab = expression(Observed ~ ('-'~'log'[10] ~ italic(p))))
qqline(original_pvals$log_p, col = 2)

original_pvals$bonferroni = p.adjust(original_pvals$p_value, method = 'b')
original_pvals$holm = p.adjust(original_pvals$p_value, method = 'holm')

original_pvals$hochberg = p.adjust(original_pvals$p_value, method = 'hochberg')
original_pvals$hommel = p.adjust(original_pvals$p_value, method = 'hommel')

original_pvals$BH = p.adjust(original_pvals$p_value, method = 'BH')
original_pvals$BY = p.adjust(original_pvals$p_value, method = 'BY')

original_pvals$fdr = p.adjust(original_pvals$p_value, method = 'fdr')

original_pvals$log_fdr = -log10(original_pvals$fdr)
original_pvals$log_bonferroni = -log10(original_pvals$bonferroni)
original_pvals$log_holm = -log10(original_pvals$holm)
original_pvals$log_hochberg = -log10(original_pvals$hochberg)
original_pvals$log_BH = -log10(original_pvals$BH)
original_pvals$log_BY = -log10(original_pvals$BY)

#p_values_fp2 = '../../2018WGS/Data/MismatchEncoded/DR_genotype_encoding/p_values/'
p_values_fp2 = 'Documents/NGSProject/2018WGS/Data/GWASH/unimputed/cleaned_BMT_cases/DR_genotypes/'

write.csv(original_pvals, file = paste0(p_values_fp2, '/aGVHD34_all_corrected_pvals_DR_genotype_encoding.csv'), row.names = FALSE)
#save(original_pvals, file=paste0(p_values_fp2, 'aGVHD24_all_corrected_pvals.RData'))

# write.csv(original_pvals, file = 'Presabs_corrected_pvals.csv')
# save(original_pvals, file = 'Presabs_corrected_pvals.RData')

## 
corrected_pvals_fp = '../../2018WGS/Data/pooled_GWAS/PresAbs_based/All_chr_corrected_pvals_Presabsence.rdata'
load(corrected_pvals_fp)
##
X = original_pvals$p_value
Y = cbind(original_pvals$bonferroni,
          original_pvals$holm,
          original_pvals$hochberg,
          original_pvals$BH,
          original_pvals$BY,
          original_pvals$fdr)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:7,
        lty=1,
        lwd=2)

legend('bottomright', 
       legend = c("Bonferroni",  "Holm", "Hochberg", "Hommel", "BH", "BY", "fdr"), 
       col = 1:7, 
       cex = 1,    
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

write.csv(original_pvals, file = "../../2018WGS/Data/pooled_GWAS/PresAbs_based///Presabs_corrected_pvals.csv")
