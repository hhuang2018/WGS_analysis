source('util.R', echo = FALSE)

Gene_stats_table <- read.csv("../WW_MiHA/All.IDH1.csv", header = F)
colnames(Gene_stats_table) <- c("GroupType", "SubjectID", "CHROM", "POS", "DD", "REF", "ALT", "QUAL", "VariantType", "MODIFIER", "GeneName", "AssembleID", "TranscriptID")

Num_subjects <- length(unique(Gene_stats_table$SubjectID))

mc_stats <- data.frame( GroupType = character(Num_subjects), 
                        SubjectID = character(Num_subjects),
                        
                        three_prime_UTR_variant = numeric(Num_subjects),
                        five_prime_UTR_premature = numeric(Num_subjects),
                        five_prime_UTR_variant =  numeric(Num_subjects),
                        downstream_gene_variant =  numeric(Num_subjects),
                        intron_variant = numeric(Num_subjects),
                        missense_variant = numeric(Num_subjects),
                        splice_region_variant = numeric(Num_subjects),
                        synonymous_variant = numeric(Num_subjects),
                        upstream_gene_variant = numeric(Num_subjects),
                        
                        all_variants =  numeric(Num_subjects),
                        stringsAsFactors = F
                        )
unique_subjectID <- unique(Gene_stats_table$SubjectID)
for(id in 1:Num_subjects){
  
  subIndex <- which(Gene_stats_table$SubjectID %in% unique_subjectID[id])
  
  mc_stats$GroupType[id] <- as.character(Gene_stats_table$GroupType[subIndex[1]])
  mc_stats$SubjectID[id] <- as.character(unique_subjectID[id])
  mc_stats$all_variants[id] <- length(subIndex)
  
  variants_stats <- as.data.frame(table(Gene_stats_table$VariantType[subIndex]))
  
  mc_stats[id, 3:11] <- variants_stats$Freq
  
}


library(beeswarm)
library(ggplot2)
library(plyr)
beeswarm_all_var <- beeswarm(all_variants ~ GroupType, data = mc_stats, 
                            method = 'swarm', # 'square', 'hex', 'center', 'swarm'
                            spacing = .5)[, c(1, 2, 6)]
colnames(beeswarm_all_var) <- c("x", "y", "GVHD") 

ggplot(beeswarm_all_var, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#All mutated variants in IDH1")) + 
  geom_point(aes(colour = GVHD), size = 2) +
  scale_colour_manual(values = c("#D55E00", "#0072B2")) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non GVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_all_var$x, 1, round)), outlier.shape = NA, alpha = 0)

## exon

beeswarm_missense_variant <- beeswarm(missense_variant ~ GroupType, data = mc_stats, 
                             method = 'swarm', # 'square', 'hex', 'center', 'swarm'
                             spacing = .2)[, c(1, 2, 6)]
colnames(beeswarm_missense_variant) <- c("x", "y", "GVHD") 

ggplot(beeswarm_missense_variant, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Missense mutations in IDH1")) + 
  geom_point(aes(colour = GVHD), size = 2) +
  scale_colour_manual(values = c("#D55E00", "#0072B2")) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non GVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_missense_variant$x, 1, round)), outlier.shape = NA, alpha = 0)

t.test(mc_stats$missense_variant[mc_stats$GroupType=="a"], mc_stats$missense_variant[mc_stats$GroupType=="n"])
t.test(mc_stats$intron_variant[mc_stats$GroupType=="a"], mc_stats$intron_variant[mc_stats$GroupType=="n"])
t.test(mc_stats$all_variants[mc_stats$GroupType=="a"], mc_stats$all_variants[mc_stats$GroupType=="n"])
t.test(mc_stats$five_prime_UTR_premature[mc_stats$GroupType=="a"], mc_stats$five_prime_UTR_premature[mc_stats$GroupType=="n"]) #  p-value = 0.08326
t.test(mc_stats$downstream_gene_variant[mc_stats$GroupType=="a"], mc_stats$downstream_gene_variant[mc_stats$GroupType=="n"]) # p-value = 0.9669
t.test(mc_stats$splice_region_variant[mc_stats$GroupType=="a"], mc_stats$splice_region_variant[mc_stats$GroupType=="n"]) # p-value = 0.3175

t.test(mc_stats$upstream_gene_variant[mc_stats$GroupType=="a"], mc_stats$upstream_gene_variant[mc_stats$GroupType=="n"]) # p-value = 0.086
