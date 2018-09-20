source("utils/util.R")
# library(BSgenome.Hsapiens.UCSC.hg38)
chr.len <- c(248956422, # chr1
                242193529, # chr2
                198295559, # chr3
                190214555, # chr4
                181538259, # chr5
                170805979, # chr6
                159345973, # chr7
                145138636, # chr8
                138394717, # chr9
                133797422, # chr10
                135086622, # chr11
                133275309, # chr12
                114364328, # chr13
                107043718, # chr14
                101991189, # chr15
                90338345,  # chr16
                83257441,  # chr17
                80373285,  # chr18
                58617616,  # chr19
                64444167,  # chr20
                46709983,  # chr21
                50818468   # chr22
)

#source('utils/util.R', echo = FALSE)

# ID table
# load("../Data/ID_table.RData")
HLA_typings <- read.csv(file = "../../2018WGS/R/R_Output/available_cases_HLA.csv", stringsAsFactors = F)
Outcome <- read.csv(file = "../../2018WGS/R/R_Output/available_cases.csv")

load("../Data/GRCh38_gene_list.RData")

IBD_file_dir <- "../../2018WGS/Data/GWASH/IBD/"
IBD_reformat_output <- "../../2018WGS/Data/GWASH/IBD/"

# IBD_files <- list.files(IBD_file_dir, pattern = "\\.gt.ibd$")
# num_files <- length(IBD_files)
IBD_files <- "chr6ibd.ibd"
fid <- 1
# for(fid in 1:num_files){
  ptm <- proc.time() 
  
  IBD_table <- read.table(file = paste0(IBD_file_dir, IBD_files[fid]), stringsAsFactors = F)
  
  colnames(IBD_table) <- c("SampleID1", "HapID1", "SampleID2", "HapID2", "Chr", "StartInd", "EndInd", "LOD")
  
  num_rows <- dim(IBD_table)[1]
  
  new_IBD_table <- data.frame(SampleID1 = character(num_rows),
                              HapID1 = numeric(num_rows),
                              SampleID2 = character(num_rows),
                              HapID2 = numeric(num_rows),
                              Chr = character(num_rows),
                              StartID = numeric(num_rows),
                              EndID = numeric(num_rows),
                              LOD = numeric(num_rows),
                              Matched = numeric(num_rows),
                              GeneNames = character(num_rows),
                              stringsAsFactors = FALSE)
  
  for (id in 1: num_rows){
    
    #SampleID <- paste0(ID_table[which(ID_table$SeqID == IBD_table$SampleID1[id]), 1:3], collapse = ".")
    new_IBD_table$SampleID1[id] <- IBD_table$SampleID1[id]
    new_IBD_table$HapID1[id] <- IBD_table$HapID1[id]
    
    #SampleID <- paste0(ID_table[which(ID_table$SeqID == IBD_table$SampleID2[id]), 1:3], collapse = ".")
    new_IBD_table$SampleID2[id] <- IBD_table$SampleID2[id]
    new_IBD_table$HapID2[id] <- IBD_table$HapID2[id]
    
    new_IBD_table$Chr[id] <- as.character(IBD_table$Chr[id])
    new_IBD_table$StartID[id] <- IBD_table$StartInd[id]
    new_IBD_table$EndID[id] <- IBD_table$EndInd[id]
    new_IBD_table$LOD[id] <- IBD_table$LOD[id]
    
    if(nchar(new_IBD_table$SampleID1[id]) == 11 & nchar(new_IBD_table$SampleID2[id]) == 9){# donor
      
      idx = which(HLA_typings$nmdp_id == as.numeric(gsub("-","",new_IBD_table$SampleID1[id])))
      
      if (HLA_typings$nmdp_rid[idx] == as.numeric(gsub("-","",new_IBD_table$SampleID2[id]))){
        new_IBD_table$Matched[id] = 1
      }
      
    }else if(nchar(new_IBD_table$SampleID1[id]) == 9 & nchar(new_IBD_table$SampleID2[id]) == 11) { # recipient
      
      idx = which(HLA_typings$nmdp_rid == as.numeric(gsub("-","",new_IBD_table$SampleID1[id])))
      
      if (HLA_typings$nmdp_id[idx] == as.numeric(gsub("-","",new_IBD_table$SampleID2[id]))){
        new_IBD_table$Matched[id] = 1
      }
      
    }
      
    
    new_IBD_table$GeneNames[id] <- get.GeneNames(chrom = new_IBD_table$Chr[id], 
                                                 startPos = new_IBD_table$StartID[id],
                                                 endPos = new_IBD_table$EndID[id],
                                                 GeneList = GRCh38_gene_list)
    
    
    
  }
  
  file.out.RData <- gsub(".gt.ibd", "_IBD.RData", IBD_files[fid])
  file.out.csv <- gsub(".gt.ibd", "_IBD.csv", IBD_files[fid])
  save(new_IBD_table, file = paste0(IBD_reformat_output, file.out.RData))
  write.csv(new_IBD_table, file = paste0(IBD_reformat_output, file.out.csv), row.names = FALSE)
  
  cat(paste0("IBD table for ", IBD_files[fid]," has been converted and saved!\n"))
  print(proc.time() - ptm)
  cat("\n")
# }

############ Plot
library(ggplot2)
library(grid)

load("../Output/IBDseq/summary/ARS_MHC_chr6_WGS_IBD_04062017.RData")
load("../Output/IBDseq/OriginalVCFIBD_chr6/ARS_MHC_chr6_WGS_IBD_0501.RData")
# load("../Output/IBDseq/summary/all_chromosome_stats_1030.RData")
# load(paste0(IBD_reformated_dir, "summaryIBDseq_summary_chr_6.RData"))

shared_DNA_theory <- data.frame(percent = c(0.5, 0.25, 0.125),
                                groups = c("parent-child/siblings", "uncle/aunt-niece/nephew", "1st Cousin"))

aGVHD_IDs_ARS = grep('a.', Matched_ARS_length$SampleID1)
aGVHD_IDs_MHC = grep('a.', Matched_group_length$SampleID1)
aGVHD_IDs_chr6 = grep('a.', Matched_high_pert$SampleID1)
aGVHD_IDs_wgs = grep('a.', all_matched_Percent$SampleID1)

nGVHD_IDs_ARS = grep('n.', Matched_ARS_length$SampleID1)
nGVHD_IDs_MHC = grep('n.', Matched_group_length$SampleID1)
nGVHD_IDs_chr6 = grep('n.', Matched_high_pert$SampleID1)
nGVHD_IDs_wgs = grep('n.', all_matched_Percent$SampleID1)


mc_length <- data.frame(Proportion = c(Matched_ARS_length$Percent[aGVHD_IDs_ARS], Matched_ARS_length$Percent[nGVHD_IDs_ARS],     # ARS
                                       Matched_group_length$Percent[aGVHD_IDs_MHC], Matched_group_length$Percent[nGVHD_IDs_MHC], # MHC
                                       Matched_high_pert$Percent[aGVHD_IDs_chr6], Matched_high_pert$Percent[nGVHD_IDs_chr6], # Chr 6
                                       all_matched_Percent$IBDPercent[aGVHD_IDs_wgs], all_matched_Percent$IBDPercent[nGVHD_IDs_wgs]), # WGS
                        Group = c(rep("aGVHD group - ARS", length(aGVHD_IDs_ARS)), rep("nonGVHD group - ARS", length(nGVHD_IDs_ARS)),
                                  rep("aGVHD group - MHC", length(aGVHD_IDs_ARS)), rep("nonGVHD group - MHC", length(nGVHD_IDs_ARS)),
                                  rep("aGVHD group - Chr6", length(aGVHD_IDs_ARS)), rep("nonGVHD group - Chr6", length(nGVHD_IDs_ARS)),
                                  rep("aGVHD group - WGS", length(aGVHD_IDs_wgs)), rep("nonGVHD group - WGS", length(nGVHD_IDs_wgs))),
                        stringsAsFactors = F)

mc_length$Group <- factor(mc_length$Group, 
                          levels = c("aGVHD group - ARS", "nonGVHD group - ARS", 
                                     "aGVHD group - MHC", "nonGVHD group - MHC", 
                                     "aGVHD group - Chr6", "nonGVHD group - Chr6", 
                                     "aGVHD group - WGS", "nonGVHD group - WGS"))

mainplot <- ggplot(mc_length, aes(x = Group , y = Proportion, fill = Group)) + 
  geom_boxplot(outlier.colour = "gray", outlier.size = 0) +
  guides(fill=FALSE)  +
  geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous("Percentage of IBD segments between pairs") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme_bw()
# coord_flip() +
# ggtitle(paste0("Length (%) of IBD segments \non Chromosome ", chr))


mc_length_chr6only <- data.frame(Proportion = c(Matched_high_pert$Percent[aGVHD_IDs_chr6], Matched_high_pert$Percent[nGVHD_IDs_chr6]),  # Chr 6
                                 Group = c(rep("aGVHD group", length(aGVHD_IDs_ARS)), rep("nonGVHD group", length(nGVHD_IDs_ARS))),
                                 stringsAsFactors = F)
subplot_chr6 <- ggplot(mc_length_chr6only, aes(x = Group , y = Proportion, fill = Group)) + 
  geom_boxplot() +
  guides(fill=FALSE) +
  # geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous(element_blank()) + scale_x_discrete(element_blank()) + 
  ggtitle('Chromosome 6') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  theme_bw()

vp_chr6 <- viewport(width = 0.2, height = 0.35, 
                    x = 0.73,
                    y = unit(0.75, "npc"), just = "right")

mc_length_WGSonly <- data.frame(Proportion = c(all_matched_Percent$IBDPercent[aGVHD_IDs_wgs], all_matched_Percent$IBDPercent[nGVHD_IDs_wgs]), # WGS
                                Group = c(rep("aGVHD group", length(all_matched_Percent$IBDPercent)), rep("nonGVHD group", length(all_random_Percent$IBDPercent))),
                                stringsAsFactors = F)

subplot_wgs <- ggplot(mc_length_WGSonly, aes(x = Group , y = Proportion, fill = Group)) + 
  geom_boxplot() +
  guides(fill=FALSE) +
  # geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous(element_blank()) + scale_x_discrete(element_blank()) +
  ggtitle('Whole Genome') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  theme_bw()

vp_wgs <- viewport(width = 0.2, height = 0.35, 
                   x = 0.95,
                   y = unit(0.75, "npc"), just = "right")

print(mainplot)

print(subplot_chr6, vp = vp_chr6)
print(subplot_wgs, vp = vp_wgs)
# theme_set(theme_bw())

################
# ARS and MHC
################
mc_length_ARSonly <- data.frame(Proportion = c(Matched_ARS_length$Percent[aGVHD_IDs_ARS], Matched_ARS_length$Percent[nGVHD_IDs_ARS]),  # ARS
                                Group = c(rep("aGVHD group - ARS", length(aGVHD_IDs_ARS)), rep("nonGVHD group - ARS", length(nGVHD_IDs_ARS))),
                                stringsAsFactors = F)
ggplot(mc_length_ARSonly, aes(x = Group , y = Proportion, fill = Group)) + 
  geom_boxplot(outlier.colour = 'gray', outlier.alpha = 0.2) +
  guides(fill=FALSE) +
  #geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous(limits = c(0, 1)) + scale_x_discrete(element_blank()) +
  theme_light()

mc_length_MHConly <- data.frame(Proportion = c(Matched_group_length$Percent[aGVHD_IDs_MHC], Matched_group_length$Percent[nGVHD_IDs_MHC]),  # MHC
                                Group = c(rep("aGVHD group - MHC", length(aGVHD_IDs_MHC)), rep("nonGVHD group - MHC", length(nGVHD_IDs_MHC))),
                                stringsAsFactors = F)
ggplot(mc_length_MHConly, aes(x = Group , y = Proportion, fill = Group)) + 
  geom_boxplot(outlier.colour = 'gray', outlier.alpha = 0.3) +
  guides(fill=FALSE) +
  #geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous(element_blank()) + scale_x_discrete(element_blank()) + 
  theme_light()

######
# chr6 and WGS
###### 
mc_length_chr6only <- data.frame(Proportion = c(Matched_high_pert$Percent[aGVHD_IDs_chr6], Matched_high_pert$Percent[nGVHD_IDs_chr6]),  # Chr 6
                                 Group = c(rep("aGVHD group - Chr6", length(aGVHD_IDs_chr6)), rep("nonGVHD group - Chr6", length(nGVHD_IDs_chr6))),
                                 stringsAsFactors = F)
ggplot(mc_length_chr6only, aes(x = Group , y = Proportion, fill = Group)) + 
  geom_boxplot(outlier.colour = 'gray', outlier.alpha = 0.2) +
  guides(fill=FALSE) +
  # geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous(limits = c(0, .52)) + scale_x_discrete(element_blank()) + 
  theme_light()

mc_length_WGSonly <- data.frame(Proportion = c(all_matched_Percent$IBDPercent[aGVHD_IDs_wgs], all_matched_Percent$IBDPercent[nGVHD_IDs_wgs]), # WGS
                                Group = c(rep("aGVHD group - WGS", length(aGVHD_IDs_wgs)), rep("nonGVHD group - WGS", length(nGVHD_IDs_wgs))),
                                stringsAsFactors = F)

ggplot(mc_length_WGSonly, aes(x = Group , y = Proportion, fill = Group)) + 
  geom_boxplot(outlier.colour = 'gray', outlier.alpha = 0.2) +
  guides(fill=FALSE) +
  # geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous(limits = c(0, .52)) + scale_x_discrete(element_blank()) + 
  theme_light()
##########
# Remove outliers
##########
# compute lower and upper whiskers
# ylim1 = boxplot.stats(df$y)$stats[c(1, 5)]

mc_length_ARSonly <- data.frame(Proportion = c(Matched_ARS_length$Percent, Random_ARS_length$Percent),  # ARS
                                Group = c(rep("Matched ARS", length(Matched_ARS_length$Percent)), rep("Random ARS", length(Random_ARS_length$Percent))),
                                stringsAsFactors = F)
ggplot(mc_length_ARSonly, aes(x = Group , y = Proportion, fill = Group)) + 
  geom_boxplot(outlier.size = NA) +
  guides(fill=FALSE) +
  #geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous(limits = c(0, 1)) + scale_x_discrete(element_blank()) +
  theme_light()

#ylim1 = boxplot.stats(mc_length_ARSonly$Proportion)$stats[c(1, 5)]

#ARS_only_pp + coord_cartesian(ylim = ylim1*1.05)
#
mc_length_MHConly <- data.frame(Proportion = c(Matched_group_length$Percent, Random_group_length$Percent),  # MHC
                                Group = c(rep("Matched MHC", length(Matched_group_length$Percent)), rep("Random MHC", length(Random_group_length$Percent))),
                                stringsAsFactors = F)
ggplot(mc_length_MHConly, aes(x = Group , y = Proportion, fill = Group)) + 
  geom_boxplot(outlier.size = NA) +
  guides(fill=FALSE) +
  # geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous(element_blank()) + scale_x_discrete(element_blank()) +
  theme_light()

#ylim2 = boxplot.stats(mc_length_MHConly$Proportion)$stats[c(1, 5)]

#MHConly_pp + coord_cartesian(ylim = ylim2*1.05)

######
# chr6 and WGS
###### 
mc_length_chr6only <- data.frame(Proportion = c(Matched_high_pert$Percent, Random_high_pert$Percent),  # Chr 6
                                 Group = c(rep("Matched Chr6", length(Matched_high_pert$Percent)), rep("Random Chr6", length(Random_high_pert$Percent))),
                                 stringsAsFactors = F)
ggplot(mc_length_chr6only, aes(x = Group , y = Proportion, fill = Group)) + 
  #geom_boxplot() +
  geom_boxplot(outlier.size = NA)+
  guides(fill=FALSE) +
  #geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous(limits = c(0, .06)) + scale_x_discrete(element_blank()) +
  theme_light()


mc_length_WGSonly <- data.frame(Proportion = c(all_matched_Percent$IBDPercent, all_random_Percent$IBDPercent), # WGS
                                Group = c(rep("Matched WGS", length(all_matched_Percent$IBDPercent)), rep("Random WGS", length(all_random_Percent$IBDPercent))),
                                stringsAsFactors = F)

ggplot(mc_length_WGSonly, aes(x = Group , y = Proportion, fill = Group)) + 
  geom_boxplot(outlier.size = NA) +
  guides(fill=FALSE) +
  #geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  # geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous(limits = c(0, .02)) + scale_x_discrete(element_blank()) +
  theme_light() +
  scale_fill_manual(values=c("#af8dc3", "#7fbf7b"))

