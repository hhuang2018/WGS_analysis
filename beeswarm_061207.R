library(ggplot2)
# miha_table <- read.table("../Boxplot205pairs.txt", header = T) ### Old data -- wrong
known_miha_freq <- read.delim("../WW_MiHA/Restricted_known_MiHAs.txt", header = F, stringsAsFactors = F)
colnames(known_miha_freq) <- c("GroupType", "GroupID",  "HLA_type", "SNP", "CHROM", "REF", "ALT")

known_unrestricted_MiHA_table <- read.csv("../ClinVar/unRestrictedMiHAs.csv", header = F, stringsAsFactors = F)
colnames(known_unrestricted_MiHA_table) <- c("GroupType", "GroupID",  "HLA_type", "SNP", "CHROM", "REF", "ALT")

class(known_miha_freq$GroupID) <- "charater"
known_Miha_stats <- as.data.frame(table(known_miha_freq$GroupID), stringsAsFactors = F)
colnames(known_Miha_stats) <- c("GroupID", "KnownRestrictedMiHAs")

class(known_unrestricted_MiHA_table$GroupID) <- "charater"
known_unrestricted_stats <- as.data.frame(table(known_unrestricted_MiHA_table$GroupID),stringsAsFactors = F)
colnames(known_unrestricted_stats) <- c("GroupID", "KnownUnRestrictedMiHAs")

############ check D-R missense mismatches
load("../Data/HLI_available_pairs_dis_table.RData")
Missense_mismatche_fp <- "../ClinVar/D_R_Missense_Mismatches/"
Missense_mismatche_files <- list.files(Missense_mismatche_fp)
num_files <- length(Missense_mismatche_files)
# sink(file = "../ClinVar/D_R_Missense_Mismatches/Missense_numbers.txt")
missense_table <- data.frame()
All_stats <- data.frame()
# row_names <- character(num_files)
for(id in 1:num_files){
  
  missense_table <- read.table(file = paste0(Missense_mismatche_fp, Missense_mismatche_files[id]))
  missense_table$V1 <- factor(missense_table$V1, levels = sapply(c(1:22, "X", "Y"), function(x) paste0("chr", x)))
  temp_stats <- as.data.frame(table(missense_table$V1))
  colnames(temp_stats) <- c("CHROM", Missense_mismatche_files[id])
  if(dim(All_stats)[1] ==0){
    
    All_stats <- temp_stats
    
  }else{
    
    All_stats <- merge(All_stats, temp_stats, by = "CHROM")
    
  }
  # row_names[id] <- Missense_mismatche_files[id]
  # print(Missense_mismatche_files[id])
  # print(table(missense_table$V1))
  
}
# colnames(All_stats) <- c("CHROM", row_names)
WGS_wXY <- colSums(All_stats[, -1])
WGS_wXY <- as.data.frame(WGS_wXY)
WGS_woXY <- colSums(All_stats[-c(23,24), -1])
WGS_woXY <- as.data.frame(WGS_woXY)
#### Make miha_table: num_misense_mismatch, num_unrestricted_MiHA, num_Restricted_MiHa
load("../Data/ID_table.RData")
HLI_metadata <- read.csv("../HLI_hla_mg_v3.csv", stringsAsFactors = F)

ID_table$caseNumber <- 0 
for(id in 1:dim(ID_table)[1]){
  
  if(ID_table$subjectType[id] == "R"){ # recipients
    
    r_ind <- which(HLI_metadata$nmdp_rid %in% ID_table$R_D_ID[id])
    if(length(r_ind) > 1) cat("id = ", id, "; r_ind = ", r_ind, "\n")
    
    ID_table$caseNumber[id] <- HLI_metadata$bmt_case_num[r_ind] 
    
  }else{  # donors
    
    d_ind <- which(HLI_metadata$nmdp_id %in% ID_table$R_D_ID[id])
    if(length(d_ind) > 1) cat("id = ", id, "; d_ind = ", d_ind, "\n")
    
    ID_table$caseNumber[id] <- HLI_metadata$bmt_case_num[d_ind] 
    
  }
  
}

n_cases <- 205
miha_table <- data.frame(groupID = character(n_cases),
                         GVHD = character(n_cases),
                         SEX = character(n_cases),
                         NumVar = numeric(n_cases),
                         HLA.Unrestricted = numeric(n_cases),
                         HLA.Restricted = numeric(n_cases),
                         stringsAsFactors = F)

for(id in 1:n_cases){
  
  miha_table$groupID[id] <- known_unrestricted_stats$GroupID[id]
  miha_table$HLA.Unrestricted[id] <- known_unrestricted_stats$KnownUnRestrictedMiHAs[id]
  
  temp_ID_table_ind <- which(ID_table$GroupID %in% miha_table$groupID[id])
  temp_HLI_ind <- which(HLI_metadata$bmt_case_num %in% unique(ID_table$caseNumber[temp_ID_table_ind]))
  
  miha_table$SEX[id] <- paste0(HLI_metadata$donor_sex[temp_HLI_ind], ">", HLI_metadata$rid_sex[temp_HLI_ind])
  
  miha_table$GVHD[id] <- if(unique(ID_table$Group[temp_ID_table_ind]) == "a") "aGVHD" else "nonGVHD"
  
  miha_table$NumVar[id] <- WGS_woXY[paste0(miha_table$groupID[id], ".txt"), ]
  
  temp_restricted_table_ind <- which(known_Miha_stats$GroupID %in% miha_table$groupID[id])
  if(length(temp_restricted_table_ind) > 0){
    
    miha_table$HLA.Restricted[id] <- known_Miha_stats$KnownRestrictedMiHAs[temp_restricted_table_ind]
    
  }
  
}

miha_table$SEX <- as.factor(miha_table$SEX)
miha_table$GVHD <- as.factor(miha_table$GVHD)
####### plots

p1 <- ggplot(miha_table, aes(x = factor(GVHD), y = HLA.Restricted))
p1 + geom_boxplot(aes(fill = GVHD)) + 
  ggtitle(paste0("IBD score (LOD) distribution on ", Region)) +
  # geom_density(stat = "density") +
  labs(x="Group", y = "LOD")

aGVHD_index <- intersect(which(miha_table$GVHD=="aGVHD"), which(miha_table$HLA.Restricted < 12))
nonGVHD_index <- intersect(which(miha_table$GVHD=="nonGVHD"), which(miha_table$HLA.Restricted < 15))
ttest <- t.test(miha_table[aGVHD_index, 5], miha_table[nonGVHD_index, 5])
ttest$p.value 

sort(miha_table[miha_table$GVHD=="aGVHD", 6])
sort(miha_table[miha_table$GVHD=="nonGVHD", 6])

library(beeswarm)

beeswarm(HLA.Restricted ~ GVHD, data = miha_table, pch = 16, cex = 1.5, col=c("#31a354", "#0072B2"),
         main = "HLA Restricted MiHA", xlab="")

beeswarm(HLA.Unrestricted ~ GVHD, data = miha_table, pch = 16, cex = 1.5, col=c("#31a354", "#0072B2"),
         main = "HLA Restricted MiHA", xlab="")

beeswarm(NumVar ~ GVHD, data = miha_table, pch = 16, cex = 1.5, col=c("#31a354", "#0072B2"),
         main = "HLA Restricted MiHA", xlab="", ylab="Number of missense variants")


beeswarm(HLA.Restricted ~ SEX, data = miha_table, pch = 16, cex = 1.5, col=c("#009E73", "#0072B2", "#D55E00", "#CC79A7"),
         main = "HLA Restricted MiHA", xlab="")

beeswarm(HLA.Unrestricted ~ SEX, data = miha_table, pch = 16, cex = 1.5, col=c("#009E73", "#0072B2", "#D55E00", "#CC79A7"),
         main = "HLA Unestricted MiHA", xlab="")


beeswarm(NumVar ~ SEX, data = miha_table, pch = 16, cex = 1.5, col=c("#009E73", "#0072B2", "#D55E00", "#CC79A7"),
         main = "R-D gender pair vs missense mismatches", xlab="", ylab="Number of missense variants")

#####3
# beeswarm <- beeswarm(time_survival ~ event_survival, 
#                      data = breast, method = 'swarm', 
#                      pwcol = ER)[, c(1, 2, 4, 6)]
# colnames(beeswarm) <- c("x", "y", "ER", "event_survival") 
library(ggplot2)
library(plyr)
# miha_table <- as.data.frame(miha_table)
beeswarm_gender <- beeswarm(NumVar ~ SEX, data = miha_table, 
                            method = 'swarm', # 'square', 'hex', 'center', 'swarm'
                            spacing = .5,
                            pwcol = GVHD)[, c(1, 2, 4, 6)]
colnames(beeswarm_gender) <- c("x", "y", "GVHD", "Sex") 

ggplot(beeswarm_gender, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Missense mismatches")) + 
  geom_point(aes(colour = GVHD), size = 2) +
  scale_colour_manual(values = c("#D55E00", "#0072B2")) + 
  scale_x_continuous(breaks = c(1:4), 
                     labels = c("F -> F", "F -> M", "M -> F", "M -> M"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_gender$x, 1, round)), outlier.shape = NA, alpha = 0)

#### HLA Restricted
beeswarm_resMiHA <- beeswarm(HLA.Restricted ~ SEX, data = miha_table, 
                             method = 'swarm', # 'square', 'hex', 'center', 'swarm'
                             spacing = .4,
                             pwcol = GVHD)[, c(1, 2, 4, 6)]
colnames(beeswarm_resMiHA) <- c("x", "y", "GVHD", "Sex") 

ggplot(beeswarm_resMiHA, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#HLA Restricted MiHA SNPs")) + 
  # geom_point(aes(colour = GVHD), size = 3, alpha = 0.8) +
  geom_jitter(width = 0.3, aes(color = GVHD), size = 2.5, alpha = 0.7) +
  scale_colour_manual(values = c("#D55E00", "#0072B2")) + 
  scale_x_continuous(breaks = c(1:4), 
                     labels = c("F -> F", "F -> M", "M -> F", "M -> M"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_resMiHA$x, 1, round)), outlier.shape = NA, alpha = 0)


#### HLA Unrestricted
beeswarm_unresMiHA <- beeswarm(HLA.Unrestricted ~ SEX, data = miha_table, 
                               method = 'swarm', # 'square', 'hex', 'center', 'swarm'
                               spacing = .7,
                               pwcol = GVHD)[, c(1, 2, 4, 6)]
colnames(beeswarm_unresMiHA) <- c("x", "y", "GVHD", "Sex") 

ggplot(beeswarm_unresMiHA, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#HLA Unrestricted MiHA SNPs")) + 
  # geom_point(aes(colour = GVHD), size = 3, alpha = 0.7) +
  geom_jitter(width = 0.3, aes(color = GVHD), size = 2.5, alpha = 0.7) +
  scale_colour_manual(values = c("#D55E00", "#0072B2")) + 
  scale_x_continuous(breaks = c(1:4), 
                     labels = c("F -> F", "F -> M", "M -> F", "M -> M"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_unresMiHA$x, 1, round)), outlier.shape = NA, alpha = 0)

##### Restricted MiHA
beeswarm_ResMiHA_all <- beeswarm(HLA.Restricted ~ GVHD, data = miha_table, 
                                 method = 'swarm',
                                 spacing = 0.5)[, c(1,2,6)]
colnames(beeswarm_ResMiHA_all) <- c("x", "y", "GVHD") 

ggplot(beeswarm_ResMiHA_all, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Mismatched HLA Restricted Known MiHA SNPs")) + 
  # geom_point(aes(colour = GVHD), size = 3, alpha = 0.5) +
  geom_jitter(width = 0.2, aes(color = GVHD), size = 3, alpha = 0.8) +
  scale_colour_manual(values = c("#D55E00", "#0072B2")) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non aGVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_ResMiHA_all$x, 1, round)), outlier.shape = NA, alpha = 0) +
  theme(legend.position = "none")

t.test(beeswarm_ResMiHA_all$y[which(beeswarm_ResMiHA_all$GVHD == "aGVHD")], beeswarm_ResMiHA_all$y[which(beeswarm_ResMiHA_all$GVHD == "nonGVHD")])

wilcox.test(beeswarm_ResMiHA_all$y[which(beeswarm_ResMiHA_all$GVHD == "aGVHD")], beeswarm_ResMiHA_all$y[which(beeswarm_ResMiHA_all$GVHD == "nonGVHD")])

wilcox.test(y ~ GVHD, beeswarm_ResMiHA_all)
###### Unrestricted MiHA

beeswarm_UnResMiHA_all <- beeswarm(HLA.Unrestricted ~ GVHD, data = miha_table,
                                   method = 'swarm',
                                   spacing = 1)[, c(1,2,6)]
colnames(beeswarm_UnResMiHA_all) <- c("x", "y", "GVHD") 

ggplot(beeswarm_UnResMiHA_all, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Mismatched HLA Unrestricted Known MiHA SNPs")) + 
  #geom_point(aes(colour = GVHD), size = 3, alpha = 0.7) +
  geom_jitter(width = 0.2, aes(color = GVHD), size = 3, alpha = 0.8) +
  scale_colour_manual(values = c("#D55E00", "#0072B2")) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non aGVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_UnResMiHA_all$x, 1, round)), outlier.shape = NA, alpha = 0) + 
  theme(legend.position = "none")

t.test(beeswarm_UnResMiHA_all$y[which(beeswarm_UnResMiHA_all$GVHD == "aGVHD")], beeswarm_UnResMiHA_all$y[which(beeswarm_UnResMiHA_all$GVHD == "nonGVHD")])

wilcox.test(beeswarm_UnResMiHA_all$y[which(beeswarm_UnResMiHA_all$GVHD == "aGVHD")], beeswarm_UnResMiHA_all$y[which(beeswarm_UnResMiHA_all$GVHD == "nonGVHD")])

wilcox.test(y ~ GVHD, beeswarm_UnResMiHA_all)
##### Number of Missense Mismatch
beeswarm_MissenseMis <-  beeswarm(NumVar ~ GVHD, data = miha_table,
                                  method = 'swarm',
                                  spacing = 0.6)[, c(1,2,6)]
colnames(beeswarm_MissenseMis) <- c("x", "y", "GVHD") 

ggplot(beeswarm_MissenseMis, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Missense Mismathced SNPs")) + 
  geom_point(aes(colour = GVHD), size = 3, alpha = 0.6) +
  scale_colour_manual(values = c("#D55E00", "#0072B2")) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non aGVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_MissenseMis$x, 1, round)), outlier.shape = NA, alpha = 0)+ 
  theme(legend.position = "none")


t.test(beeswarm_MissenseMis$y[which(beeswarm_MissenseMis$GVHD == "aGVHD")], beeswarm_MissenseMis$y[which(beeswarm_MissenseMis$GVHD == "nonGVHD")])

wilcox.test(beeswarm_MissenseMis$y[which(beeswarm_MissenseMis$GVHD == "aGVHD")], beeswarm_MissenseMis$y[which(beeswarm_MissenseMis$GVHD == "nonGVHD")])

wilcox.test(y ~ GVHD, beeswarm_MissenseMis)




######

########################
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()                  ## R version 3.0 or later
biocLite("Gviz")
biocLite("rtracklayer")
biocLite("trackViewer")
install.packages("devtools")

devtools::install_github("griffithlab/GenVisR")
#######################
## mutation plot
# library(Gviz)
# library(rtracklayer)
# library(trackViewer)
# library(GenVisR)
chrLengths <- c(248956422, # chr1
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

library(ggplot2)
known_MiHA_table <- read.table("../WW_MiHA/Known_MiHA_coordinates.txt", header = T)

load("../WW_MiHA/summary/Data/missense_summary_DtoR.RData")

all_MiHA_SNP <- data.frame()
Snp_pplots <- list()
for(id in 1:22){
  chr <- paste0("chr", id)
  
  aGVHD_Mutation <- aGVHD_SNP_list[[chr]]
  nGVHD_Mutation <- nGVHD_SNP_list[[chr]]
  
  known_MiHA_chr <- known_MiHA_table[which(known_MiHA_table$Chr %in% chr), ]
  known_MiHA_chr$NumDiff <- rep(0, length(known_MiHA_chr$Pos))
  names(known_MiHA_chr)[5] <- "POS"
  # known_MiHA_chr$POS <- factor(known_MiHA_chr$POS)
  # MiHA_SNP <- known_MiHA_chr$POS
  # known_MiHA_chr$group <- "aGVHD"
  
  mc <- data.frame(position = c(aGVHD_Mutation$POS, nGVHD_Mutation$POS),
                   frequency = c(aGVHD_Mutation$NumDiff, -nGVHD_Mutation$NumDiff),
                   group = c(rep("aGVHD", length(aGVHD_Mutation$POS)), rep("nGVHD", length(nGVHD_Mutation$POS))), 
                   stringsAsFactors = F)
  # mc$frequency[mc$position == aGVHD_Mutation$POS] <- aGVHD_Mutation$NumDiff
  # mc$frequency[mc$position == nGVHD_Mutation$POS] <- -nGVHD_Mutation$NumDiff
  # ggplot(aGVHD_Mutation, aes(x = POS, ymax = NumDiff, ymin = 0)) + 
  #   geom_linerange()
  if(min(mc$position) > 1){
    mc <- rbind(mc, data.frame(position = 1,
                               frequency = 0,
                               group = "nGVHD"))
  }
  
  if(max(mc$position) < chrLengths[id]){
    mc <- rbind(mc, data.frame(position = chrLengths[id],
                               frequency = 0,
                               group = "nGVHD"))
  }
  
  # mc$group[which(mc$position %in% known_MiHA_chr$POS)] <- "MiHA"
  
  # print(
  # p1 <- ggplot(mc, aes(x = position, ymax = frequency, ymin = 0, color=factor(group))) +
  #   geom_linerange() +
  #   geom_hline(yintercept = 0) +
  #   # ggtitle(chr) +
  #   theme_bw() +
  #   scale_x_continuous(breaks = round(seq(1, chrLengths[id], by = 2e7), 2),
  #                      labels= paste0(round(seq(1, chrLengths[id], by = 2e7)/1e6, 2), " Mb")) +
  #   scale_color_manual(values = c("#D55E00", "#0072B2")) +
  #   theme(axis.title.x=element_blank(),
  #         axis.text.x=element_blank(),
  #         axis.ticks.x=element_blank(),
  #         legend.position = "none")
  # #   ) 
  # Snp_pplots <- append(Snp_pplots, list(p1))
  ## MiHA SNPs by chromosome
  if(length(known_MiHA_chr$POS) >= 1){
    index <- which(mc$position %in% known_MiHA_chr$POS)
    if(length(index) > 0){
      MiHA_SNP <- lapply(1:length(index), function(x) cbind(mc[index[x], ], known_MiHA_chr[which(known_MiHA_chr$POS == mc$position[index[x]]), ]))
      MiHA_SNP <- do.call(rbind.data.frame, MiHA_SNP)
      
      names(MiHA_SNP)[3] <- "Group"
      
      all_MiHA_SNP <- rbind(all_MiHA_SNP, MiHA_SNP)
      # print(
      #   ggplot(MiHA_SNP, aes(x = MiHAs, ymax = frequency, ymin = 0, color = Group)) +
      #     geom_linerange(size = 3) +
      #     geom_hline(yintercept = 0) +
      #     ggtitle(paste0(chr, " - MiHAs")) +
      #     theme_bw() +
      #     xlab("") +
      # scale_fill_discrete(name = "Group")
      # )
    }
  }
  
}

multiplot(plotlist = Snp_pplots, cols = 1)

print(
  ggplot(all_MiHA_SNP, aes(x = MiHAs, ymax = frequency, ymin = 0, color = Group)) +
    geom_linerange(size = 3) +
    geom_hline(yintercept = 0) +
    ggtitle("Known MiHAs") +
    theme_bw() +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
)

write.csv(all_MiHA_SNP, file = "../Data/All_known_MiHAs_counts.csv", row.names = F)

unique_MiHAs <- as.character(unique(all_MiHA_SNP$MiHAs))
num_MiHAs <- length(unique_MiHAs)
LLR_MiHAs <- data.frame(name = character(num_MiHAs),
                        hla = character(num_MiHAs),
                        LLR = numeric(num_MiHAs),
                        stringsAsFactors = F)
all_MiHA_SNP$MiHAs <- as.character(all_MiHA_SNP$MiHAs)
all_MiHA_SNP$Gene <- as.character(all_MiHA_SNP$Gene)
all_MiHA_SNP$SNPs <- as.character(all_MiHA_SNP$SNPs)
for(miha_id in 1:num_MiHAs){
  
  index <- which(all_MiHA_SNP$MiHAs %in% unique_MiHAs[miha_id])
  
  aGVHD_count <- all_MiHA_SNP[index[which(all_MiHA_SNP[index, "Group"] == "aGVHD")], "frequency"]
  nGVHD_count <- abs(all_MiHA_SNP[index[which(all_MiHA_SNP[index, "Group"] == "nGVHD")], "frequency"])
  Gene_MiHA_name <- paste0(as.character(all_MiHA_SNP[index[which(all_MiHA_SNP[index, "Group"] == "aGVHD")], c("Gene", "MiHAs", "SNPs")]), collapse = "*")
  
  if(aGVHD_count > 0 |nGVHD_count >0){
    
    LLR_MiHAs$LLR[miha_id] <- log10(aGVHD_count/nGVHD_count)
    LLR_MiHAs$name[miha_id] <- Gene_MiHA_name
  }
  
}

mc_LLR <- LLR_MiHAs[order(LLR_MiHAs$LLR, decreasing = TRUE), ]
mc_LLR <- within(mc_LLR, name <- factor(name, levels=factor(mc_LLR$name)))

pp1 <- ggplot(all_MiHA_SNP, aes(x = MiHAs, ymax = frequency, ymin = 0, color = Group)) +
  geom_linerange(size = 3) +
  geom_hline(yintercept = 0) +
  ggtitle("Known MiHAs") +
  theme_bw() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pp2 <- ggplot(mc_LLR, aes(x = name, ymax = LLR, ymin = 0, color = "#D55E00")) +
  geom_linerange(size = 3) +
  geom_hline(yintercept = 0) +
  # ggtitle("Known MiHAs") +
  theme_bw() +
  xlab("") +
  # coord_flip() +
  # scale_colour(values = "#D55E00") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") 

multiplot(pp1, pp2, cols = 1)

# x <- sample.int(length(aGVHD_Mutation$POS), length(aGVHD_Mutation$POS))
# sample.gr <- GRanges("chr1", IRanges(aGVHD_Mutation$POS, width=1, names=paste0("POS", aGVHD_Mutation$POS)))
# features <- GRanges("chr1", IRanges(c(1, 501, 1001), 
#                                     width=c(120, 400, 405),
#                                     names=paste0("block", 1:3)))
# lolliplot(sample.gr, features)

##################################
# log ratio
##################################
library(ggplot2)
source("util.R")
known_MiHA_table <- read.table("../WW_MiHA/Known_MiHA_coordinates.txt", header = T)

load("../WW_MiHA/summary/Data/missense_summary_DtoR.RData")

head(aGVHD_SNP_list[[1]])

aGVHD_missense_snps <- do.call("rbind", aGVHD_SNP_list)
nGVHD_missense_snps <- do.call("rbind", nGVHD_SNP_list)

write.csv(aGVHD_missense_snps, file = "../Output/Missense_variant_stats/aGVHD_group_missense_mismatches.csv", row.names = F)
write.csv(nGVHD_missense_snps, file = "../Output/Missense_variant_stats/non_aGVHD_group_missense_mismatches.csv", row.names = F)



all_MiHA_SNP <- data.frame()
pplots <- list()
pplots2 <- list()
for(id in 1:22){
  chr <- paste0("chr", id)
  
  aGVHD_Mutation <- aGVHD_SNP_list[[chr]]
  nGVHD_Mutation <- nGVHD_SNP_list[[chr]]
  
  known_MiHA_chr <- known_MiHA_table[which(known_MiHA_table$Chr %in% chr), ]
  known_MiHA_chr$NumDiff <- rep(0, length(known_MiHA_chr$Pos))
  names(known_MiHA_chr)[5] <- "POS"
  # known_MiHA_chr$POS <- factor(known_MiHA_chr$POS)
  # MiHA_SNP <- known_MiHA_chr$POS
  # known_MiHA_chr$group <- "aGVHD"
  
  mc <- data.frame(position = c(aGVHD_Mutation$POS, nGVHD_Mutation$POS),
                   frequency = c(aGVHD_Mutation$NumDiff, -nGVHD_Mutation$NumDiff),
                   group = c(rep("aGVHD", length(aGVHD_Mutation$POS)), rep("nGVHD", length(nGVHD_Mutation$POS))), 
                   stringsAsFactors = F)
  inter_pos <- sort(intersect(aGVHD_Mutation$POS, nGVHD_Mutation$POS), decreasing = F)
  
  log_ratio <- data.frame(position = c(aGVHD_Mutation$POS, nGVHD_Mutation$POS), 
                          log_ratio = numeric(length(aGVHD_Mutation$POS) + length(nGVHD_Mutation$POS)),
                          group = character(length(aGVHD_Mutation$POS) + length(nGVHD_Mutation$POS)),
                          stringsAsFactors = F)
  aGVHD_only_id <- which(!(aGVHD_Mutation$POS %in% inter_pos))
  nGVHD_only_id <- which(!(nGVHD_Mutation$POS %in% inter_pos))
  
  log_ratio$log_ratio[aGVHD_only_id] <- aGVHD_Mutation$NumDiff[aGVHD_only_id]
  log_ratio$group[aGVHD_only_id] <- "aGVHD_only"
  log_ratio$log_ratio[nGVHD_only_id + length(aGVHD_Mutation$POS)] <- -nGVHD_Mutation$NumDiff[nGVHD_only_id]
  log_ratio$group[nGVHD_only_id + length(aGVHD_Mutation$POS)] <- "nGVHD_only"
  # intersect(aGVHD_Mutation$POS[aGVHD_only_id],  nGVHD_Mutation$POS[nGVHD_only_id])
  
  for(inter_id in 1:length(inter_pos)){
    
    aGVHD_inter_id <- which(aGVHD_Mutation$POS %in% inter_pos[inter_id])
    nGVHD_inter_id <- which(nGVHD_Mutation$POS %in% inter_pos[inter_id])
    
    log_ratio$log_ratio[aGVHD_inter_id] <- log10(aGVHD_Mutation$NumDiff[aGVHD_inter_id] / nGVHD_Mutation$NumDiff[nGVHD_inter_id])
    log_ratio$group[aGVHD_inter_id] <- "log_odds"
  }
  
  # mc$frequency[mc$position == aGVHD_Mutation$POS] <- aGVHD_Mutation$NumDiff
  # mc$frequency[mc$position == nGVHD_Mutation$POS] <- -nGVHD_Mutation$NumDiff
  # ggplot(aGVHD_Mutation, aes(x = POS, ymax = NumDiff, ymin = 0)) + 
  #   geom_linerange()
  if(min(log_ratio$position) > 1){
    log_ratio <- rbind(log_ratio, data.frame(position = 1,
                                             log_ratio = 0,
                                             group = "nGVHD_only", 
                                             stringsAsFactors = F))
  }
  
  if(max(log_ratio$position) < chrLengths[id]){
    log_ratio <- rbind(log_ratio, data.frame(position = chrLengths[id],
                                             log_ratio = 0,
                                             group = "nGVHD_only",
                                             stringsAsFactors = F))
  }
  
  # mc$group[which(mc$position %in% known_MiHA_chr$POS)] <- "MiHA"
  
  # print(
  p1 <- ggplot(log_ratio[which(log_ratio$group == "nGVHD_only" | log_ratio$group == "aGVHD_only"), ], aes(x = position, y = log_ratio, color=factor(group))) +
    geom_point(aes(shape = factor(group)), alpha = 1/2) +
    geom_hline(yintercept = 0) +
    #ggtitle(chr) +
    theme_bw() +
    scale_x_continuous(breaks = round(seq(1, chrLengths[id], by = 2e7), 2),
                       labels= paste0(round(seq(1, chrLengths[id], by = 2e7)/1e6, 2), " Mb")) +
    scale_color_manual(values = c("#009E73", "#0072B2")) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")
  
  p2 <-  ggplot(log_ratio[which(log_ratio$group == "log_odds"), ], aes(x = position, y = log_ratio)) +
    geom_point(shape = 16, color = "#D55E00", alpha = 1/4) +
    geom_hline(yintercept = 0) + 
    theme_bw() +
    scale_x_continuous(breaks = round(seq(1, chrLengths[id], by = 2e7), 2),
                       labels= paste0(round(seq(1, chrLengths[id], by = 2e7)/1e6, 2), " Mb")) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")
  
  pplots <- append(pplots, list(p1))
  pplots2 <- append(pplots2, list(p2))
  #   
  # if(length(known_MiHA_chr$POS) >= 1){
  #   index <- which(mc$position %in% known_MiHA_chr$POS)
  #   if(length(index) > 0){
  #     MiHA_SNP <- lapply(1:length(index), function(x) cbind(mc[index[x], ], known_MiHA_chr[which(known_MiHA_chr$POS == mc$position[index[x]]), ]))
  #     MiHA_SNP <- do.call(rbind.data.frame, MiHA_SNP)
  #     
  #     names(MiHA_SNP)[3] <- "Group"
  #     
  #     all_MiHA_SNP <- rbind(all_MiHA_SNP, MiHA_SNP)
  #     # print(
  #       ggplot(MiHA_SNP, aes(x = MiHAs, ymax = frequency, ymin = 0, color = Group)) +
  #         geom_linerange(size = 3) +
  #         geom_hline(yintercept = 0) +
  #         xlab("") 
  # 
  #   }
  #   
  # }
}

#eval(parse(text = paste0("multiplot( ", sapply(1:22, function(x) paste0("pplots[[", x,"]],"))," cols = 1)")))

multiplot(plotlist = pplots, cols = 3)

multiplot(plotlist = pplots2, cols = 3)

print(
  ggplot(all_MiHA_SNP, aes(x = MiHAs, ymax = frequency, ymin = 0, color = Group)) +
    geom_linerange(size = 3) +
    geom_hline(yintercept = 0) +
    ggtitle("Known MiHAs") +
    theme_bw() +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
)

