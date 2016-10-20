load("../Data/ID_table.RData")

known_miha_freq <- read.csv("../WW_MiHA/HLA_restricted_knownMiHAs.csv", header = F)
colnames(known_miha_freq) <- c("GroupType", "GroupID", "CHROM", "REF", "ALT", "HLA_type", "SNP")

table(known_miha_freq[c("HLA_type","SNP")])

known_miha_freq$HLA_SNP <- sapply(1:dim(known_miha_freq)[1], 
                                  function(x) paste0(known_miha_freq$HLA_type[x], "-", known_miha_freq$SNP[x]))

single_MiHA_stats <- as.data.frame(table(known_miha_freq[c("GroupType", "HLA_SNP")]))

GroupID_MiHA_stats <- as.data.frame(table(known_miha_freq[c("GroupID", "HLA_SNP")]))
GroupID_MiHA_stats <- GroupID_MiHA_stats[GroupID_MiHA_stats$Freq>0, ]

length(unique(single_MiHA_stats$HLA_SNP))
length(unique(GroupID_MiHA_stats$HLA_SNP))
length(unique(known_miha_freq$SNP)) # 37 MiHA SNPs

unique_IDs <- unique(known_miha_freq$GroupID) # 131 groups

SNP_mat <- matrix(data = 0, nrow = length(unique_IDs), ncol = length(unique(GroupID_MiHA_stats$HLA_SNP)) + 2)
colnames(SNP_mat) <- c("GroupID", "GroupType", as.character(unique(GroupID_MiHA_stats$HLA_SNP)))
SNP_mat <- as.data.frame(SNP_mat)

for(id in 1:length(unique_IDs)){
  
  selected_group <- known_miha_freq[known_miha_freq$GroupID == unique_IDs[id],]
  SNP_mat$GroupID[id] <- unique_IDs[id]
  SNP_mat$GroupType[id] <- ID_table$Group[which(ID_table$GroupID == unique_IDs[id])[1]]
  
  if(length(unique(selected_group$HLA_SNP)) == length(selected_group$HLA_SNP)){ # no duplicated SNPs for one group
    
    SNP_mat[id, selected_group$HLA_SNP] <- 1
    
  }else{ # one homozygous SNP (duplicated)
    
    SNP_counts <- as.data.frame(table(selected_group$HLA_SNP), stringsAsFactors = F)
    SNP_mat[id, SNP_counts[, 1]] <- SNP_counts[, 2]
    
  }
 
}
write.csv(SNP_mat, file = "../WW_MiHA/SNP_mat.csv", row.names = FALSE)

sort(colSums(SNP_mat[, -c(1,2)]), decreasing = T)
sort(colSums(SNP_mat[which(SNP_mat$GroupType == "n"), -c(1,2)]), decreasing = T)
sort(colSums(SNP_mat[which(SNP_mat$GroupType == "a"), -c(1,2)]), decreasing = T)

# presence-absence 
presence_SNP_mat <- SNP_mat
presence_SNP_mat[presence_SNP_mat[, -c(1,2)]>0] <- 1




# A*02:01
A_SNP <- SNP_mat[, c(1:2, 4:13)]
A_SNP <- A_SNP[-which(rowSums(A_SNP[, -c(1,2)])==0), ]
A_a_SNP <- A_SNP[A_SNP$GroupType == "a", ]
A_n_SNP <- A_SNP[A_SNP$GroupType == "n", ]

library(corrgram)
corrgram(A_a_SNP[, -1], order = T, lower.panel = panel.shade,
         upper.panel = panel.pie)

corrgram(A_a_SNP[, -c(1,2)], order = T, lower.panel = panel.ellipse,text.panel=panel.txt,
         upper.panel = panel.pts, diag.panel = panel.minmax)

col.corrgram <- function(ncol){   
  colorRampPalette(c("darkgoldenrod4", "burlywood1",
                     "darkkhaki", "darkgreen"))(ncol)}
nonZeros_rows <- which(rowSums(A_a_SNP[, -c(1,2)]) !=0)
corrgram(A_a_SNP[, -c(1,2)], order=F, lower.panel=panel.pts, 
         upper.panel=panel.cor, diag.panel = panel.density, text.panel=panel.txt, 
         main="Correlogram of HLA-A*02:01 restricted MiHAs \n in aGVHD groups", 
         cor.method = "spearman") #,
         # col.regions = colorRampPalette(c("darkgoldenrod4", "burlywood1",
         #                                  "darkkhaki", "darkgreen")))


corrgram(A_n_SNP[, -c(1,2)], order = F, lower.panel = panel.shade,
         upper.panel = panel.pie, text.panel=panel.txt, 
         main="Correlogram of HLA-A*02:01 restricted MiHAs \n in non-aGVHD groups") 


###
# B*07:02
B_SNP <- SNP_mat[, c(1:2, 17:26)]
B_SNP <- B_SNP[-which(rowSums(B_SNP[, -c(1,2)])==0), ]
B_a_SNP <- B_SNP[B_SNP$GroupType == "a", ]
B_n_SNP <- B_SNP[B_SNP$GroupType == "n", ]

corrgram(B_a_SNP[, -c(1,2)], order= F, lower.panel=panel.shade, 
         upper.panel=panel.pie, text.panel=panel.txt, 
         main="Correlogram of HLA-B*07:02 restricted MiHAs \n in aGVHD groups") #,
# col.regions = colorRampPalette(c("darkgoldenrod4", "burlywood1",
#                                  "darkkhaki", "darkgreen")))


corrgram(B_n_SNP[, -c(1,2)], order = F, lower.panel = panel.shade,
         upper.panel = panel.pie, text.panel=panel.txt, 
         main="Correlogram of HLA-B*07:02 restricted MiHAs \n in non-GVHD groups") 