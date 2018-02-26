library(corrgram)
load("../Data/ID_table.RData")

# known_miha_freq <- read.csv("../WW_MiHA/HLA_restricted_knownMiHAs.csv", header = F)
known_miha_freq <- read.delim("../WW_MiHA/Restricted_known_MiHAs.txt", header = F)
colnames(known_miha_freq) <- c("GroupType", "GroupID",  "HLA_type", "SNP", "CHROM", "REF", "ALT")


# table(known_miha_freq[c("HLA_type","SNP")])

known_miha_freq$HLA_SNP <- sapply(1:dim(known_miha_freq)[1], 
                                  function(x) paste0(known_miha_freq$HLA_type[x], "-", known_miha_freq$SNP[x]))

#################
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
presence_SNP_mat <- SNP_mat[, -c(1,2)]
presence_SNP_mat[presence_SNP_mat>0] <- 1
presence_SNP_mat <- cbind(SNP_mat[, c(1,2)], presence_SNP_mat)

nonZeros_rows <- which(rowSums(presence_SNP_mat[, -c(1,2)]) !=0)
corrgram(presence_SNP_mat[nonZeros_rows, -c(1,2)], order=F, lower.panel=panel.shade, 
         upper.panel=panel.pie, diag.panel = NULL, text.panel=panel.txt, 
         main="Correlogram of all MiHAs", 
         cor.method = "spearman")

# A*02:01
A_SNP <- SNP_mat[, c(1:2, 4:13)]
A_SNP <- A_SNP[-which(rowSums(A_SNP[, -c(1,2)])==0), ]
A_a_SNP <- A_SNP[A_SNP$GroupType == "a", ]
A_n_SNP <- A_SNP[A_SNP$GroupType == "n", ]

# library(corrgram)
corrgram(A_a_SNP[, -1], order = T, lower.panel = panel.shade,
         upper.panel = panel.pie)

# corrgram(A_a_SNP[, -c(1,2)], order = T, lower.panel = panel.ellipse,text.panel=panel.txt,
#          upper.panel = panel.pts, diag.panel = panel.minmax)

col.corrgram <- function(ncol){   
  colorRampPalette(c("darkgoldenrod4", "burlywood1",
                     "darkkhaki", "darkgreen"))(ncol)}
nonZeros_rows <- which(rowSums(A_a_SNP[, -c(1,2)]) !=0)
new_A_a_SNP <- A_a_SNP[,-c(1,2)]
new_A_a_SNP[new_A_a_SNP>1] <- 1

new_A_a_SNP <- cbind(A_a_SNP[, c(1,2)], new_A_a_SNP)
Labels_SNP <- gsub("A\\*02:01-", "",colnames(new_A_a_SNP[nonZeros_rows, -c(1,2)]))
corrgram(new_A_a_SNP[nonZeros_rows, -c(1,2)], order=F, lower.panel=panel.shade, 
         upper.panel=panel.pie, diag.panel = NULL, text.panel=panel.txt, 
         labels = Labels_SNP, #label.pos = c(0,0), label.srt = 45,
         main="HLA-A*02:01 restricted MiHAs \n in aGVHD groups", 
         cor.method = "spearman")

corrgram(A_a_SNP[nonZeros_rows, -c(1,2)], order=F, lower.panel=panel.shade, 
         upper.panel=panel.pie, diag.panel = NULL, text.panel=panel.txt, 
         main="HLA-A*02:01 restricted MiHAs \n in aGVHD groups", 
         cor.method = "spearman") #,
         # col.regions = colorRampPalette(c("darkgoldenrod4", "burlywood1",
         #                                  "darkkhaki", "darkgreen")))

nonZeros_rows <- which(rowSums(A_n_SNP[, -c(1,2)]) !=0)
p_nGVHD <- corrgram(A_n_SNP[nonZeros_rows, -c(1,2)], order = F, lower.panel = panel.shade,
         upper.panel = panel.pie, text.panel=panel.txt, 
         labels = Labels_SNP,
         main="Correlogram of HLA-A*02:01 restricted MiHAs \n in non-aGVHD groups",
         cor.method = "spearman") 

new_A_n_SNP <- A_n_SNP[,-c(1,2)]
new_A_n_SNP[new_A_n_SNP>1] <- 1

new_A_n_SNP <- cbind(A_n_SNP[, c(1,2)], new_A_n_SNP)
Labels_SNP <- gsub("A\\*02:01-", "",colnames(new_A_n_SNP[nonZeros_rows, -c(1,2)]))
corrgram(new_A_n_SNP[nonZeros_rows, -c(1,2)], order=F, lower.panel=panel.shade, 
         upper.panel=panel.pie, diag.panel = NULL, text.panel=panel.txt, 
         labels = Labels_SNP,
         main="HLA-A*02:01 restricted MiHAs \n in non-aGVHD groups", 
         cor.method = "spearman")

##### Log Ratio
known_MiHA_coordinates <- read.table("../WW_MiHA/Known_MiHA_coordinates.txt", header = T, stringsAsFactors = F)
num_rows <- dim(known_miha_freq)[1] 
known_miha_freq$MiHA_HLA_SNP_Gene <- sapply(1:num_rows, function(x) {ind <- which(known_MiHA_coordinates$SNPs == known_miha_freq$SNP[x])
paste0(known_MiHA_coordinates$MiHAs[ind], "<>", known_miha_freq$HLA_SNP[x], "<>", known_MiHA_coordinates$Gene[ind])}) 

unique_MiHAs <- unique(known_miha_freq$MiHA_HLA_SNP_Gene)
num_unique_MiHAs <- length(unique_MiHAs)

MiHA_LLR <- data.frame(MiHA_HLA_SNP_Gene = character(num_unique_MiHAs),
                        aGVHD_count = numeric(num_unique_MiHAs),
                        nGVHD_count = numeric(num_unique_MiHAs),
                        LLR = numeric(num_unique_MiHAs),
                       stringsAsFactors = F)

for(id in 1:num_unique_MiHAs){
  
  index <- which(known_miha_freq$MiHA_HLA_SNP_Gene %in% unique_MiHAs[id])
  
  count_summary <- as.data.frame(table(known_miha_freq$GroupType[index]))
  
  MiHA_LLR$MiHA_HLA_SNP_Gene[id] <- unique_MiHAs[id]
  MiHA_LLR$aGVHD_count[id] <- count_summary$Freq[which(count_summary$Var1[index] == "a")]
  MiHA_LLR$nGVHD_count[id] <- count_summary$Freq[which(count_summary$Var1[index] == "n")]
  MiHA_LLR$LLR[id] <- log10(MiHA_LLR$aGVHD_count[id] / MiHA_LLR$nGVHD_count[id])
  
}




# multiplot(p_aGVHD, p_nGVHD, cols = 2)
###
# B*07:02
# B_SNP <- SNP_mat[, c(1:2, 17:26)]
B_SNP <- SNP_mat[, c(1:2, 18:27)]
B_SNP <- B_SNP[-which(rowSums(B_SNP[, -c(1,2)])==0), ]
B_a_SNP <- B_SNP[B_SNP$GroupType == "a", ]
B_n_SNP <- B_SNP[B_SNP$GroupType == "n", ]

nonZeros_rows <- which(rowSums(B_a_SNP[, -c(1,2)]) !=0)
new_B_a_SNP <- B_a_SNP[,-c(1,2)]
new_B_a_SNP[new_B_a_SNP>1] <- 1
new_B_a_SNP <- cbind(B_a_SNP[, c(1,2)], new_B_a_SNP)

Labels_SNP <- gsub("B\\*07:02-", "",colnames(B_a_SNP[nonZeros_rows, -c(1,2)]))
corrgram(B_a_SNP[nonZeros_rows, -c(1,2)], order= F, lower.panel=panel.shade, 
         upper.panel=panel.pie, text.panel=panel.txt, 
         labels = Labels_SNP,
         main="HLA-B*07:02 restricted MiHAs \n in aGVHD groups",
         cor.method = "spearman") #,
# col.regions = colorRampPalette(c("darkgoldenrod4", "burlywood1",
#                                  "darkkhaki", "darkgreen")))
corrgram(new_B_a_SNP[nonZeros_rows, -c(1,2)], order= F, lower.panel=panel.shade, 
         upper.panel=panel.pie, text.panel=panel.txt, 
         labels = Labels_SNP,
         main="Correlogram of HLA-B*07:02 restricted MiHAs \n in aGVHD groups",
         cor.method = "spearman")

new_B_n_SNP <- B_n_SNP[,-c(1,2)]
new_B_n_SNP[new_B_n_SNP>1] <- 1
new_B_n_SNP <- cbind(B_n_SNP[, c(1,2)], new_B_n_SNP)

Labels_SNP <- gsub("B\\*07:02-", "",colnames(new_B_n_SNP[, -c(1,2)]))
corrgram(new_B_n_SNP[, -c(1,2)], order = F, lower.panel = panel.shade,
         upper.panel = panel.pie, text.panel=panel.txt, 
         labels = Labels_SNP,
         main="HLA-B*07:02 restricted MiHAs \n in non-GVHD groups",
         cor.method = "spearman") 

corrgram(B_n_SNP[, -c(1,2)], order = F, lower.panel = panel.shade,
         upper.panel = panel.pie, text.panel=panel.txt, 
         labels = Labels_SNP,
         main="Correlogram of HLA-B*07:02 restricted MiHAs \n in non-GVHD groups",
         cor.method = "spearman") 


##########################
# co-occurence matrix
##########################
SNP_mat <- read.csv(file = "../WW_MiHA/SNP_mat.csv")

# presence-absence 
presence_SNP_mat <- SNP_mat[, -c(1,2)]
presence_SNP_mat[presence_SNP_mat>0] <- 1
# presence_SNP_mat <- cbind(SNP_mat[, c(1,2)], presence_SNP_mat)
presence_SNP_mat <- as.matrix(presence_SNP_mat)
rownames(presence_SNP_mat) <- SNP_mat[, 1]

SNP_adjacent_mat <- presence_SNP_mat %*% t(presence_SNP_mat)  ## adjacency matrix

library(igraph)
SNP_graph <- graph.adjacency(SNP_adjacent_mat, weighted=T, mode = "undirected") # build a graph
SNP_graph <- simplify(SNP_graph) # remove loops
V(SNP_graph)$label <- V(SNP_graph)$name  # set labels of vertices
V(SNP_graph)$degree <- degree(SNP_graph) # set degress of vertices
# set seed to make the layout reproducible
set.seed(3952)
layout1 <- layout.fruchterman.reingold(SNP_graph)
plot(SNP_graph, layout=layout1)
plot(SNP_graph, layout=layout.kamada.kawai)
tkplot(SNP_graph, layout=layout.kamada.kawai)

# 
V(SNP_graph)$label.cex <- 2.2 * V(SNP_graph)$degree / max(V(SNP_graph)$degree)+ .2
V(SNP_graph)$label.color <- rgb(0, 0, .2, .8)
V(SNP_graph)$frame.color <- NA
egam <- (log(E(SNP_graph)$weight)+.4) / max(log(E(SNP_graph)$weight)+.4)
E(SNP_graph)$color <- rgb(.5, .5, 0, egam)
E(SNP_graph)$width <- egam
# plot the graph in layout1
plot(SNP_graph, layout=layout1)