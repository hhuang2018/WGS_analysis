allele_affinity_table <- read.csv("../ExpressionAnalysis/AllBInderFinal_052017/F-M_6_alleles_affinites.csv", stringsAsFactors = F)

num_groups <- dim(allele_affinity_table)[1]

binder_table <- data.frame(PID = numeric(num_groups),
                           Outcome = character(num_groups),
                           A.aff = numeric(num_groups),
                           B.aff = numeric(num_groups),
                           C.aff = numeric(num_groups),
                           all.aff = numeric(num_groups),
                           stringsAsFactors = F)
for(id in 1:num_groups){
  
  binder_table$PID[id] <- allele_affinity_table$PID[id]
  binder_table$Outcome[id] <- allele_affinity_table$Outcome[id]
  binder_table$A.aff[id] <- min(allele_affinity_table[id, c(13, 15)])
  binder_table$B.aff[id] <- min(allele_affinity_table[id, c(17, 19)])
  binder_table$C.aff[id] <- min(allele_affinity_table[id, c(21, 23)])
  binder_table$all.aff[id] <- min(binder_table[id, 3:5])
  
}
# aa <- paste0("HLA-A > ",binder_table$Outcome, "GVHD")

mc_length <- data.frame(Affinity = c(binder_table$A.aff, binder_table$B.aff, binder_table$C.aff, binder_table$all.aff),
                        Group = c(paste0("HLA-A > ", binder_table$Outcome, "GVHD"), 
                                  paste0("HLA-B > ", binder_table$Outcome, "GVHD"),
                                  paste0("HLA-C > ", binder_table$Outcome, "GVHD"),
                                  paste0("ALL Allele > ", binder_table$Outcome, "GVHD")),
                        Group2 <- c(binder_table$Outcome, 
                                    binder_table$Outcome,
                                    binder_table$Outcome,
                                    binder_table$Outcome),
                        # rep("Matched WGS", length(all_matched_Percent$IBDPercent)), rep("Random WGS", length(all_random_Percent$IBDPercent))),
                        stringsAsFactors = F)

mc_length$Group <- factor(mc_length$Group, 
                          levels = c("HLA-A > aGVHD", 
                                     "HLA-A > nGVHD",
                                     "HLA-B > aGVHD",
                                     "HLA-B > nGVHD",
                                     "HLA-C > aGVHD",
                                     "HLA-C > nGVHD",
                                     "ALL Allele > aGVHD",
                                     "ALL Allele > nGVHD"))

ggplot(mc_length, aes(x = Group , y = Affinity, fill = Group2)) + 
  geom_boxplot()+ #outlier.colour = "gray", outlier.size = 0) +
  guides(fill=FALSE)  +
  scale_y_continuous("Affinity Score")
