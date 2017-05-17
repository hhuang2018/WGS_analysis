all_binders_table <- read.csv("../ExpressionAnalysis/aws_MiHAIP-v1.4.5_analysis/All_BinderOut_trimed_tau_mata.txt", header = T)

colnames(all_binders_table) <- c("GroupID", "HetHom", "CHROM", "POS", 
                                 "DonorNT", "DonorNTfreq", "RecipientNT", "RecipientNTfreq",
                                 "ENST", "HLATyping", "peptide", "AffinityScore",
                                 "normRank", "ENST.1", "ENSG", "Gene")

MAF_binder_table <- read.csv("../ExpressionAnalysis/binder/BinderMAF_Joined_rmSingle.txt")

unique(MAF_binder_table$HighExpTissue)
table(MAF_binder_table$HighExpTissue)

Blood_related_tissues <- c("Whole.Blood", "Cells...EBV.transformed.lymphocytes")
Male_specific_tissues <- c("Testis", "Prostate")         # male specific
Female_specific_tissues <- c("Ovary", "Fallopian.Tube")  # Female specific

MAF_binder_table$TissueSpec <- sapply(1:dim(MAF_binder_table)[1], 
                                      function(x) if(is.element(MAF_binder_table$HighExpTissue[x], Blood_related_tissues)){
                                        "Blood_related"} else if(is.element(MAF_binder_table$HighExpTissue[x], Male_specific_tissues)){
                                          "Male_specific"} else if(is.element(MAF_binder_table$HighExpTissue[x], Female_specific_tissues)){
                                            "Female_specific"}else{
                                              "General"})

predicted_peptides_all <- sapply(1:dim(MAF_binder_table)[1], function(x) paste0(MAF_binder_table$Pep[x], "-", MAF_binder_table$Gene[x]))
predicted_peptides <- unique(predicted_peptides_all)
num_peptides <- length(predicted_peptides)
peptide_summary_table <- data.frame(peptide = predicted_peptides, 
                                    nGroup = numeric(num_peptides),
                                    aGroup = numeric(num_peptides),
                                    Tissues = character(num_peptides),
                                    stringsAsFactors = F)
for(id in 1:num_peptides){
  
  pp_ids <- which(predicted_peptides_all %in% predicted_peptides[id])
  
  stats <- as.data.frame(table(MAF_binder_table$GVHD[pp_ids]))
  if(stats$Var1[1] == "a"){
    peptide_summary_table$aGroup[id] <- stats$Freq[1]
    if(dim(stats)[1] == 2){
      peptide_summary_table$nGroup[id] <- stats$Freq[2]
    }
  }else{
    peptide_summary_table$nGroup[id] <- stats$Freq[1]
    if(dim(stats)[1] == 2){
      peptide_summary_table$aGroup[id] <- stats$Freq[2]
    }
  }
  
  Tissue_spec <- unique(MAF_binder_table$TissueSpec[pp_ids])
  if(length(Tissue_spec) == 1) peptide_summary_table$Tissues[id] <- Tissue_spec else print(id)
  
}

library(ggplot2)
ggplot(peptide_summary_table, aes(x=nGroup, y=aGroup)) + geom_point(aes(color = factor(Tissues)))


####################### 04/06/2017
expression_table <- read.csv(file = "../ExpressionAnalysis/aws_MiHAIP-v1.4.5_analysis/New_all_binder_joined.csv", 
                             header = T, sep = "\t")
cleanup_table <- read.table(file = "../ExpressionAnalysis/aws_MiHAIP-v1.4.5_analysis/New_all_binder_cleanup.txt", 
                          header = T, sep = "\t")

Blood_related_tissues <- c("Whole.Blood", "Cells...EBV.transformed.lymphocytes")
Male_specific_tissues <- c("Testis", "Prostate")         # male specific
Female_specific_tissues <- c("Ovary", "Fallopian.Tube")  # Female specific

expression_table$TissueSpec <- sapply(1:dim(expression_table)[1], 
                                      function(x) if(is.element(expression_table$HighExpTissue[x], Blood_related_tissues)){
                                        "Blood_related"} else if(is.element(expression_table$HighExpTissue[x], Male_specific_tissues)){
                                          "Male_specific"} else if(is.element(expression_table$HighExpTissue[x], Female_specific_tissues)){
                                            "Female_specific"}else{
                                              "General"})

predicted_peptides_all <- sapply(1:dim(cleanup_table)[1], function(x) paste0(cleanup_table$Pep[x], "-", cleanup_table$Gene_Symbol[x]))
predicted_peptides <- unique(predicted_peptides_all)
num_peptides <- length(predicted_peptides)
peptide_summary_table <- data.frame(peptide = predicted_peptides, 
                                    nGroup = numeric(num_peptides),
                                    aGroup = numeric(num_peptides),
                                    Tissues = character(num_peptides),
                                    stringsAsFactors = F)
for(id in 1:num_peptides){
  
  pp_ids <- which(predicted_peptides_all %in% predicted_peptides[id])
  
  stats <- as.data.frame(table(cleanup_table$Outcome[pp_ids]))
  if(stats$Var1[1] == "a"){
    peptide_summary_table$aGroup[id] <- stats$Freq[1]
    if(dim(stats)[1] == 2){
      peptide_summary_table$nGroup[id] <- stats$Freq[2]
    }
  }else{
    peptide_summary_table$nGroup[id] <- stats$Freq[1]
    if(dim(stats)[1] == 2){
      peptide_summary_table$aGroup[id] <- stats$Freq[2]
    }
  }
  
  # Tissue_spec <- unique(MAF_binder_table$TissueSpec[pp_ids])
  # if(length(Tissue_spec) == 1) peptide_summary_table$Tissues[id] <- Tissue_spec else print(id)
  
}

library(ggplot2)

ggplot(peptide_summary_table, aes(x = nGroup, y = aGroup)) +
  geom_point(shape = 19, alpha = 1/4) +
  # geom_smooth(method = lm)
  geom_abline(slope = 1, intercept = 0, colour = "#0072B2")


peptide_summary_table$LOD <- log10((peptide_summary_table$aGroup + 1e-8) / (peptide_summary_table$nGroup + 1e-8))
ordered_peptide_summary <-  peptide_summary_table[order(peptide_summary_table$LOD, decreasing = T), ] 
ordered_peptide_summary$group <- sapply(1:dim(ordered_peptide_summary)[1], 
                                        function(x) if(ordered_peptide_summary$LOD[x] > 5){"aGVHD-dominant"
                                        }else if(ordered_peptide_summary$LOD[x] < -5){"nGVHD-dominant"
                                        }else {"Normal"})
ordered_peptide_summary$peptide <- factor(ordered_peptide_summary$peptide, levels = ordered_peptide_summary$peptide)
ggplot(ordered_peptide_summary, aes(x = peptide, y = LOD, color = factor(group))) + geom_point(size = 1, alpha = 0.5) +
  theme(axis.text.x = element_blank())

ggplot(ordered_peptide_summary[which(ordered_peptide_summary$group == "Normal"),], aes(x = peptide, y=LOD)) +
  geom_point(size = 1, alpha = 0.5, colour = "steelblue") +
  theme(axis.text.x = element_blank())


#####
new_expression_table <- read.csv(file = "../ExpressionAnalysis/aws_MiHAIP-v1.4.5_analysis/New_all_binder_joined_rearranged.txt", 
                             header = T, sep = "\t")
# new_expression_table <- read.csv("/mnt/cloudbiodata_nfs_1/hli_scratch/wwang/AWS_MiHAIP-1.4.5_output/New_all_binder_joined_rearranged.txt", header = T, sep = "\t")
Blood_related_tissues <- c("Whole.Blood", "Cells...EBV.transformed.lymphocytes")
Male_specific_tissues <- c("Testis", "Prostate")         # male specific
Female_specific_tissues <- c("Ovary", "Fallopian.Tube")  # Female specific

new_expression_table$TissueSpec <- sapply(1:dim(new_expression_table)[1], 
                                      function(x) if(is.element(new_expression_table$HighExpTissue[x], Blood_related_tissues)){
                                        "Blood_related"} else if(is.element(new_expression_table$HighExpTissue[x], Male_specific_tissues)){
                                          "Male_specific"} else if(is.element(new_expression_table$HighExpTissue[x], Female_specific_tissues)){
                                            "Female_specific"}else{
                                              "General"})

predicted_peptides_all <- sapply(1:dim(new_expression_table)[1], function(x) paste0(new_expression_table$Pep[x], "-", new_expression_table$Gene_Symbol[x]))
predicted_peptides <- unique(predicted_peptides_all)
num_peptides <- length(predicted_peptides)
peptide_summary_table <- data.frame(peptide = predicted_peptides, 
                                    nGroup = numeric(num_peptides),
                                    aGroup = numeric(num_peptides),
                                    Tissues = character(num_peptides),
                                    stringsAsFactors = F)
for(id in 1:num_peptides){
  
  pp_ids <- which(predicted_peptides_all %in% predicted_peptides[id])
  
  stats <- as.data.frame(table(new_expression_table$Outcome[pp_ids]))
  if(stats$Var1[1] == "a"){
    peptide_summary_table$aGroup[id] <- stats$Freq[1]
    if(dim(stats)[1] == 2){
      peptide_summary_table$nGroup[id] <- stats$Freq[2]
    }
  }else{
    peptide_summary_table$nGroup[id] <- stats$Freq[1]
    if(dim(stats)[1] == 2){
      peptide_summary_table$aGroup[id] <- stats$Freq[2]
    }
  }
  
  Tissue_spec <- unique(MAF_binder_table$TissueSpec[pp_ids])
  if(length(Tissue_spec) == 1) peptide_summary_table$Tissues[id] <- Tissue_spec else print(id)
  
}

