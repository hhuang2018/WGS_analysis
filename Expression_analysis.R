expression_table <- read.table(file = "../ExpressionAnalysis/all_expression.txt", 
                               header = T, sep = "\t")
colnames(expression_table)[1:2] <- c("ENSG.Number", "Gene")

# affExp_table <- read.csv("../ExpressionAnalysis/AffExp.csv", header = T)
# rm(affExp_table)

num_tissues <- dim(expression_table)[2] - 2
num_genes <- dim(expression_table)[1]
# expression * 10 >=0 --> log10()
converted_exp_table <- expression_table[, -c(1,2)] * 10
converted_exp_table[converted_exp_table<1] <- 1
converted_exp_table <- log10(converted_exp_table)

row_max <- sapply(1:num_genes, function(x) max(converted_exp_table[x, ]))
rm_id <- which(row_max == 0)
converted_exp_table <- converted_exp_table[-rm_id, ]
ENSG_names <- expression_table$ENSG.Number[-rm_id]
Gene_names <- expression_table$Gene[-rm_id]
num_genes <- dim(expression_table)[1] - length(rm_id)

x_hat_table <- t(sapply(1:num_genes, 
                      function(x) converted_exp_table[x, ]/max(converted_exp_table[x, ])))

Tau_list <- data.frame(EnsembleID = as.character(num_genes),
                       Gene = character(num_genes),
                       Tau = numeric(num_genes),
                       HighExpTissue = character(num_genes),
                       stringsAsFactors = F)

for(id in 1:num_genes){
  
  Tau_list$Tau[id] <- sum(1- unlist(x_hat_table[id, ])) / (num_tissues-1) 
  Tau_list$EnsembleID[id] <- as.character(ENSG_names[id])
  Tau_list$Gene[id] <- as.character(Gene_names[id])
  Tau_list$HighExpTissue[id] <- names(which.max(x_hat_table[id,]))
  
}

# high_spec_ids <- which(Tau_list$Tau == 1)
# names(which.max(x_hat_table[high_spec_ids[1],]))

Tau_list <- Tau_list[order(Tau_list$Tau, decreasing = T), ]
Tau_list <- within(Tau_list, EnsembleID <- factor(EnsembleID, levels=factor(Tau_list$EnsembleID)))

write.csv(Tau_list, file = "../ExpressionAnalysis/Tau_list.csv", row.names = F)

length(which(Tau_list$Tau==1))

library(ggplot2)
ggplot(Tau_list[1:100, ], aes(x=Gene, y=Tau)) +
  geom_bar(stat="identity", position=position_dodge(), fill = "#56B4E9") +
  ggtitle("Tau") +
  # ggtitle(Typing_name) +
  # xlab(paste0(Typing_name, " typing")) + ylab("Count of groups")  + 
  xlab(Typing_name) + ylab("Tau")  + 
  theme(axis.text.x = element_text(angle = Angle, hjust = H))