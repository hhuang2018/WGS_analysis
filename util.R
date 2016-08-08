#### Get gene names within a given position range
get.GeneNames <- function(chrom, startPos, endPos, GeneList){

  chrom_list <- GeneList[which(GeneList$chrom %in% chrom), ]
  
  if(dim(chrom_list)[1] > 0){
    
    regionStart <- intersect(which(chrom_list$txStart <= startPos), which(chrom_list$txEnd >= startPos))
    if(length(regionStart) == 0) regionStart <- which(chrom_list$txEnd >= startPos)[1]
    if(is.na(regionStart)) regionStart <- startPos
        
    regionEnd <- intersect(which(chrom_list$txEnd >= endPos), which(chrom_list$txStart <= endPos))
    if(length(regionEnd) == 0) regionEnd <- which(chrom_list$txStart <= endPos)[1]
    if(is.na(regionEnd)) regionEnd <- 0
    
    geneRange <- c(min(regionStart),max(regionEnd))
    
    if(geneRange[1] <= geneRange[2]){
      
      geneNames <- as.character(chrom_list$GeneName[geneRange[1]:geneRange[2]])
      geneNames[which(chrom_list$strand[geneRange[1]:geneRange[2]] %in% "-")] <- paste0(geneNames[which(chrom_list$strand[geneRange[1]:geneRange[2]] %in% "-")], "(-)")
      
      GeneNames <- paste0(unique(geneNames), collapse = ",")
      
    } else GeneNames <- "OutOfTxRange"
    
  } else GeneNames <- chrom
  
  return(GeneNames)
}

#####################
##
#####################
is.same.gt <- function(vcf_gt){
  
  gt_1 <- sort(as.numeric(unlist(strsplit(vcf_gt[1], "/"))))
  gt_2 <- sort(as.numeric(unlist(strsplit(vcf_gt[2], "/"))))
  
  uniq_gt_1 <- unique(gt_1)
  uniq_gt_2 <- unique(gt_2)
  
  if(length(uniq_gt_1) == 1){ # if the first one is homozygous
    
    if(length(uniq_gt_2) == 1){
      # both are homozygous
#       if(uniq_gt_1 == uniq_gt_2){
#         # both alleles are the same
#         genotype_table$GT[rind] <- 0
#       }else genotype_table$GT[rind] <- 2 # both alleles are different
      GT <- 2
    }else{
      # if the second one is heterozygous
      if(length(unique(c(uniq_gt_1, uniq_gt_2))) == 2 ){
        # one allele is the same
        GT <- 1
      }else GT <- 2 # both alleles are different
      
    }
  }else{ # if the first one is heterozygous #####
    
    if(length(uniq_gt_2) == 1){
      # if the second one is homozygous
      if(length(unique(c(uniq_gt_1, uniq_gt_2))) == 2 ){
        # one allele is the same
        GT <- 1
      }else GT <- 2 # both alleles are different
      
    }else{
      # if both are heterozygous
      # shared_allele_num <- length(intersect(uniq_gt_1, uniq_gt_2)) 
      switch(as.character(length(intersect(uniq_gt_1, uniq_gt_2)) ),
             "0" = GT <- 2, # both alleles are different
             "1" = GT <- 1, # one allele is different
             "2" = GT <- 0) # both alleles are the same 
    }
  }
  return(GT)
}

############
gene_variant_stats <- function(geneInfo, geneVariants, promoterRegion = 1000){
  #
  # returns the variant statistics by functional regions (introns, exons)
  #
  
  startIndex <- as.integer(unlist(strsplit(as.character(geneInfo$exonStarts), ",")))
  endIndex <- as.integer(unlist(strsplit(as.character(geneInfo$exonEnds), ",")))
  
  gene_length <- geneInfo$txEnd - geneInfo$txStart + 1
  exon_length <- sum(endIndex - startIndex) + geneInfo$exonCount
  if(geneInfo$exonCount > 1){
    intron_length <- sum(startIndex[2:geneInfo$exonCount] - endIndex[1:(geneInfo$exonCount-1)]) - (geneInfo$exonCount - 1) # only intron region, without UTR
  } else intron_length <- 0
  # intron_length <- sum(startIndex[2:geneInfo$exonCount] - endIndex[1:(geneInfo$exonCount-1)]) - (geneInfo$exonCount - 1) # only intron region, without UTR
  utr_length <- gene_length - exon_length - intron_length
  
  RegionIndex <- vector(mode = "numeric", length = 2*geneInfo$exonCount)
  RegionIndex[seq(from = 1, to = 2*geneInfo$exonCount, by = 2)] <- startIndex
  RegionIndex[seq(from = 2, to = 2*geneInfo$exonCount, by = 2)] <- endIndex
  
  total_variants <- geneVariants[(geneVariants$POS >= geneInfo$txStart) & (geneVariants$POS <= geneInfo$txEnd), ]
  total_variants <- total_variants[which(total_variants$FILTER == "PASS"), ]
  total_num <- dim(total_variants)[1]
  
  RegionNumber <- findInterval(total_variants$POS, RegionIndex)
  total_variants$ID <- sapply(RegionNumber, function(x) RegionNumber2RegionName(x, length(RegionIndex)))
  
  promoter_variants <- geneVariants[(geneVariants$POS >= (geneInfo$txStart - promoterRegion)) & (geneVariants$POS <= (geneInfo$txStart-1)), ]
  promoter_variants <- promoter_variants[which(promoter_variants$FILTER == "PASS"), ]
  promoter_num <- dim(promoter_variants)[1]
  
  utr_intron_variant_index <- which(RegionNumber %% 2 == 0)
  utr_variant_index <- which(RegionNumber == 0 | RegionNumber == geneInfo$exonCount*2)
  if(length(utr_variant_index) != 0){ 
    intron_variant_index <- utr_intron_variant_index[-which(utr_intron_variant_index %in% utr_variant_index)]
  }else intron_variant_index <- utr_intron_variant_index
  exon_variant_index <- which(RegionNumber %% 2 == 1)
  insert_variant_index <- which(nchar(total_variants$REF) < nchar(total_variants$ALT))
  deletion_variant_index <- which(nchar(total_variants$REF) > nchar(total_variants$ALT))
  
  intron_variants <- total_variants[intron_variant_index, ]
  exon_variants <- total_variants[exon_variant_index, ]
  utr_variants <- total_variants[utr_variant_index, ]
  insert_variants <- total_variants[insert_variant_index, ]
  deletion_variants <- total_variants[deletion_variant_index, ]
  
  return(list(promoter_variants = promoter_variants,
              all_variants = total_variants, 
              intron_variants = intron_variants,
              exon_variants = exon_variants,
              utr_variants = utr_variants,
              insert_variants = insert_variants,
              deletion_variants = deletion_variants,
              gene_length = gene_length,
              exon_length = exon_length,
              intron_length = intron_length,
              utr_length = utr_length,
              promoter_length = promoterRegion,
              variants_stats = data.frame(total_num = total_num,
                                          intron_variants_num = length(intron_variant_index),
                                          exon_variants_num = length(exon_variant_index), 
                                          utr_variants_num = length(utr_variant_index),
                                          promoter_variants_num = promoter_num, 
                                          insert_variants_num = length(insert_variant_index),
                                          deletion_variants_num = length(deletion_variant_index))
              ))
}

############
RegionNumber2RegionName <- function(RegionNumber, num_reigons){
  # regionIndex <- findInterval(Site, startIndex)
  if(RegionNumber == 0){
    regionName <- "fiveUTR" #, Site - startIndex[1] + 1)
  }else if(RegionNumber == num_reigons){
    regionName <- "threeUTR"
  } else{
    if(RegionNumber %% 2 == 1){ # Exon
      regionName <- paste0("Exon", (RegionNumber+1)/2)
    }else { # Intron
      regionName <- paste0("Intron", RegionNumber/2)
    }
  }
  return(regionName)
}

####################
# Compare Variant 1 (Donor) and Variant 2 (Recipient) - stats of the differences
####################
compare_variants <- function(Variant1, Variant2, feat_length){
  
  if(dim(Variant1)[1] > 0 && dim(Variant2)[1] > 0 && feat_length > 0){
    DiffPos_variant_1 <- Variant1[-which(Variant1$POS %in% Variant2$POS), ]
    DiffPos_variant_2 <- Variant2[-which(Variant2$POS %in% Variant1$POS), ]
    
    SamePos_varaint_1 <- Variant1[which(Variant1$POS %in% Variant2$POS), ]
    SamePos_varaint_2 <- Variant2[which(Variant2$POS %in% Variant1$POS), ]
    
    SamePos_change_ind <- which(unlist(sapply(1:dim(SamePos_varaint_1)[1], function(x) SamePos_varaint_1$ALT[x] != SamePos_varaint_2$ALT[SamePos_varaint_2$POS == SamePos_varaint_1$POS[x]])))
    # if(length(aa) != 0){}
    
    weight_diffPos <- (dim(DiffPos_variant_1)[1] - dim(DiffPos_variant_2)[1]) / feat_length
    weight_allPos <- (dim(DiffPos_variant_1)[1] + dim(DiffPos_variant_2)[1] + length(SamePos_change_ind)) / feat_length
    
    if(weight_diffPos > 0 ){
      log_weight_diffPos <- -log(weight_diffPos)
    }else if(weight_diffPos < 0){
      log_weight_diffPos <- log(-weight_diffPos)
    }else log_weight_diffPos <- NA
    
    log_weight_allPos <- ifelse(weight_allPos == 0, NA, log(weight_allPos))
      
    num_diffPos_change <- dim(DiffPos_variant_1)[1] + dim(DiffPos_variant_2)[1]
    num_samePos_change <- length(SamePos_change_ind)
    
  }else{
    weight_allPos <- NA
    weight_diffPos <- NA
    log_weight_allPos <- NA
    log_weight_diffPos <- NA
    num_diffPos_change <- NA
    num_samePos_change <- NA
    
  }
   
  return(list(allPos_score = weight_allPos,
              diffPos_score = weight_diffPos,
              log_allPos_score = log_weight_allPos,
              log_diffPos_score = log_weight_diffPos,
              diffPos_change = num_diffPos_change,
              samePos_change = num_samePos_change))
}

###############
# Extract SnpEff Annotation from VCF file
###############
extract_annotation <- function(annotation, Chrom_meta_info, Chrom_variants_info1){
  
  if(is.na(Chrom_variants_info1) || is.na(Chrom_meta_info)) stop("Empty info!")
  annotation_id <- which(grepl(annotation, Chrom_meta_info))
  ann_colnames <- gsub(" ", "", unlist(strsplit(unlist(strsplit(Chrom_meta_info[annotation_id],"'"))[2], "\\|")))
  ann_colnames <- c("PredictionID", ann_colnames)
  num_columns <- length(ann_colnames)
  
  num_rows <- length(Chrom_variants_info1)
  
  Chrom_ANN_info <- as.data.frame(matrix(data = "", nrow = num_rows, ncol = num_columns), stringsAsFactors = FALSE)
  colnames(Chrom_ANN_info) <- ann_colnames
  
#  for(id in 1:num_rows){
    all_info <- unlist(strsplit(Chrom_variants_info1, ";"))
    all_info <- gsub("\\(", "", all_info)
    all_info <- gsub("\\)", "", all_info)
    ann_index <- which(grepl(annotation, all_info))
    
    allele_info <- unlist(strsplit(all_info[ann_index], ","))
    
    if(length(allele_info) == 1){ # if there are only one allele info
      ann_info <- c(1, unlist(strsplit(unlist(strsplit(all_info[ann_index], split = paste0(annotation,"=")))[2], "\\|")))
      if(length(ann_info) < num_columns) {
        ann_info <- c(ann_info, rep("", num_columns - length(ann_info)))
      }
      names(ann_info) <- ann_colnames
      
      Chrom_ANN_info <- t(as.data.frame(ann_info, row.names = NULL, stringsAsFactors = FALSE))
    }else{ # if there are more than one allele info
      ann_info <- c(1, unlist(strsplit(unlist(strsplit(allele_info[1], split = paste0(annotation,"=")))[2], "\\|")))
      if(length(ann_info) < num_columns) {
        ann_info <- c(ann_info, rep("", num_columns - length(ann_info)))
      }
      ann_info <- t(as.data.frame(ann_info, row.names = NULL, stringsAsFactors = FALSE))
      for(kd in 2:length(allele_info)){
        temp_info <- t(as.data.frame(c(kd, unlist(strsplit(allele_info[kd], "\\|"))), 
                                     row.names = NULL, stringsAsFactors = FALSE))
        if(dim(temp_info)[2] < num_columns) {
          ann_info <- rbind(ann_info, cbind(temp_info, rep("", num_columns - dim(temp_info)[2])))
        } else ann_info <- rbind(ann_info, temp_info)
      }
      colnames(ann_info) <- ann_colnames
      rownames(ann_info) <- NULL
      Chrom_ANN_info <- ann_info # as.data.frame(ann_info, row.names = NULL, stringsAsFactors = FALSE)
    }
    # Chrom_ANN_info[id, ] <- ann_info
    
#  }
  return(Chrom_ANN_info)
}

##################
# Parse Info from VCF files
##################
parse_meta_info <- function(Chrom_meta_info, Chrom_variants_info){
  
  num_rows <- length(Chrom_variants_info)
  num_info_ID <- length(Chrom_meta_info)
  
  INFO_colnames <- sapply(1:(num_info_ID -3), function(x) unlist(strsplit(unlist(strsplit(Chrom_meta_info[x], ","))[1],"="))[3])
    
  original_INFO <- as.data.frame(matrix(data = NA, nrow = num_rows, ncol = num_info_ID - 3))
  colnames(original_INFO) <- INFO_colnames
  
  ANNotation_info <- as.data.frame(matrix("", nrow = num_rows*3, ncol = 18), 
                                   row.names = NULL, stringsAsFactors = FALSE)
  LOF_info <- NULL
  NMD_info <- NULL
  ann_counter <- 0
  for(id in 1:num_rows){
    all_info <- unlist(strsplit(Chrom_variants_info[id], ";"))
    ann_index <-which(grepl("ANN=", all_info))
    ann_counter <- ann_counter + 1
    
    if(ann_index > 1){ # if there is original info, in addition to ANNotation
      
      for(jd in 1:(ann_index-1)){
        temp_info <- unlist(strsplit(all_info[jd], "="))
        original_INFO[id, which(INFO_colnames %in% temp_info[1])] <- temp_info[2]
      }
      
      temp_info <- extract_annotation("ANN", Chrom_meta_info, Chrom_variants_info[id])
      #ANNotation_info <- rbind(ANNotation_info, cbind(index = rep(id, dim(temp_info)[1]), temp_info))
      repeat_num <- dim(temp_info)[1]
      if((ann_counter+repeat_num) < dim(ANNotation_info)[1]){
        ANNotation_info[ann_counter:(ann_counter+repeat_num), ] <- cbind(index = rep(id, repeat_num), temp_info)
      }else if(ann_counter < dim(ANNotation_info)[1]){
        within_ids <-  dim(ANNotation_info)[1] - ann_counter
        ANNotation_info[ann_counter:dim(ANNotation_info)[1], ] <- cbind(index = rep(id, within_ids), temp_info[1:within_ids, ])
        ANNotation_info <- rbind(ANNotation_info, cbind(index = rep(id, (repeat_num-within_ids), temp_info[(within_ids+1):repeat_num,])))
        ann_counter <- dim(ANNotation_info)[1]
      }else{
        ANNotation_info <- rbind(ANNotation_info, cbind(index = rep(id, repeat_num), temp_info))
        ann_counter <- dim(ANNotation_info)[1]
      }
      
      if(ann_index < length(all_info)){ # if there is LOF information
        temp_info <- cbind(id, extract_annotation("LOF", Chrom_meta_info, Chrom_variants_info[id]))
        LOF_info <- rbind(LOF_info, cbind(index = rep(id, dim(temp_info)[1]), temp_info))
      }
      if(length(all_info) - ann_index == 2){ # if there is NMD information
        temp_info <- cbind(id, extract_annotation("NMD", Chrom_meta_info, Chrom_variants_info[id]))
        NMD_info <- rbind(NMD_info, cbind(index = rep(id, dim(temp_info)[1]), temp_info))
      }
      
    }else if(ann_index == 1){ # if there is no original info
      
      temp_info <- extract_annotation("ANN", Chrom_meta_info, Chrom_variants_info[id])
      # ANNotation_info <- rbind(ANNotation_info, cbind(index = rep(id, dim(temp_info)[1]), temp_info))
      repeat_num <- dim(temp_info)[1]
      if((ann_counter+repeat_num) < dim(ANNotation_info)[1]){
        ANNotation_info[ann_counter:(ann_counter+repeat_num), ] <- cbind(index = rep(id, repeat_num), temp_info)
      }else if(ann_counter < dim(ANNotation_info)[1]){
        within_ids <-  dim(ANNotation_info)[1] - ann_counter
        ANNotation_info[ann_counter:dim(ANNotation_info)[1], ] <- cbind(index = rep(id, within_ids), temp_info[1:within_ids, ])
        ANNotation_info <- rbind(ANNotation_info, cbind(index = rep(id, (repeat_num-within_ids), temp_info[(within_ids+1):repeat_num,])))
        ann_counter <- dim(ANNotation_info)[1]
      }else{
        ANNotation_info <- rbind(ANNotation_info, cbind(index = rep(id, repeat_num), temp_info))
        ann_counter <- dim(ANNotation_info)[1]
      }
      
      if(ann_index < length(all_info)){ # if there is LOF information
        temp_info <- cbind(id, extract_annotation("LOF", Chrom_meta_info, Chrom_variants_info[id]))
        LOF_info <- rbind(LOF_info, cbind(index = rep(id, dim(temp_info)[1]), temp_info))
      }
      if(length(all_info) - ann_index == 2){ # if there is NMD information
        temp_info <- cbind(id, extract_annotation("NMD", Chrom_meta_info, Chrom_variants_info[id]))
        NMD_info <- rbind(NMD_info, cbind(index = rep(id, dim(temp_info)[1]), temp_info))
      }
    } #else{ # if there's no ANNotation info
      
   # }
    
  }
  return(list(original_INFO = original_INFO,
              Annotation_INFO = ANNotation_info,
              LossOfFunction_IFO = LOF_info,
              NonsenseMedDecay_INFO = NMD_info))
}

#############
# Visualization
#############

##########
# Sort genes by decreasing/increasing order of the designated scores
##########
sort_chr_genes <- function(chr_stats, feature_region, score_type, decrease = TRUE, rm.na = TRUE){
  
  if(grepl("_", feature_region)){
    feature_in_interest <- feature_region
  }else{
    feature_in_interest <- paste0(feature_region, "_", score_type)
  }
  
  feat_ind <- which(colnames(chr_stats) %in% feature_in_interest)
  
  # sort by decreasing order
  eval(parse(text = paste0("ordered_chr_stats <- chr_stats[order(chr_stats$", feature_in_interest, ", decreasing = ", decrease, ", na.last = TRUE), c(1,", feat_ind,")]")))
 
  if(rm.na){# remove genes with no changes
     eval(parse(text = paste0("ordered_chr_stats <- ordered_chr_stats[!is.na(ordered_chr_stats$", feature_in_interest, "), ]")))
     
  }
  
  # reorder GeneNames
  ordered_chr_stats$GeneName <- factor(ordered_chr_stats$GeneName, levels = factor(ordered_chr_stats$GeneName))
  
  if(rm.na){
    
  }
  
  return(ordered_chr_stats)
}

###################
# plot 
###################
plot.chrVariants <- function(diff_num, feature_region, all_diff="diff"){
  require(ggplot2)
  
  feature_in_interest <- paste0(feature_region, "_", all_diff)
  feat_ind <- which(colnames(diff_num) %in% feature_in_interest)
  
  # sort by decreasing order
  eval(parse(text = paste0("ordered_RD_pair <- diff_num[order(diff_num$", feature_in_interest, ", decreasing = TRUE, na.last = TRUE), c(1,", feat_ind,")]")))
  # remove genes with no changes
  eval(parse(text = paste0("ordered_RD_pair <- ordered_RD_pair[!is.na(ordered_RD_pair$", feature_in_interest, "), ]")))
  # reorder GeneNames
  ordered_RD_pair$GeneName <- factor(ordered_RD_pair$GeneName, levels = factor(ordered_RD_pair$GeneName))
  # ggplot parameters
  eval(parse(text = paste0("p_RD_pair_diff <- ggplot(data = ordered_RD_pair, aes(x = GeneName, y = ", feature_in_interest, "))")))# + 
  if(length(ordered_RD_pair$GeneName) > 0){
    p_RD_pair_diff <- p_RD_pair_diff + geom_point()  + 
      theme(text = element_text(size=3), axis.text.x = element_text(angle = 60, hjust = 1)) +
      # ggtitle("All genes with variants") +
      annotate("text", x = length(ordered_RD_pair$GeneName)/2, y = 11, size = 10, label = paste0("All genes with variants (", feature_region, ")")) +
      ylab("Normalized log-scale score")
  } else{
    p_RD_pair_diff <- p_RD_pair_diff + geom_point()  + 
      theme(text = element_text(size=3), axis.text.x = element_text(angle = 60, hjust = 1)) +
      # ggtitle("All genes with variants") +
      annotate("text", x = 1, y = 11, size = 10, label = paste0("All genes with variants (", feature_region, ")")) +
      annotate("text", x = 1, y = 9, size = 7, label = paste0("(No different variants detected)")) +
      ylab("Normalized log-scale score")
  }
  
  
  ## Donor dominant gene variants
  eval(parse(text = paste0("D_prime_genes <- ordered_RD_pair[ordered_RD_pair$", feature_in_interest," > 0 ,]")))
  D_prime_genes$GeneName <- factor(D_prime_genes$GeneName, levels = factor(D_prime_genes$GeneName))
  eval(parse(text = paste0("p_D_prime_diff <- ggplot(data = D_prime_genes, aes(x = GeneName, y = ", feature_in_interest, "))")))# + 
  if(length(D_prime_genes$GeneName) >0){
    p_D_prime_diff <- p_D_prime_diff + geom_point()  + 
      theme(text = element_text(size=4.5), axis.text.x = element_text(angle = 60, hjust = 1)) +
      # ggtitle("Genes with donor-dominant variants") +
      annotate("text", x = length(D_prime_genes$GeneName)/2, y = 11, size = 10, label = paste0("Genes with donor-dominant variants (", feature_region, ")")) +
      ylab("Log-scale score")
  }else{
    p_D_prime_diff <- p_D_prime_diff + geom_point()  + 
      theme(text = element_text(size=4.5), axis.text.x = element_text(angle = 60, hjust = 1)) +
      # ggtitle("Genes with donor-dominant variants") +
      annotate("text", x = 1, y = 11, size = 10, label = paste0("Genes with donor-dominant variants (", feature_region, ")")) +
      annotate("text", x = 1, y = 9, size = 7, label = paste0("(No different variants detected)")) +
      ylab("Log-scale score")
  }
  
  
  ## Recipient dominant gene variants
  eval(parse(text = paste0("R_prime_genes <- ordered_RD_pair[ordered_RD_pair$", feature_in_interest, " < 0 ,]")))
  # change the score to positive - will reflect the correct order
  eval(parse(text = paste0("R_prime_genes$total_diff <- -R_prime_genes$", feature_in_interest)))
  eval(parse(text = paste0("R_prime_genes <- R_prime_genes[order(R_prime_genes$", feature_in_interest, ", decreasing = T, na.last = T), ]")))
  R_prime_genes$GeneName <- factor(R_prime_genes$GeneName, levels = factor(R_prime_genes$GeneName))
  eval(parse(text = paste0("p_R_prime_diff <- ggplot(data = R_prime_genes, aes(x = GeneName, y = ", feature_in_interest, "))"))) # + 
  if(length(R_prime_genes$GeneName) > 0){
    p_R_prime_diff <- p_R_prime_diff + geom_point()  + 
      theme(text = element_text(size=4.5), axis.text.x = element_text(angle = 60, hjust = 1)) +
      #ggtitle("Genes with Recipient-dominant variants") +
      annotate("text", x = length(R_prime_genes$GeneName)/2, y = 11, size = 10, label = paste0("Genes with Recipient-dominant variants (", feature_region, ")")) +
      ylab("Log-scale score")
  }else{
    p_R_prime_diff <- p_R_prime_diff + geom_point()  + 
      theme(text = element_text(size=4.5), axis.text.x = element_text(angle = 60, hjust = 1)) +
      #ggtitle("Genes with Recipient-dominant variants") +
      annotate("text", x = 1, y = 11, size = 10, label = paste0("Genes with Recipient-dominant variants (", feature_region, ")")) +
      annotate("text", x = 1, y = 9, size = 7, label = paste0("(No different variants detected)")) +
      ylab("Log-scale score")
  }
  
  
#   multiplot(p_total_diff, pd_diff, pr_diff, cols = 1)
#   dev.off()
  
#   p <- ggplot(diff_num, aes(factor(cyl), mpg))
#   
#   p + geom_boxplot()
  return(list(ordered_R_D_diff_variants = ordered_RD_pair, 
              ordered_D_prime_genes = D_prime_genes,
              ordered_R_prime_genes = R_prime_genes,
              p_RD_pair = p_RD_pair_diff,
              p_D_prime = p_D_prime_diff,
              p_R_prime = p_R_prime_diff))
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}