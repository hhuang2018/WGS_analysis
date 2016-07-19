gene_variant_stats <- function(geneInfo, geneVariants){
  #
  # returns the variant statistics by functional regions (introns, exons)
  #
  
  startIndex <- as.integer(unlist(strsplit(as.character(geneInfo$exonStarts), ",")))
  endIndex <- as.integer(unlist(strsplit(as.character(geneInfo$exonEnds), ",")))
  
  RegionIndex <- vector(mode = "numeric", length = 2*length(startIndex))
  RegionIndex[seq(from = 1, to = 2*length(startIndex), by = 2)] <- startIndex
  RegionIndex[seq(from = 2, to = 2*length(startIndex), by = 2)] <- endIndex
  
  total_variants <- geneVariants[(geneVariants$POS >= geneInfo$txStart) & (geneVariants$POS <= geneInfo$txEnd), ]
  total_variants <- total_variants[which(total_variants$FILTER == "PASS"), ]
  total_num <- dim(total_variants)[1]
  
  RegionNumber <- findInterval(total_variants$POS, RegionIndex)
  total_variants$ID <- sapply(RegionNumber, function(x) RegionNumber2RegionName(x, length(RegionIndex)))
  
  intron_variant_index <- which(RegionNumber %% 2 == 0)
  exon_variant_index <- which(RegionNumber %% 2 == 1)
  insert_variant_index <- which(nchar(total_variants$REF) < nchar(total_variants$ALT))
  deletion_variant_index <- which(nchar(total_variants$REF) > nchar(total_variants$ALT))
  
  return(list(all_variants = total_variants, 
              intron_variants = intron_variant_index,
              exon_variants = exon_variant_index,
              insert_variants = insert_variant_index,
              deletion_variants = deletion_variant_index,
              variants_stats = data.frame(total_num = total_num,
                                          intron_variants_num = length(intron_variant_index),
                                          exon_variants_num = length(exon_variant_index), 
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

###############
# Extract SnpEff Annotation from VCF file
###############
extract_annotation <- function(annotation, Chrom_meta_info, Chrom_variants_info){
  
  if(is.na(Chrom_variants_info) || is.na(Chrom_meta_info)) stop("Empty info!")
  annotation_id <- which(grepl(annotation, Chrom_meta_info))
  ann_colnames <- gsub(" ", "", unlist(strsplit(unlist(strsplit(Chrom_meta_info[annotation_id],"'"))[2], "\\|")))
  ann_colnames <- c("PredictionID", ann_colnames)
  num_columns <- length(ann_colnames)
  
  num_rows <- length(Chrom_variants_info)
  
  Chrom_ANN_info <- as.data.frame(matrix(data = "", nrow = num_rows, ncol = num_columns), stringsAsFactors = FALSE)
  colnames(Chrom_ANN_info) <- ann_colnames
  
#  for(id in 1:num_rows){
    all_info <- unlist(strsplit(Chrom_variants_info, ";"))
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
    }else{ # if there are more than one allele info
      ann_info <- c(1, unlist(strsplit(unlist(strsplit(allele_info[1], split = paste0(annotation,"=")))[2], "\\|")))
      ann_info <- rbind(ann_info, t(sapply(2:length(allele_info), function(x) c(x, unlist(strsplit(allele_info[x], "\\|"))))))
      if(dim(ann_info)[2] < num_columns) {
        ann_info <- cbind(ann_info, rep("", num_columns -dim(ann_info)[2]))
      }
      colnames(ann_info) <- ann_colnames
    }
    # Chrom_ANN_info[id, ] <- ann_info
    Chrom_ANN_info <- ann_info
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
  
  ANNotation_info <- NULL
  LOF_info <- NULL
  NMD_info <- NULL
  
  for(id in 1:num_rows){
    all_info <- unlist(strsplit(Chrom_variants_info[id], ";"))
    ann_index <-which(grepl("ANN=", all_info))
    
    if(ann_index > 1){ # if there is original info, in addition to ANNotation
      
      for(jd in 1:(ann_index-1)){
        temp_info <- unlist(strsplit(all_info[jd], "="))
        original_INFO[id, which(INFO_colnames %in% temp_info[1])] <- temp_info[2]
      }
      
      temp_info <- extract_annotation("ANN", Chrom_meta_info, Chrom_variants_info[id])
      ANNotation_info <- rbind(ANNotation_info, cbind(index = rep(id, dim(temp_info)[1]), temp_info))

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
      ANNotation_info <- rbind(ANNotation_info, cbind(index = rep(id, dim(temp_info)[1]), temp_info))
      
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
####
# plot.genome <- function(){
#   require(ggplot2)
#   
#   p <- ggplot(diff_num, aes(factor(cyl), mpg))
#   
#   p + geom_boxplot()
#   
# }