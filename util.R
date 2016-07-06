gene_variant_stats <- function(geneInfo, geneVariants){
  
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