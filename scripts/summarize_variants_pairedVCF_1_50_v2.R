source('util.R', echo = FALSE)

source("http://bioconductor.org/biocLite.R") 
biocLite("VariantAnnotation") #install the package
library("VariantAnnotation") #load the package


load("../Data/ID_table.RData")

if (!require("vcfR", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("vcfR", dependencies = TRUE)
  library("vcfR", verbose=F, warn.conflicts =F)
}

vcf_fp <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/paired_vcf/"
vcf_file <- list.files(path = vcf_fp, pattern = "\\.vcf.gz$")
vcf_summary_output <- "/mnt/cloudbiodata_nfs_2/users/hhuang/paired_variants_summary/"
num_files <- length(vcf_file)

for(fid in 1:50){
  
  vcf_info <- read.vcfR(paste0(vcf_fp, vcf_file[fid]), verbose = FALSE)
  
  if(fid == 1){
    meta_data <- vcf_info@meta
    annotation_columns <- unlist(strsplit(unlist(strsplit(meta_data[which(grepl("ID=ANN", meta_data))], "\'"))[2],"\\|"))
    annotation_columns <- gsub(" ", "", annotation_columns)
    annotation_columns <- gsub("/", "__", annotation_columns)
    LOF_columns <- unlist(strsplit(unlist(strsplit(meta_data[which(grepl("ID=LOF", meta_data))], "\'"))[2],"\\|"))
    LOF_columns <- gsub(" ", "", LOF_columns)
    NMD_columns <- unlist(strsplit(unlist(strsplit(meta_data[which(grepl("ID=NMD", meta_data))], "\'"))[2],"\\|"))
    NMD_columns <- gsub(" ", "", NMD_columns)
    rm(meta_data)
  }
  
  for(chrom_id in 1:22){
    
    CHROM <- paste0("chr", chrom_id)
    
    chr_original_variants <- as.data.frame(vcf_info@fix[vcf_info@fix[, "CHROM"] == CHROM, c("CHROM", "POS", "REF", "ALT")], stringsAsFactors = F)
    chr_INFO_all <- vcf_info@fix[vcf_info@fix[, "CHROM"] == CHROM,"INFO"]
    chr_INFO_ANN <- strsplit(chr_INFO_all, ";")
    
    chr_gt <- as.data.frame(vcf_info@gt[vcf_info@fix[, "CHROM"] == CHROM, 2:3], stringsAsFactors = F)
    if(ID_table[which(ID_table$SeqID == as.numeric(colnames(chr_gt)[1])), "subjectType"] != "R"){ # first column should be Recipient
      temp_gt <- chr_gt
      chr_gt[, 1] <- temp_gt[, 2]
      chr_gt[, 2] <- temp_gt[, 1]
    }
    
    #### Annotation 
    ANN_indices <- which(grepl("ANN=",chr_INFO_ANN))
    if(length(ANN_indices) == dim(chr_original_variants)[1]){
      chr_variants <- chr_original_variants
    }else{
      chr_variants <- chr_original_variants[ANN_indices, ]
    }
    
    chr_ANN_possible <- sapply(ANN_indices, function(x) 
      strsplit(gsub("ANN=","",chr_INFO_ANN[[x]][grep("ANN=", chr_INFO_ANN[[x]])]), "\\,"))
    
    chr_ANN <- sapply(1:length(chr_ANN_possible), function(x) strsplit(chr_ANN_possible[[x]], "\\|"))
    
    #chr_all_info <- sapply(1:length(chr_ANN), function(x) if(length(chr_ANN[[x]][[1]])==16) cbind(chr_variants[x, ], chr_gt[x, ], do.call("rbind", chr_ANN[[x]])) else cbind(chr_variants[x, ], chr_gt[x, ], do.call("rbind", chr_ANN[[x]]), "16"=character(1))) 
    # chr_gt_reform <- do.call("rbind", sapply(1:dim(chr_gt)[1], function(x) do.call("cbind", strsplit(as.character(chr_gt[x,]), ":"))[1, ]))
    # chr_info_reformat <- do.call("rbind", chr_all_info)
    # colnames(chr_info_reformat[]) <- c(, annotation_columns)
    chr_info_reformat <- as.data.frame(matrix(data = NA, nrow = , ncol = 22), stringsAsFactors = F)
    colnames(chr_info_reformat) <- c("CHROM", "POS", "REF", "ALT", colnames(chr_gt), annotation_columns)
    for(ANN_id in 1:length(chr_ANN)){
      
      if(class(chr_ANN[[ANN_id]]) == "list"){
        temp_var <- cbind(chr_variants[ANN_id, ], chr_gt[ANN_id, ], do.call("rbind", chr_ANN[[ANN_id]]))
      }else temp_var <- cbind(chr_variants[ANN_id, ], chr_gt[ANN_id, ],chr_ANN[[ANN_id]])
      
      if(dim(temp_var)[2] < 22){
        temp_var <- cbind(temp_var, character(1))
      }
      colnames(temp_var) <- c("CHROM", "POS", "REF", "ALT", colnames(chr_gt), annotation_columns)
      
      chr_info_reformat <- rbind(chr_info_reformat, temp_var)
      
    }
    rm(temp_var, chr_variants, ANN_indices)
    
    ### LOF
    LOF_indices <- which(grepl("LOF=",chr_INFO_ANN))
    if(length(LOF_indices) == dim(chr_original_variants)[1]){
      chr_LOF_variants <- chr_original_variants
    }else{
      chr_LOF_variants <- chr_original_variants[LOF_indices, ]
    }
    
    chr_LOF_possible <- sapply(LOF_indices, function(x) 
      strsplit(gsub("LOF=","",chr_INFO_ANN[[x]][grep("LOF=", chr_INFO_ANN[[x]])]), "\\,"))
    
    chr_LOF <- sapply(1:length(chr_LOF_possible), function(x) strsplit(chr_LOF_possible[[x]], "\\|"))
    
    chr_LOF_reformat <- as.data.frame(matrix(data = NA, nrow =0, ncol = 10), stringsAsFactors = F)
    colnames(chr_LOF_reformat) <- c("CHROM", "POS", "REF", "ALT", colnames(chr_gt), LOF_columns)
    for(LOF_id in 1:length(chr_LOF)){
      
      if(class(chr_LOF[[LOF_id]]) == "list"){
        temp_var <- cbind(chr_LOF_variants[LOF_id, ], chr_gt[LOF_id, ], do.call("rbind", chr_LOF[[LOF_id]]))
      }else temp_var <- cbind(chr_LOF_variants[LOF_id, ], chr_gt[LOF_id, ], t(chr_LOF[[LOF_id]]))
      
      if(dim(temp_var)[2] < 10){
        temp_var <- cbind(temp_var, character(1))
      }
      colnames(temp_var) <- c("CHROM", "POS", "REF", "ALT", colnames(chr_gt), LOF_columns)
      
      chr_LOF_reformat <- rbind(chr_LOF_reformat, temp_var)
      
    }
    rm(temp_var, chr_LOF_variants, LOF_indices)
    ######## NMD
    NMD_indices <- which(grepl("NMD=",chr_INFO_ANN))
    if(length(NMD_indices) == dim(chr_original_variants)[1]){
      chr_NMD_variants <- chr_original_variants
    }else{
      chr_NMD_variants <- chr_original_variants[NMD_indices, ]
    }
    
    chr_NMD_possible <- sapply(NMD_indices, function(x) 
      strsplit(gsub("NMD=","",chr_INFO_ANN[[x]][grep("NMD=", chr_INFO_ANN[[x]])]), "\\,"))
    
    chr_NMD <- sapply(1:length(chr_NMD_possible), function(x) strsplit(chr_NMD_possible[[x]], "\\|"))
    
    chr_NMD_reformat <- as.data.frame(matrix(data = NA, nrow =0, ncol = 10), stringsAsFactors = F)
    colnames(chr_NMD_reformat) <- c("CHROM", "POS", "REF", "ALT", colnames(chr_gt), NMD_columns)
    for(NMD_id in 1:length(chr_NMD)){
      
      if(class(chr_NMD[[NMD_id]]) == "list"){
        temp_var <- cbind(chr_NMD_variants[NMD_id, ], chr_gt[NMD_id, ], do.call("rbind", chr_NMD[[NMD_id]]))
      }else temp_var <- cbind(chr_NMD_variants[NMD_id, ], chr_gt[NMD_id, ], t(chr_NMD[[NMD_id]]))
      
      if(dim(temp_var)[2] < 10){
        temp_var <- cbind(temp_var, character(1))
      }
      colnames(temp_var) <- c("CHROM", "POS", "REF", "ALT", colnames(chr_gt), NMD_columns)
      
      chr_NMD_reformat <- rbind(chr_NMD_reformat, temp_var)
      
    }
    rm(temp_var, chr_NMD_variants, NMD_indices)
    
    save(chr_info_reformat, chr_LOF_reformat, chr_NMD_reformat, file = paste0(vcf_summary_output, gsub(".vcf.gz", "",vcf_file[fid]), "_", CHROM, ".RData"))
    
    cat(paste0(gsub(".vcf.gz", "",vcf_file[fid]), " - ", CHROM, " finished!\n"))
    
  } ### autosomal 1~22
  
  for(chrom_id in c("X", "Y", "M")){
    
    CHROM <- paste0("chr", chrom_id)
    
    chr_original_variants <- as.data.frame(vcf_info@fix[vcf_info@fix[, "CHROM"] == CHROM, c("CHROM", "POS", "REF", "ALT")], stringsAsFactors = F)
    chr_INFO_all <- vcf_info@fix[vcf_info@fix[, "CHROM"] == CHROM,"INFO"]
    chr_INFO_ANN <- strsplit(chr_INFO_all, ";")
    
    chr_gt <- as.data.frame(vcf_info@gt[vcf_info@fix[, "CHROM"] == CHROM, 2:3], stringsAsFactors = F)
    if(ID_table[which(ID_table$SeqID == as.numeric(colnames(chr_gt)[1])), "subjectType"] != "R"){ # first column should be Recipient
      temp_gt <- chr_gt
      chr_gt[, 1] <- temp_gt[, 2]
      chr_gt[, 2] <- temp_gt[, 1]
    }
    
    #### Annotation 
    ANN_indices <- which(grepl("ANN=",chr_INFO_ANN))
    if(length(ANN_indices) == dim(chr_original_variants)[1]){
      chr_variants <- chr_original_variants
    }else{
      chr_variants <- chr_original_variants[ANN_indices, ]
    }
    
    chr_ANN_possible <- sapply(ANN_indices, function(x) 
      strsplit(gsub("ANN=","",chr_INFO_ANN[[x]][grep("ANN=", chr_INFO_ANN[[x]])]), "\\,"))
    
    chr_ANN <- sapply(1:length(chr_ANN_possible), function(x) strsplit(chr_ANN_possible[[x]], "\\|"))
    
    #chr_all_info <- sapply(1:length(chr_ANN), function(x) if(length(chr_ANN[[x]][[1]])==16) cbind(chr_variants[x, ], chr_gt[x, ], do.call("rbind", chr_ANN[[x]])) else cbind(chr_variants[x, ], chr_gt[x, ], do.call("rbind", chr_ANN[[x]]), "16"=character(1))) 
    # chr_gt_reform <- do.call("rbind", sapply(1:dim(chr_gt)[1], function(x) do.call("cbind", strsplit(as.character(chr_gt[x,]), ":"))[1, ]))
    # chr_info_reformat <- do.call("rbind", chr_all_info)
    # colnames(chr_info_reformat[]) <- c(, annotation_columns)
    chrXYM_info_reformat <- as.data.frame(matrix(data = NA, nrow =0, ncol = 22), stringsAsFactors = F)
    colnames(chrXYM_info_reformat) <- c("CHROM", "POS", "REF", "ALT", colnames(chr_gt), annotation_columns)
    for(ANN_id in 1:length(chr_ANN)){
      
      if(class(chr_ANN[[ANN_id]]) == "list"){
        temp_var <- cbind(chr_variants[ANN_id, ], chr_gt[ANN_id, ], do.call("rbind", chr_ANN[[ANN_id]]))
      }else temp_var <- cbind(chr_variants[ANN_id, ], chr_gt[ANN_id, ],chr_ANN[[ANN_id]])
      
      if(dim(temp_var)[2] < 22){
        temp_var <- cbind(temp_var, character(1))
      }
      colnames(temp_var) <- c("CHROM", "POS", "REF", "ALT", colnames(chr_gt), annotation_columns)
      
      chrXYM_info_reformat <- rbind(chrXYM_info_reformat, temp_var)
      
    }
    rm(temp_var, chr_variants, ANN_indices)
    
    ### LOF
    LOF_indices <- which(grepl("LOF=",chr_INFO_ANN))
    if(length(LOF_indices) == dim(chr_original_variants)[1]){
      chr_LOF_variants <- chr_original_variants
    }else{
      chr_LOF_variants <- chr_original_variants[LOF_indices, ]
    }
    
    chr_LOF_possible <- sapply(LOF_indices, function(x) 
      strsplit(gsub("LOF=","",chr_INFO_ANN[[x]][grep("LOF=", chr_INFO_ANN[[x]])]), "\\,"))
    
    chr_LOF <- sapply(1:length(chr_LOF_possible), function(x) strsplit(chr_LOF_possible[[x]], "\\|"))
    
    chrXYM_LOF_reformat <- as.data.frame(matrix(data = NA, nrow =0, ncol = 22), stringsAsFactors = F)
    colnames(chrXYM_LOF_reformat) <- c("CHROM", "POS", "REF", "ALT", colnames(chr_gt), LOF_columns)
    for(LOF_id in 1:length(chr_LOF)){
      
      if(class(chr_LOF[[LOF_id]]) == "list"){
        temp_var <- cbind(chr_LOF_variants[LOF_id, ], chr_gt[LOF_id, ], do.call("rbind", chr_LOF[[LOF_id]]))
      }else temp_var <- cbind(chr_LOF_variants[LOF_id, ], chr_gt[LOF_id, ], t(chr_LOF[[LOF_id]]))
      
      if(dim(temp_var)[2] < 10){
        temp_var <- cbind(temp_var, character(1))
      }
      colnames(temp_var) <- c("CHROM", "POS", "REF", "ALT", colnames(chr_gt), LOF_columns)
      
      chrXYM_LOF_reformat <- rbind(chrXYM_LOF_reformat, temp_var)
      
    }
    rm(temp_var, chr_LOF_variants, LOF_indices)
    
    
    ######## NMD
    NMD_indices <- which(grepl("NMD=",chr_INFO_ANN))
    if(length(NMD_indices) == dim(chr_original_variants)[1]){
      chr_NMD_variants <- chr_original_variants
    }else{
      chr_NMD_variants <- chr_original_variants[NMD_indices, ]
    }
    
    chr_NMD_possible <- sapply(NMD_indices, function(x) 
      strsplit(gsub("NMD=","",chr_INFO_ANN[[x]][grep("NMD=", chr_INFO_ANN[[x]])]), "\\,"))
    
    chr_NMD <- sapply(1:length(chr_NMD_possible), function(x) strsplit(chr_NMD_possible[[x]], "\\|"))
    
    chrXYM_NMD_reformat <- as.data.frame(matrix(data = NA, nrow =0, ncol = 22), stringsAsFactors = F)
    colnames(chrXYM_NMD_reformat) <- c("CHROM", "POS", "REF", "ALT", colnames(chr_gt), NMD_columns)
    for(NMD_id in 1:length(chr_NMD)){
      
      if(class(chr_NMD[[NMD_id]]) == "list"){
        temp_var <- cbind(chr_NMD_variants[NMD_id, ], chr_gt[NMD_id, ], do.call("rbind", chr_NMD[[NMD_id]]))
      }else temp_var <- cbind(chr_NMD_variants[NMD_id, ], chr_gt[NMD_id, ], t(chr_NMD[[NMD_id]]))
      
      if(dim(temp_var)[2] < 10){
        temp_var <- cbind(temp_var, character(1))
      }
      colnames(temp_var) <- c("CHROM", "POS", "REF", "ALT", colnames(chr_gt), NMD_columns)
      
      chrXYM_NMD_reformat <- rbind(chrXYM_NMD_reformat, temp_var)
      
    }
    rm(temp_var, chr_NMD_variants, NMD_indices)
    
    save(chrXYM_info_reformat, chrXYM_LOF_reformat, chrXYM_NMD_reformat,
         file = paste0(vcf_summary_output, gsub(".vcf.gz", "",vcf_file[fid]), "_", CHROM, ".RData"))
    
    cat(paste0(gsub(".vcf.gz", "",vcf_file[fid]), " - ", CHROM, " finished!\n"))
  } ### Chr
  
  
  cat("\n\n")
} ### file