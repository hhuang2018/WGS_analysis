source('util.R', echo = FALSE)
library(vcfR)

# file_path <- "../Output/snpEff_gene_stats/"
# output_dir <- "../Output/"
load("../Data/hg38_known_gene_symbols_HGNC.RData")

file_path <- "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_annotated_RefSeq/"
output_dir <- file_path

all_files <- list.files(file_path, pattern = "\\.vcf.gz$")


vcf_info <- read.vcfR("../Data/a_36_D_52668563_RefSeq_annotated_canon.vcf.gz", verbose = FALSE)

vcf_file_D <- paste0(VCF_file_dir, "n_197_D_annotated.vcf.gz")
vcf_D <- read.vcfR(vcf_file_D, verbose = FALSE)

vcf_variants <- as.data.frame(vcf_info@fix, stringsAsFactors = FALSE)
vcf_variants$POS <- as.integer(vcf_variants$POS)

# vcf_gt <- as.data.frame(vcf_info@gt, stringsAsFactors = FALSE)
# vcf_gt <- vcf_gt[(vcf_variants$FILTER == "PASS"), ] ## QC

vcf_variants <- vcf_variants[(vcf_variants$FILTER == "PASS"), ] # only look at high-quality variants

chr_meta <- vcf_info@meta
info_id <- which(grepl("INFO=", chr_meta))

# Annotation
annotation_info <- parse_meta_info(chr_meta[info_id], vcf_variants$INFO)


genetype_info <- vcf_info@fix

## Effects by impact
all_Effects_table_205_recipients <- data.frame(GroupID = numeric(0), 
                                               GroupType = character(0),
                                               bmt_case_num = numeric(0),
                                               nmdp_rid = numeric(0),
                                               stringsAsFactors = F) 
all_Effects_table_205_donor <- data.frame(GroupID = numeric(0), 
                                          GroupType = character(0),
                                          bmt_case_num = numeric(0),
                                          nmdp_did = numeric(0),
                                          stringsAsFactors = F) 

all_Effects_table_recipients <- data.frame(GroupID = numeric(0), 
                                           GroupType = character(0),
                                           bmt_case_num = numeric(0),
                                           nmdp_rid = numeric(0),
                                           stringsAsFactors = F) 
all_Effects_table_donor <- data.frame(GroupID = numeric(0), 
                                      GroupType = character(0),
                                      bmt_case_num = numeric(0),
                                      nmdp_did = numeric(0),
                                      stringsAsFactors = F) 

## Effects by functional class
all_Effects_fun_table_205_recipients <- data.frame(GroupID = numeric(0), 
                                                   GroupType = character(0),
                                                   bmt_case_num = numeric(0),
                                                   nmdp_rid = numeric(0),
                                                   stringsAsFactors = F) 
all_Effects_fun_table_205_donor <- data.frame(GroupID = numeric(0), 
                                              GroupType = character(0),
                                              bmt_case_num = numeric(0),
                                              nmdp_did = numeric(0),
                                              stringsAsFactors = F) 

all_Effects_fun_table_recipients <- data.frame(GroupID = numeric(0), 
                                               GroupType = character(0),
                                               bmt_case_num = numeric(0),
                                               nmdp_rid = numeric(0),
                                               stringsAsFactors = F) 
all_Effects_fun_table_donor <- data.frame(GroupID = numeric(0), 
                                          GroupType = character(0),
                                          bmt_case_num = numeric(0),
                                          nmdp_did = numeric(0),
                                          stringsAsFactors = F) 

## Count by effects
all_Count_Effects_table_205_recipients <- data.frame(GroupID = numeric(0), 
                                                     GroupType = character(0),
                                                     bmt_case_num = numeric(0),
                                                     nmdp_rid = numeric(0),
                                                     stringsAsFactors = F) 
all_Count_Effects_table_205_donor <- data.frame(GroupID = numeric(0), 
                                                GroupType = character(0),
                                                bmt_case_num = numeric(0),
                                                nmdp_did = numeric(0),
                                                stringsAsFactors = F) 

all_Count_Effects_table_recipients <- data.frame(GroupID = numeric(0), 
                                                 GroupType = character(0),
                                                 bmt_case_num = numeric(0),
                                                 nmdp_rid = numeric(0),
                                                 stringsAsFactors = F) 
all_Count_Effects_table_donor <- data.frame(GroupID = numeric(0), 
                                            GroupType = character(0),
                                            bmt_case_num = numeric(0),
                                            nmdp_did = numeric(0),
                                            stringsAsFactors = F) 

## Count by genomic region
all_Count_Region_table_205_recipients <- data.frame(GroupID = numeric(0), 
                                                    GroupType = character(0),
                                                    bmt_case_num = numeric(0),
                                                    nmdp_rid = numeric(0),
                                                    stringsAsFactors = F) 
all_Count_Region_table_205_donor <- data.frame(GroupID = numeric(0), 
                                               GroupType = character(0),
                                               bmt_case_num = numeric(0),
                                               nmdp_did = numeric(0),
                                               stringsAsFactors = F) 

all_Count_Region_table_recipients <- data.frame(GroupID = numeric(0), 
                                                GroupType = character(0),
                                                bmt_case_num = numeric(0),
                                                nmdp_rid = numeric(0),
                                                stringsAsFactors = F) 
all_Count_Region_table_donor <- data.frame(GroupID = numeric(0), 
                                           GroupType = character(0),
                                           bmt_case_num = numeric(0),
                                           nmdp_did = numeric(0),
                                           stringsAsFactors = F) 

for(fid in 1:length(all_files)){
  
  # id <- 1
  file_content <- scan(paste0(file_path, all_files[fid]), what = "character", sep = "\n")
  
  comment_lines <- which(grepl("#", file_content))
  
  ### Effects by impact
  index1 <- which(grepl("# Effects by impact", file_content[comment_lines]))
  
  if(length(index1) >0){
    
    Effects_colnames <-  gsub(" ", "", unlist(strsplit(file_content[comment_lines[index1] + 1], ",")))
    
    Effects_table_origin <- file_content[(comment_lines[index1] + 2):(comment_lines[index1+1]-1)]
    
    Effects_table_origin_transform <- sapply(1:length(Effects_table_origin), function(x) gsub(" ", "", unlist(strsplit(Effects_table_origin[x], ","))))
    Effects_table_origin_transform <- t(Effects_table_origin_transform)
    
    Effects_table <- as.data.frame(Effects_table_origin_transform, stringsAsFactors = F)
    colnames(Effects_table) <- Effects_colnames
    Effects_table$Count <- as.numeric(Effects_table$Count)
    Effects_table$Percent <- as.numeric(gsub("%", "", Effects_table$Percent))
    
  }
  
  ### Effects by functional class
  index2 <- which(grepl("# Effects by functional class", file_content[comment_lines]))
  
  if(length(index2) >0){
    
    Effects_fun_colnames <-  gsub(" ", "", unlist(strsplit(file_content[comment_lines[index2] + 1], ",")))
    
    Effects_fun_table_origin <- file_content[(comment_lines[index2] + 2):(comment_lines[index2+1]-1)]
    
    Effects_fun_table_origin_transform <- sapply(1:length(Effects_fun_table_origin), function(x) gsub(" ", "", unlist(strsplit(Effects_fun_table_origin[x], ","))))
    if(class(Effects_fun_table_origin_transform) == "list")   Effects_fun_table_origin_transform <- sapply(1:(length(Effects_fun_table_origin)-1), function(x) gsub(" ", "", unlist(strsplit(Effects_fun_table_origin[x], ","))))
    
    Effects_fun_table_origin_transform <- t(Effects_fun_table_origin_transform)
    
    Effects_fun_table <- as.data.frame(Effects_fun_table_origin_transform, stringsAsFactors = F)
    colnames(Effects_fun_table) <- Effects_fun_colnames
    Effects_fun_table$Count <- as.numeric(Effects_fun_table$Count)
    Effects_fun_table$Percent <- as.numeric(gsub("%", "", Effects_fun_table$Percent))
    
  }
  
  ### Count by effects
  index3 <- which(grepl("# Count by effects", file_content[comment_lines]))
  
  if(length(index3) >0){
    
    Count_Effects_colnames <-  gsub(" ", "", unlist(strsplit(file_content[comment_lines[index3] + 1], ",")))
    
    Count_Effects_table_origin <- file_content[(comment_lines[index3] + 2):(comment_lines[index3+1]-1)]
    
    Count_Effects_table_origin_transform <- sapply(1:length(Count_Effects_table_origin), function(x) gsub(" ", "", unlist(strsplit(Count_Effects_table_origin[x], ","))))
    # if(class(Count_Effects_table_origin_transform) == "list")   Count_Effects_table_origin_transform <- sapply(1:(length(Count_Effects_table_origin)-1), function(x) gsub(" ", "", unlist(strsplit(Count_Effects_table_origin[x], ","))))
    
    Count_Effects_table_origin_transform <- t(Count_Effects_table_origin_transform)
    
    Count_Effects_table <- as.data.frame(Count_Effects_table_origin_transform, stringsAsFactors = F)
    colnames(Count_Effects_table) <- Count_Effects_colnames
    Count_Effects_table$Count <- as.numeric(Count_Effects_table$Count)
    Count_Effects_table$Percent <- as.numeric(gsub("%", "", Count_Effects_table$Percent))
    
  }
  
  ### Count by genomic region
  index4 <- which(grepl("# Count by genomic region", file_content[comment_lines]))
  
  if(length(index4) >0){
    
    Count_Region_colnames <-  gsub(" ", "", unlist(strsplit(file_content[comment_lines[index4] + 1], ",")))
    
    Count_Region_table_origin <- file_content[(comment_lines[index4] + 2):(comment_lines[index4+1]-1)]
    
    Count_Region_table_origin_transform <- sapply(1:length(Count_Region_table_origin), function(x) gsub(" ", "", unlist(strsplit(Count_Region_table_origin[x], ","))))
    # if(class(Count_Region_table_origin_transform) == "list")   Count_Region_table_origin_transform <- sapply(1:(length(Count_Region_table_origin)-1), function(x) gsub(" ", "", unlist(strsplit(Count_Region_table_origin[x], ","))))
    
    Count_Region_table_origin_transform <- t(Count_Region_table_origin_transform)
    
    Count_Region_table <- as.data.frame(Count_Region_table_origin_transform, stringsAsFactors = F)
    colnames(Count_Region_table) <- Count_Region_colnames
    Count_Region_table$Count <- as.numeric(Count_Region_table$Count)
    Count_Region_table$Percent <- as.numeric(gsub("%", "", Count_Region_table$Percent))
    
  }
  
  ############# Add to the total table
  sub_types <- unlist(strsplit(all_files[fid], "_"))
  
  if(sub_types[3] == "D"){ ## donor
    
    avail_did_index <- which(avaialble_IDs_HLA_typing$nmdp_id == as.numeric(sub_types[4]))
    
    if(length(avail_did_index) > 0){ # within 205 pairs
      
      temp_preface_table <- data.frame(GroupID = as.numeric(sub_types[2]),
                                       GroupType = sub_types[1],
                                       bmt_case_num = avaialble_IDs_HLA_typing$bmt_case_num[avail_did_index],
                                       nmdp_did = as.numeric(sub_types[4]),
                                       stringsAsFactors = F)
      
      ## Effects_table
      temp_Effects_table <- cbind(temp_preface_table, t(Effects_table$Count))
      colnames(temp_Effects_table)[5:(4 +  length(Effects_table$Type))] <- Effects_table$Type
      
      all_Effects_table_205_donor <- rbind(all_Effects_table_205_donor, temp_Effects_table)
      all_Effects_table_donor <- rbind(all_Effects_table_donor, temp_Effects_table)
      
      ## Effects_fun_table
      temp_Effects_fun_table <- cbind(temp_preface_table, t(Effects_fun_table$Count))
      colnames(temp_Effects_fun_table)[5:(4 +  length(Effects_fun_table$Type))] <- Effects_fun_table$Type
      
      all_Effects_fun_table_205_donor <- rbind(all_Effects_fun_table_205_donor, temp_Effects_fun_table)
      all_Effects_fun_table_donor <- rbind(all_Effects_fun_table_donor, temp_Effects_fun_table)
      
      # ## Count_Effects_table
      # temp_Count_Effects_table <- cbind(temp_preface_table, t(Count_Effects_table$Count))
      # colnames(temp_Count_Effects_table)[5:(4 + length(Count_Effects_table$Type))] <- Count_Effects_table$Type
      # 
      # all_Count_Effects_table_205_donor <- rbind(all_Count_Effects_table_205_donor, temp_Count_Effects_table)
      # all_Count_Effects_table_donor <- rbind(all_Count_Effects_table_donor, temp_Count_Effects_table)
      
      ## Count_Region_table
      temp_Count_Region_table <- cbind(temp_preface_table, t(Count_Region_table$Count))
      colnames(temp_Count_Region_table)[5:(4 + length(Count_Region_table$Type))] <- Count_Region_table$Type
      
      all_Count_Region_table_205_donor <- rbind(all_Count_Region_table_205_donor, temp_Count_Region_table)
      all_Count_Region_table_donor <- rbind(all_Count_Region_table_donor, temp_Count_Region_table)
      
      
    }else{ # not within the 205 pairs
      
      # avail_did_index <- which(ID_table$R_D_ID == as.numeric(sub_types[4]))
      
      temp_preface_table <- data.frame(GroupID = as.numeric(sub_types[2]),
                                       GroupType = sub_types[1],
                                       bmt_case_num = 0,
                                       nmdp_did = as.numeric(sub_types[4]),
                                       stringsAsFactors = F)
      
      ## Effects_table
      temp_Effects_table <- cbind(temp_preface_table, t(Effects_table$Count))
      colnames(temp_Effects_table)[5:(4 +  length(Effects_table$Type))] <- Effects_table$Type
      
      all_Effects_table_donor <- rbind(all_Effects_table_donor, temp_Effects_table)
      
      ## Effects_fun_table
      temp_Effects_fun_table <- cbind(temp_preface_table, t(Effects_fun_table$Count))
      colnames(temp_Effects_fun_table)[5:(4 +  length(Effects_fun_table$Type))] <- Effects_fun_table$Type
      
      all_Effects_fun_table_donor <- rbind(all_Effects_fun_table_donor, temp_Effects_fun_table)
      
      ## Count_Effects_table
      # temp_Count_Effects_table <- cbind(temp_preface_table, t(Count_Effects_table$Count))
      # colnames(temp_Count_Effects_table)[5:(4 + length(Count_Effects_table$Type))] <- Count_Effects_table$Type
      # 
      # if(identical(colnames(all_Count_Effects_table_donor)[-(1:4)], Count_Effects_table$Type)){
      #   
      #   all_Count_Effects_table_donor <- rbind(all_Count_Effects_table_donor, temp_Count_Effects_table)
      #   
      # }else{
      #   
      #   num_types <- length(Count_Effects_table$Type)
      #   for(tid in 1:num_types){
      #     
      #     eval(parse(text = paste0("all_Count_Effects_table_donor$")))
      #     
      #   }
      #   
      #   temp_colname_ID <- which(Count_Effects_table$Type %in% colnames(all_Count_Effects_table_donor)[-(1:4)])
      #   all_table_colname_ID <- which(colnames(all_Count_Effects_table_donor)[-(1:4)] %in% Count_Effects_table$Type)
      #   
      #   temp_Count_Effects_table_ordered <- temp_Count_Effects_table[c(1:4, temp_colname_ID+4)]
      #   
      #   if(length(temp_colname_ID) < length(Count_Effects_table$Type)){ # if the temp table have new columns, then add them to the end
      #     
      #     temp_new_col_table <- cbind(temp_preface_table, t(Count_Effects_table$Count))
      #     colnames(temp_Count_Effects_table)[5:(4 + length(Count_Effects_table$Type))] <- Count_Effects_table$Type
      #     
      #     temp_Count_Effects_table_ordered <- cbind(temp_Count_Effects_table_ordered, 
      #                                               t(data.frame(Count_Effects_table$Count[-temp_colname_ID], row.names = Count_Effects_table$Type[-temp_colname_ID], stringsAsFactors = F)))
      #     
      #     
      #   }
      #   
      #   all_Count_Effects_table_donor <- rbind(all_Count_Effects_table_donor, temp_Count_Effects_table)
      #   
      # }
      # 
      
      
      
      
      ## Count_Region_table
      temp_Count_Region_table <- cbind(temp_preface_table, t(Count_Region_table$Count))
      colnames(temp_Count_Region_table)[5:(4 + length(Count_Region_table$Type))] <- Count_Region_table$Type
      
      all_Count_Region_table_donor <- rbind(all_Count_Region_table_donor, temp_Count_Region_table)
      
    }
    
  }else{ ## recipient
    
    avail_rid_index <- which(avaialble_IDs_HLA_typing$nmdp_rid == as.numeric(sub_types[4]))
    
    if(length(avail_rid_index) > 0){ # within 205 pairs
      
      temp_preface_table <- data.frame(GroupID = as.numeric(sub_types[2]),
                                       GroupType = sub_types[1],
                                       bmt_case_num = avaialble_IDs_HLA_typing$bmt_case_num[avail_rid_index],
                                       nmdp_did = as.numeric(sub_types[4]),
                                       stringsAsFactors = F)
      
      ## Effects_table
      temp_Effects_table <- cbind(temp_preface_table, t(Effects_table$Count))
      colnames(temp_Effects_table)[5:(4 +  length(Effects_table$Type))] <- Effects_table$Type
      
      all_Effects_table_205_recipients <- rbind(all_Effects_table_205_recipients, temp_Effects_table)
      all_Effects_table_recipients <- rbind(all_Effects_table_recipients, temp_Effects_table)
      
      ## Effects_fun_table
      temp_Effects_fun_table <- cbind(temp_preface_table, t(Effects_fun_table$Count))
      colnames(temp_Effects_fun_table)[5:(4 +  length(Effects_fun_table$Type))] <- Effects_fun_table$Type
      
      all_Effects_fun_table_205_recipients <- rbind(all_Effects_fun_table_205_recipients, temp_Effects_fun_table)
      all_Effects_fun_table_recipients <- rbind(all_Effects_fun_table_recipients, temp_Effects_fun_table)
      
      # ## Count_Effects_table
      # temp_Count_Effects_table <- cbind(temp_preface_table, t(Count_Effects_table$Count))
      # colnames(temp_Count_Effects_table)[5:(4 + length(Count_Effects_table$Type))] <- Count_Effects_table$Type
      # 
      # all_Count_Effects_table_205_recipients <- rbind(all_Count_Effects_table_205_recipients, temp_Count_Effects_table)
      # all_Count_Effects_table_recipients <- rbind(all_Count_Effects_table_recipients, temp_Count_Effects_table)
      
      ## Count_Region_table
      temp_Count_Region_table <- cbind(temp_preface_table, t(Count_Region_table$Count))
      colnames(temp_Count_Region_table)[5:(4 + length(Count_Region_table$Type))] <- Count_Region_table$Type
      
      all_Count_Region_table_205_recipients <- rbind(all_Count_Region_table_205_recipients, temp_Count_Region_table)
      all_Count_Region_table_recipients <- rbind(all_Count_Region_table_recipients, temp_Count_Region_table)
      
      
    }else{ # not within the 205 pairs
      
      # avail_rid_index <- which(ID_table$R_D_ID == as.numeric(sub_types[4]))
      
      temp_preface_table <- data.frame(GroupID = as.numeric(sub_types[2]),
                                       GroupType = sub_types[1],
                                       bmt_case_num = 0,
                                       nmdp_did = as.numeric(sub_types[4]),
                                       stringsAsFactors = F)
      
      ## Effects_table
      temp_Effects_table <- cbind(temp_preface_table, t(Effects_table$Count))
      colnames(temp_Effects_table)[5:(4 +  length(Effects_table$Type))] <- Effects_table$Type
      
      all_Effects_table_recipients <- rbind(all_Effects_table_recipients, temp_Effects_table)
      
      ## Effects_fun_table
      temp_Effects_fun_table <- cbind(temp_preface_table, t(Effects_fun_table$Count))
      colnames(temp_Effects_fun_table)[5:(4 +  length(Effects_fun_table$Type))] <- Effects_fun_table$Type
      
      all_Effects_fun_table_recipients <- rbind(all_Effects_fun_table_recipients, temp_Effects_fun_table)
      
      # ## Count_Effects_table
      # temp_Count_Effects_table <- cbind(temp_preface_table, t(Count_Effects_table$Count))
      # colnames(temp_Count_Effects_table)[5:(4 + length(Count_Effects_table$Type))] <- Count_Effects_table$Type
      # 
      # all_Count_Effects_table_recipients <- rbind(all_Count_Effects_table_recipients, temp_Count_Effects_table)
      
      ## Count_Region_table
      temp_Count_Region_table <- cbind(temp_preface_table, t(Count_Region_table$Count))
      colnames(temp_Count_Region_table)[5:(4 + length(Count_Region_table$Type))] <- Count_Region_table$Type
      
      all_Count_Region_table_recipients <- rbind(all_Count_Region_table_recipients, temp_Count_Region_table)
      
    }
    
  }
  
}

save(all_Effects_table_donor, all_Effects_table_recipients, all_Effects_table_205_donor, all_Effects_table_205_recipients,
     all_Effects_fun_table_donor, all_Effects_fun_table_recipients, all_Effects_fun_table_205_donor, all_Effects_fun_table_205_recipients,
     # all_Count_Effects_table_donor, all_Count_Effects_table_recipients, all_Count_Effects_table_205_donor, all_Count_Effects_table_205_recipients,
     all_Count_Region_table_donor, all_Count_Region_table_recipients, all_Count_Region_table_205_donor, all_Count_Region_table_205_recipients,
     file = paste0(output_dir, "all_snpEff_stats.RData"))



######## The following is discarded - not an optimum method to read lines from a file
###### Use connection of files 
# con <- file(paste0(file_path, all_files[1]), open = "r") # open a conncetion
# 
# counter <- 0
# 
# while (length(oneLine <- readLines(con, n = 1)) > 0 & counter < 201) {
#   counter <- counter + 1
#   myLine <- unlist((strsplit(oneLine, ",")))
#   print(myLine)
#   
#   # if(grepl("# Effects by impact", myLine)){  ## Count of effects by impact
#   # 
#   #   counter <- counter + 1
#   #   oneLine <- readLines(con, n = 1)
#   #   myLine <- unlist((strsplit(oneLine, " , ")))
#   #   Effects_by_impact <- data.frame()# matrix(nrow = 5, ncol = 4)
#   #   # Effects_by_impact <- as.data.frame(Effects_by_impact)
#   #   colnames(Effects_by_impact) <- myLine
#   #   counter1 <- 0
#   #   while(length(oneLine <- readLines(con, n = 1)) > 0 & !grepl("#", oneLine)){
#   #     
#   #     counter1 <- counter1 + 1
#   #     Effects_by_impact[counter1] <- unlist((strsplit(oneLine, " , ")))
#   #     
#   #   }
#   # }
#   # 
#   # grepl("# Count by effects", myLine)   ## Count of effects by type
#   # 
#   # grepl("# Count by genomic region", myLine) ## Count of effects by region
#   # 
#   
# } 
# close(con)

###### 
library(ggplot2)
load("../Output/all_snpEff_stats.RData")

## effect by impact
all_Effects_table_205_donor$SubType <- "Donor"
all_Effects_table_205_recipients$SubType <- "Recipient"

df4 <- data.frame(Type = c(rep(all_Effects_table_205_donor$SubType, dim(all_Effects_table_205_donor)[2] -5), 
                           rep(all_Effects_table_205_recipients$SubType, dim(all_Effects_table_205_recipients)[2] -5)),
                  GenomicRegion = c(rep(colnames(all_Effects_table_205_donor)[-c(1:4, 9)], each = dim(all_Effects_table_205_donor)[1]), 
                                    rep(colnames(all_Effects_table_205_recipients)[-c(1:4, 9)], each = dim(all_Effects_table_205_recipients)[1])),
                  Count = c(as.vector(as.matrix(all_Effects_table_205_donor[, 5:8])), as.vector(as.matrix(all_Effects_table_205_recipients[, 5:8]))))

# df3$GenomicRegion <- factor(df3$GenomicRegion, levels(df3$GenomicRegion)[c(3,1,2)])

ggplot(data = df4, aes(x = GenomicRegion, y = Count, fill = Type)) + 
  geom_boxplot() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Effects by impact")

sapply(5:8, function(x) t.test(all_Effects_table_205_donor[, x], all_Effects_table_205_recipients[, x]))

###
## ## Effects by functional class
all_Effects_fun_table_205_donor$SubType <- "Donor"
all_Effects_fun_table_205_recipients$SubType <- "Recipient"

df3 <- data.frame(Type = c(rep(all_Effects_fun_table_205_donor$SubType, dim(all_Effects_fun_table_205_donor)[2] -5), 
                           rep(all_Effects_fun_table_205_recipients$SubType, dim(all_Effects_fun_table_205_recipients)[2] -5)),
                  GenomicRegion = c(rep(colnames(all_Effects_fun_table_205_donor)[-c(1:4, 8)], each = dim(all_Effects_fun_table_205_donor)[1]), 
                                    rep(colnames(all_Effects_fun_table_205_recipients)[-c(1:4, 8)], each = dim(all_Effects_fun_table_205_recipients)[1])),
                  Count = c(as.vector(as.matrix(all_Effects_fun_table_205_donor[, 5:7])), as.vector(as.matrix(all_Effects_fun_table_205_recipients[, 5:7]))))

df3$GenomicRegion <- factor(df3$GenomicRegion, levels(df3$GenomicRegion)[c(3,1,2)])

ggplot(data = df3, aes(x = GenomicRegion, y = Count, fill = Type)) + 
  geom_boxplot() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Effects by functional class")

sapply(5:7, function(x) t.test(all_Effects_fun_table_205_donor[, x], all_Effects_fun_table_205_recipients[, x]))

## plot count of variants by genomic region
all_Count_Region_table_205_donor$SubType <- "Donor"
all_Count_Region_table_205_recipients$SubType <- "Recipient"

df <- data.frame(Type = c(rep(all_Count_Region_table_205_donor$SubType, dim(all_Count_Region_table_205_donor)[2] -5), 
                          rep(all_Count_Region_table_205_recipients$SubType, dim(all_Count_Region_table_205_recipients)[2] -5)),
                 GenomicRegion = c(rep(colnames(all_Count_Region_table_205_donor)[-c(1:4,17)], each = dim(all_Count_Region_table_205_donor)[1]), 
                                   rep(colnames(all_Count_Region_table_205_recipients)[-c(1:4,17)], each = dim(all_Count_Region_table_205_recipients)[1])),
                 Count = c(as.vector(as.matrix(all_Count_Region_table_205_donor[, 5:16])), as.vector(as.matrix(all_Count_Region_table_205_recipients[, 5:16]))))

df$GenomicRegion <- factor(df$GenomicRegion, levels(df$GenomicRegion)[c(5, 4, 10, 1, 2,3,6:9, 11,12)])

ggplot(data = df, aes(x = GenomicRegion, y = Count, fill = Type)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Count of variants by genomic region")


df2 <- data.frame(Type = c(rep(all_Count_Region_table_205_donor$SubType, dim(all_Count_Region_table_205_donor)[2] -9), 
                           rep(all_Count_Region_table_205_recipients$SubType, dim(all_Count_Region_table_205_recipients)[2] -9)),
                  GenomicRegion = c(rep(colnames(all_Count_Region_table_205_donor)[-c(1:4,17, 8,9,5,14)], each = dim(all_Count_Region_table_205_donor)[1]), 
                                    rep(colnames(all_Count_Region_table_205_recipients)[-c(1:4,17, 8,9,5,14)], each = dim(all_Count_Region_table_205_recipients)[1])),
                  Count = c(as.vector(as.matrix(all_Count_Region_table_205_donor[, c(6,7,10:13, 15, 16)])), as.vector(as.matrix(all_Count_Region_table_205_recipients[, c(6,7,10:13, 15, 16)]))))

ggplot(data = df2, aes(x = GenomicRegion, y = Count, fill = Type)) + geom_boxplot()

sapply(5:16, function(x) t.test(all_Count_Region_table_205_donor[, x], all_Count_Region_table_205_recipients[, x]))
