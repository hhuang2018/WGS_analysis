source('util.R', echo = FALSE)

HLA_typings <- read.csv(file = "../HLI_hla_mg_v3.csv")

load("../Data/ID_table.RData")

available_IDs <-ID_table[ID_table$GroupID %in% ID_table$GroupID[duplicated(ID_table$GroupID)],]

avaialble_IDs_HLA_typing <- HLA_typings[(HLA_typings$nmdp_rid %in% available_IDs$R_D_ID),]
num_groups <- dim(avaialble_IDs_HLA_typing)[1]
group_type <- sapply(1:num_groups, function(x) available_IDs$Group[available_IDs$R_D_ID %in% avaialble_IDs_HLA_typing$nmdp_rid[x]])
group_type[group_type == "a"] <- "aGVHD"
group_type[group_type == "n"] <- "non-GVHD"

file_path <- "../Output/snpEff_gene_stats/"
output_dir <- "../Output/"
  
all_files <- list.files(file_path, pattern = "\\.csv$")

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
  file_content <- scan(paste0(file_path, all_files[id]), what = "character", sep = "\n")
  
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
  sub_types <- unlist(strsplit(all_files[id], "_"))
 
  if(sub_types[3] == "D"){ ## donor
    
    avail_did_index <- which(avaialble_IDs_HLA_typing$R_D_ID == as.numeric(sub_types[4]))
    
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
      
      ## Count_Effects_table
      temp_Count_Effects_table <- cbind(temp_preface_table, t(Count_Effects_table$Count))
      colnames(temp_Count_Effects_table)[5:(4 + length(Count_Effects_table$Type))] <- Count_Effects_table$Type
      
      all_Count_Effects_table_205_donor <- rbind(all_Count_Effects_table_205_donor, temp_Count_Effects_table)
      all_Count_Effects_table_donor <- rbind(all_Count_Effects_table_donor, temp_Count_Effects_table)
      
      ## Count_Region_table
      temp_Count_Region_table <- cbind(temp_preface_table, t(Count_Region_table$Count))
      colnames(temp_Count_Region_table)[5:(4 + length(Count_Region_table$Type))] <- Count_Region_table$Type
      
      all_Count_Region_table_205_donor <- rbind(all_Count_Region_table_205_donor, temp_Count_Region_table)
      all_Count_Region_table_donor <- rbind(all_Count_Region_table_donor, temp_Count_Region_table)

      
    }else{ # not within the 205 pairs
      
      
      temp_preface_table <- data.frame(GroupID = as.numeric(sub_types[2]),
                                       GroupType = sub_types[1],
                                       bmt_case_num = avaialble_IDs_HLA_typing$bmt_case_num[avail_did_index],
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
      temp_Count_Effects_table <- cbind(temp_preface_table, t(Count_Effects_table$Count))
      colnames(temp_Count_Effects_table)[5:(4 + length(Count_Effects_table$Type))] <- Count_Effects_table$Type
      
      all_Count_Effects_table_donor <- rbind(all_Count_Effects_table_donor, temp_Count_Effects_table)
      
      ## Count_Region_table
      temp_Count_Region_table <- cbind(temp_preface_table, t(Count_Region_table$Count))
      colnames(temp_Count_Region_table)[5:(4 + length(Count_Region_table$Type))] <- Count_Region_table$Type
      
      all_Count_Region_table_donor <- rbind(all_Count_Region_table_donor, temp_Count_Region_table)
      
    }
    
  }else{ ## recipient
    
    avail_rid_index <- which(avaialble_IDs_HLA_typing$R_D_ID == as.numeric(sub_types[4]))
    
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
      
      ## Count_Effects_table
      temp_Count_Effects_table <- cbind(temp_preface_table, t(Count_Effects_table$Count))
      colnames(temp_Count_Effects_table)[5:(4 + length(Count_Effects_table$Type))] <- Count_Effects_table$Type
      
      all_Count_Effects_table_205_recipients <- rbind(all_Count_Effects_table_205_recipients, temp_Count_Effects_table)
      all_Count_Effects_table_recipients <- rbind(all_Count_Effects_table_recipients, temp_Count_Effects_table)
      
      ## Count_Region_table
      temp_Count_Region_table <- cbind(temp_preface_table, t(Count_Region_table$Count))
      colnames(temp_Count_Region_table)[5:(4 + length(Count_Region_table$Type))] <- Count_Region_table$Type
      
      all_Count_Region_table_205_recipients <- rbind(all_Count_Region_table_205_recipients, temp_Count_Region_table)
      all_Count_Region_table_recipients <- rbind(all_Count_Region_table_recipients, temp_Count_Region_table)
      
      
    }else{ # not within the 205 pairs
      
      
      temp_preface_table <- data.frame(GroupID = as.numeric(sub_types[2]),
                                       GroupType = sub_types[1],
                                       bmt_case_num = avaialble_IDs_HLA_typing$bmt_case_num[avail_did_index],
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
      
      ## Count_Effects_table
      temp_Count_Effects_table <- cbind(temp_preface_table, t(Count_Effects_table$Count))
      colnames(temp_Count_Effects_table)[5:(4 + length(Count_Effects_table$Type))] <- Count_Effects_table$Type
      
      all_Count_Effects_table_recipients <- rbind(all_Count_Effects_table_recipients, temp_Count_Effects_table)
      
      ## Count_Region_table
      temp_Count_Region_table <- cbind(temp_preface_table, t(Count_Region_table$Count))
      colnames(temp_Count_Region_table)[5:(4 + length(Count_Region_table$Type))] <- Count_Region_table$Type
      
      all_Count_Region_table_recipients <- rbind(all_Count_Region_table_recipients, temp_Count_Region_table)
      
    }
    
  }
  
}

save(all_Effects_table_donor, all_Effects_table_recipients, all_Effects_table_205_donor, all_Effects_table_205_recipients,
     all_Effects_fun_table_donor, all_Effects_fun_table_recipients, all_Effects_fun_table_205_donor, all_Effects_fun_table_205_recipients,
     all_Count_Effects_table_donor, all_Count_Effects_table_recipients, all_Count_Effects_table_205_donor, all_Count_Effects_table_205_recipients,
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


