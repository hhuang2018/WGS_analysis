HLA_typings <- read.csv(file = "../../2018WGS/R/R_Output/available_cases_HLA.csv", stringsAsFactors = F)
Outcome <- read.csv(file = "../../2018WGS/R/R_Output/available_cases.csv", stringsAsFactors = F)

num2rid <- function(ID){
  
  strID <- as.character(ID)
  strID_list <- strsplit(strID, '')[[1]]
  
  if(7-nchar(strID) < 3){
    addZeros <- paste0(rep('0', 7-nchar(strID)), collapse = '')
    
  }else if(7-nchar(strID) == 3){
    addZeros <- paste0(rep('0', 7-nchar(strID)), collapse = '')
    addZeros <- paste0(addZeros, '-', collapse = '-')
    
  }else{
    addZeros <- ''
    for(idxx2 in 1:(7-nchar(strID))){
      
      if(idxx2 == 3 | idxx2 == 6){
        addZeros <- paste0(addZeros, '0-', collapse = '')
      }else{
        addZeros <- paste0(addZeros, '0', collapse = '')
      }
      
    }
    
  }
  
  
  new_rid <- ''
  for(idxx in nchar(strID):1){
    
    if(idxx == nchar(strID) - 1 | idxx == nchar(strID) - 4){
      new_rid <- paste0(strID_list[idxx], '-', new_rid, collapse = '')
    }else{
      new_rid <- paste0(strID_list[idxx], new_rid, collapse = '')
    }
    
  }
  new_rid <- paste0(addZeros, new_rid, collapse = '')
  
  return(new_rid)
}

num2did <- function(ID){
  
  total_num_len = 9
  block_digt_len = 4
  
  strID <- as.character(ID)
  strID_list <- strsplit(strID, '')[[1]]
  
  if(total_num_len-nchar(strID) < block_digt_len){
    addZeros <- paste0(rep('0', total_num_len-nchar(strID)), collapse = '')
    
  }else if(total_num_len-nchar(strID) == block_digt_len){
    addZeros <- paste0(rep('0', total_num_len-nchar(strID)), collapse = '')
    addZeros <- paste0(addZeros, '-', collapse = '-')
    
  }else{
    addZeros <- ''
    for(idxx2 in 1:(total_num_len-nchar(strID))){
      
      if(idxx2 == block_digt_len | idxx2 == block_digt_len*2){
        addZeros <- paste0(addZeros, '0-', collapse = '')
      }else{
        addZeros <- paste0(addZeros, '0', collapse = '')
      }
      
    }
    
  }
  
  new_did <- ''
  for(idxx in nchar(strID):1){
    
    if(idxx == nchar(strID) - 1 | idxx == nchar(strID) - block_digt_len-1){
      new_did <- paste0(strID_list[idxx], '-', new_did, collapse = '')
    }else{
      new_did <- paste0(strID_list[idxx], new_did, collapse = '')
    }
    
  }
  new_did <- paste0(addZeros, new_did, collapse = '')

  return(new_did)
}

#num2rid(1237356)
#num2did(24335656)
# 10  - AML
AMLpatients <- Outcome[Outcome$disease == 10, c('rid', 'agvhi24', 'agvhi34', 'cgvhi')]

AMLpatients_IDs <- sapply(AMLpatients$rid, num2rid)
AMLpatients_txOutcome <- AMLpatients$agvhi24

write.csv(AMLpatients_IDs, file = 'AMLpatients_list.csv', row.names = F)
