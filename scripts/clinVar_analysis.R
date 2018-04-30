## missense mutations
load("../Output/Missense_variant_stats/donor_missesense_stats_updated.RData")      #
load("../Output/Missense_variant_stats/recipient_missesense_stats_updated.RData")

## ClinVar AML genes
clinVar.AML <- read.delim("../ClinVar/clinvar_AML_result.txt") # 53 total related genes

## Extract Pathogenic/Likely Pathogenic mutations
AML_condition <- "Pathogenic"
selected_index <- which(grepl(AML_condition, clinVar.AML$Clinical.significance..Last.reviewed., 
                              ignore.case = T))

selected_index <- selected_index[-25]
# pathogenic - 34
# Likely pathogenic - 3
# Conflicting interpretations of pathogenicity, Affects - 1

Location <- data.frame(CHROM = as.character(clinVar.AML$Chromosome[selected_index]),
                       POS = as.character(clinVar.AML$Location[selected_index]),
                       REF = character(length(selected_index)),
                       ALT = character(length(selected_index)),
                       Genes = character(length(selected_index)),
                       Names = character(length(selected_index)),
                       stringsAsFactors = F)

for(id in 1:length(selected_index)){
  
  Location$Genes[id] <- as.character(clinVar.AML$Gene.s.[selected_index[id]])
  Location$Names[id] <- as.character(clinVar.AML$Name[selected_index[id]])
  
  clinName <- unlist(strsplit(as.character(clinVar.AML$Name[selected_index[id]]), ""))
  missense_flag <- which(clinName == ">")
  if(length(missense_flag) > 0){
    
    Location$REF[id] <- clinName[missense_flag-1]
    Location$ALT[id] <- clinName[missense_flag+1]
  }
  
  clinName <- unlist(strsplit(as.character(clinVar.AML$Name[selected_index[id]]), ":"))
  clinName <- unlist(strsplit(clinName[2], " "))
  duplication_flag <- grepl("dup", clinName[1], fixed = TRUE)
  insertion_flag <- grepl("ins", clinName[1], fixed = TRUE)
  deletion_flag <- grepl("del", clinName[1], fixed = TRUE)
  
  if(duplication_flag){  ## duplicates/repeats
    start_ind_temp <- gregexpr("dup", clinName[1])
    start_ind <- start_ind_temp[[1]][1] + attr(start_ind_temp[[1]], "match.length")
    # end_ind_temp <- gregexpr(" \\(", clinName[1])
    # if(end_ind_temp[[1]] > 0){
    end_ind <- nchar(clinName[1])  # end_ind_temp[[1]] - 1
    
    dup_nuc <- unlist(strsplit(clinName[1], ""))[start_ind:end_ind]
    if(!grepl("[0-9]", dup_nuc) && !grepl("[a-z]", dup_nuc)){
      Location$REF[id] <- paste0(dup_nuc, collapse = "")
      Location$ALT[id] <- paste0(dup_nuc, collapse = "")
    }
    
    # }
    
  }
  
  if(insertion_flag){ # insertion
    
    start_ind_temp <- gregexpr("ins", clinName[1])
    start_ind <- start_ind_temp[[1]][1] + attr(start_ind_temp[[1]], "match.length")
    # end_ind_temp <- gregexpr(" \\(", clinName[2])
    # end_ind <- end_ind_temp[[1]] - 1
    end_ind <- nchar(clinName[1]) 
    ins_nuc <- unlist(strsplit(clinName[1], ""))[start_ind:end_ind]
    if(!grepl("[0-9]", ins_nuc) && !grepl("[a-z]", ins_nuc)){
      Location$REF[id] <- ins_nuc[1]
      Location$ALT[id] <- paste0(ins_nuc, collapse = "")
    }
  }
  
  if(deletion_flag){ # deletion
    
    start_ind_temp <- gregexpr("del", clinName[1])
    start_ind <- start_ind_temp[[1]][1] + attr(start_ind_temp[[1]], "match.length")[1]
    # end_ind_temp <- gregexpr(" \\(", clinName[2])
    # end_ind <- end_ind_temp[[1]] - 1
    end_ind <- nchar(clinName[1]) 
    
    del_nuc <- unlist(strsplit(clinName[1], ""))[start_ind:end_ind]
    
    if(!grepl("[0-9]", del_nuc) && !grepl("[a-z]", del_nuc)){
      Location$REF[id] <- paste0(del_nuc, collapse = "")
      Location$ALT[id] <- del_nuc[1]
    }
  }
  
}

aml_locations <- Location[which(Location$REF!=""), ] # 37 to 32 with REF and ALT
num_loacations <- dim(aml_locations)[1] 
donor_aml_mutations <- data.frame(CHROM = character(num_loacations),
                                  POS = numeric(num_loacations),
                                  MutNum = numeric(num_loacations),
                                  TotalNum = numeric(num_loacations),
                                  stringsAsFactors = F)
recipient_aml_mutations <-  data.frame(CHROM = character(num_loacations),
                                       POS = numeric(num_loacations),
                                       MutNum = numeric(num_loacations),
                                       TotalNum = numeric(num_loacations),
                                       stringsAsFactors = F)
for(id in 1:num_loacations){
  
  if(grepl("-", aml_locations$POS[id])) {
    
    POS_temp <- unlist(strsplit(aml_locations$POS[id], "-"))[1]
    POS <- as.numeric(gsub(" ", "", POS_temp)) # remove white space
    
  }else POS <- as.numeric(aml_locations$POS[id])
  
  ####  donor
  chr_index <- which(donor_missense_stats$CHROM == paste0("chr", aml_locations$CHROM[id]))
  POS_index <- which(donor_missense_stats$POS[chr_index] == POS)
  if(length(POS_index) > 0){
    
    donor_aml_mutations$CHROM[id] <- paste0("chr", aml_locations$CHROM[id])
    donor_aml_mutations$POS[id] <- POS
    if(donor_missense_stats$REF[chr_index[POS_index]] == aml_locations$REF[id]){
      
      switch (aml_locations$ALT,
      "A" = {donor_aml_mutations$MutNum[id] <- donor_missense_stats$ALT.A[chr_index[POS_index]] },
      "T" = {donor_aml_mutations$MutNum[id] <- donor_missense_stats$ALT.T[chr_index[POS_index]] },
      "G" = {donor_aml_mutations$MutNum[id] <- donor_missense_stats$ALT.G[chr_index[POS_index]] },
      "C" = {donor_aml_mutations$MutNum[id] <- donor_missense_stats$ALT.C[chr_index[POS_index]] }
      )
      donor_aml_mutations$TotalNum[id] <- donor_missense_stats$NumDiff[chr_index[POS_index]]
      
    }else print(paste0("donor - ", id))
    
  }else{
    donor_aml_mutations$CHROM[id] <- paste0("chr", aml_locations$CHROM[id])
    donor_aml_mutations$POS[id] <- POS
    donor_aml_mutations$Num[id] <- 0
  }
  
  
  ##### recipient
  chr_index <- which(recipient_missense_stats$CHROM == paste0("chr", aml_locations$CHROM[id]))
  POS_index <- which(recipient_missense_stats$POS[chr_index] == POS)
  if(length(POS_index) > 0){ 
    
    recipient_aml_mutations$CHROM[id] <- paste0("chr", aml_locations$CHROM[id])
    recipient_aml_mutations$POS[id] <- POS
    if(recipient_missense_stats$REF[chr_index[POS_index]] == aml_locations$REF[id]){
      
      switch (aml_locations$ALT,
              "A" = {recipient_aml_mutations$MutNum[id] <- recipient_missense_stats$ALT.A[chr_index[POS_index]] },
              "T" = {recipient_aml_mutations$MutNum[id] <- recipient_missense_stats$ALT.T[chr_index[POS_index]] },
              "G" = {recipient_aml_mutations$MutNum[id] <- recipient_missense_stats$ALT.G[chr_index[POS_index]] },
              "C" = {recipient_aml_mutations$MutNum[id] <- recipient_missense_stats$ALT.C[chr_index[POS_index]] }
      )
      recipient_aml_mutations$TotalNum[id] <- recipient_missense_stats$NumDiff[chr_index[POS_index]]
      
    }else print(paste0("recipient - ", id))
    
  }else{
    recipient_aml_mutations$CHROM[id] <- paste0("chr", aml_locations$CHROM[id])
    recipient_aml_mutations$POS[id] <- POS
    recipient_aml_mutations$Num[id] <- 0
  }
}
