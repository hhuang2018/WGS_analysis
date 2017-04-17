
renamed_vcf_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_renamed/"

all_vcf_files <- list.files(renamed_vcf_dir, pattern = "\\.vcf.gz$")

output_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/chr6_MHC_ARS/"

# dir.create(file.path(output_dir), showWarnings = FALSE)

num_files <- length(all_vcf_files)

HLA_A_E2 <- "chr6:29942757-29943026"
HLA_A_E3 <- "chr6:29943268-29943543"
HLA_B_E2 <- "chr6:31271599-31271868"
HLA_B_E3 <- "chr6:31271073-31271348"
HLA_C_E2 <- "chr6:31356688-31356957"
HLA_C_E3 <- "chr6:31356167-31356442"
HLA_DRB1_E2 <- "chr6:32584109-32584378"
HLA_DQB1_E2 <- "chr6:32584109-32584378"

for(id in 1:num_files){
  ptm <- proc.time()
  
  cat("VCF file #", id, "\n")
  out.filename <- gsub(".vcf.gz", "_A_E2.vcf", all_vcf_files[id])
  system(paste0("tabix ", renamed_vcf_dir, all_vcf_files[id], " ", HLA_A_E2, " > ", output_dir, out.filename))
  cat(paste0("tabix ", paired_vcf_dir, all_vcf_files[id], " ", HLA_A_E2, " > ", output_dir, out.filename), "\n")
  
  out.filename <- gsub(".vcf.gz", "_A_E3.vcf", all_vcf_files[id])
  system(paste0("tabix ", renamed_vcf_dir, all_vcf_files[id]," ", HLA_A_E3, " > ", output_dir, out.filename))
  cat(paste0("tabix ", paired_vcf_dir, all_vcf_files[id], " ", HLA_A_E3, " > ", output_dir, out.filename), "\n")
  
  out.filename <- gsub(".vcf.gz", "_B_E2.vcf", all_vcf_files[id])
  system(paste0("tabix ", renamed_vcf_dir, all_vcf_files[id]," ", HLA_B_E2, " > ", output_dir, out.filename))
  cat(paste0("tabix ", paired_vcf_dir, all_vcf_files[id], " ", HLA_B_E2, " > ", output_dir, out.filename), "\n")
  
  out.filename <- gsub(".vcf.gz", "_B_E3.vcf", all_vcf_files[id])
  system(paste0("tabix ", renamed_vcf_dir, all_vcf_files[id], " ", HLA_B_E3, " > ", output_dir, out.filename))
  cat(paste0("tabix ", paired_vcf_dir, all_vcf_files[id], " ", HLA_B_E3, " > ", output_dir, out.filename), "\n")
  
  out.filename <- gsub(".vcf.gz", "_C_E2.vcf", all_vcf_files[id])
  system(paste0("tabix ", renamed_vcf_dir, all_vcf_files[id], " ", HLA_C_E2, " > ", output_dir, out.filename))
  cat(paste0("tabix ", paired_vcf_dir, all_vcf_files[id], " ", HLA_C_E2, " > ", output_dir, out.filename), "\n")
  
  out.filename <- gsub(".vcf.gz", "_C_E3.vcf", all_vcf_files[id])
  system(paste0("tabix ", renamed_vcf_dir, all_vcf_files[id], " ", HLA_C_E3, " > ", output_dir, out.filename))
  cat(paste0("tabix ", paired_vcf_dir, all_vcf_files[id], " ", HLA_C_E3, " > ", output_dir, out.filename), "\n")
  
  out.filename <- gsub(".vcf.gz", "_DRB1_E2.vcf", all_vcf_files[id])
  system(paste0("tabix ", renamed_vcf_dir, all_vcf_files[id], " ", HLA_DRB1_E2, " > ", output_dir, out.filename))
  cat(paste0("tabix ", paired_vcf_dir, all_vcf_files[id], " ", HLA_DRB1_E2, " > ", output_dir, out.filename), "\n")
  
  out.filename <- gsub(".vcf.gz", "_DQB1_E2.vcf", all_vcf_files[id])
  system(paste0("tabix ", renamed_vcf_dir, all_vcf_files[id], " ", HLA_DQB1_E2, " > ", output_dir, out.filename))
  cat(paste0("tabix ", paired_vcf_dir, all_vcf_files[id], " ", HLA_DQB1_E2, " > ", output_dir, out.filename), "\n")
  
  print(proc.time()-ptm)
}
cat("Extraction done! \n")

######
load("../Data/HLI_available_pairs_dis_table.RData")
# all_vcf_files <- list.files(output_dir, pattern = "\\.vcf$")

# avail_vcf_ids <- which(sapply(1:length(all_vcf_files), function(x) is.element(unlist(strsplit(all_vcf_files[x], "_"))[2], Available_paired_table$GroupID)))

mismatch_table <- data.frame(GroupID = character(205),
                             A_E2_total = numeric(205),
                             A_E2_lowScore = numeric(205),
                             A_E3_total = numeric(205),
                             A_E3_lowScore = numeric(205),
                             B_E2_total = numeric(205),
                             B_E2_lowScore = numeric(205),
                             B_E3_total = numeric(205),
                             B_E3_lowScore = numeric(205),
                             C_E2_total = numeric(205),
                             C_E2_lowScore = numeric(205),
                             C_E3_total = numeric(205),
                             C_E3_lowScore = numeric(205),
                             DRB1_E2_total = numeric(205),
                             DRB1_E2_lowScore = numeric(205),
                             DQB1_E2_total = numeric(205),
                             DQB1_E2_lowScore = numeric(205),
                             stringsAsFactors = F)
regions <- c("A_E2", "A_E3",
             "B_E2", "B_E3",
             "C_E2", "C_E3", 
             "DRB1_E2", "DQB1_E2")

for(id in 1:205){
  
  prefix <- paste0(Available_paired_table$GroupType[id], "_", Available_paired_table$GroupID[id])
  
  mismatch_table[id, 1] <- Available_paired_table$GroupID[id]
    
  for(jd in 1:8){
    Exon_donor <- read.table(paste0(output_dir, prefix, "_D_", regions[jd], ".vcf"))
    Exon_recipient <- read.table(paste0(output_dir, prefix, "_R_", regions[jd], ".vcf"))
    
    intersect_pos <- intersect(Exon_donor$V2, Exon_recipient$V2)
    
    donor_id <- which(Exon_donor$V2 %in% intersect_pos)
    recipient_id <- which(Exon_recipient$V2 %in% intersect_pos)
    
    mismatch_table[id, 2*jd] <- length(Exon_donor$V2) - length(donor_id) + 
      length(Exon_recipient$V2) - length(recipient_id)  # insertion
    
    same_pos_ALT <- intersect(Exon_donor$V5[donor_id], Exon_recipient$V5[recipient_id])
    
    donor_ALT_id <- which(Exon_donor$V5[donor_id] %in% same_pos_ALT)
    recipient_ALT_id <- which(Exon_recipient$V5[recipient_id] %in% same_pos_ALT)
    
    mismatch_table[id, 2*jd] <-  mismatch_table[id, 2*jd] + 
      length(donor_id) - length(donor_ALT_id) + 
      length(recipient_id) - length(recipient_ALT_id)  # different ALT
    
    mismatch_table[id, 2*jd+1] <- length(which(c(as.character(Exon_donor$V7), as.character(Exon_recipient$V7))!= "PASS"))
  }
}

save(mismatch_table, file = paste0(output_dir, "mimatch_table_summary.RData"))

