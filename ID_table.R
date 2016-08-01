library(vcfR)

vcf_file_dir <- "/mnt/scratch/hhuang/hli_vcf_renamed/"

all_vcf_files <- list.files(vcf_file_dir, pattern = "\\.vcf.gz$")

num_files <- length(all_vcf_files)

ID_table <- data.frame(Group = character(num_files),
                       GroupID = numeric(num_files),
                       subjectType = character(num_files),
                       R_D_ID = numeric(num_files),
                       SeqID = numeric(num_files),
                       stringsAsFactors = FALSE)

for(id in 1:num_files){
  
  file_info <- unlist(strsplit(all_vcf_files[id], "_"))
  ID_table[id, 1:3] <- file_info[1:3]
  ID_table$R_D_ID[id] <- as.numeric(unlist(strsplit(file_info[4], "\\."))[1])
  
  vcf_file <- paste0(vcf_file_dir, all_vcf_files[id])
  vcf_info <- read.vcfR(vcf_file, verbose = FALSE)
  ID_table$SeqID[id] <- as.numeric(colnames(vcf_info@gt)[2])

}

save(ID_table, file = "../Output/ID_table.RData")
write.csv(ID_table, file = "../Output/ID_table.csv", row.names = FALSE)

