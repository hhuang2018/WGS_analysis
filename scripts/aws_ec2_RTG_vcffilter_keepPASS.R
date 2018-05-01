#library(doParallel)
#require(foreach)
# paralle setting
#cl <- makeCluster(2)
#registerDoParallel(cl)

# list .vcf.gz files
original_file_dir = "~/data/HLI/"
destination_dir = "~/data/HLI_filtered/"
#reference_dir = "~/data/references/hg38_ucsc.sdf"
#combined_dir = "~/data/HLI_output/combined_cases/"

filenames <- list.files(original_file_dir, pattern = "\\.vcf.gz$")
groupnames <- sapply(1:length(filenames), function(x) paste0(unlist(strsplit(filenames[x], "_"))[1:2], collapse ="_"))

unique_groups <- unique(groupnames)

num_files <- length(unique_groups)

##### need csv format of summary
#foreach(id = 1:num_files) %dopar% {
system(paste0("echo \"Starting Time: ", Sys.time(), "\" > ", destination_dir, "README"), 
       intern = TRUE)

for(id in 1:num_files){
  index <- which(grepl(paste0(unique_groups[id], "_"), filenames))
  
  if(length(index) == 2){
    
    donor_file <- filenames[index[grepl('_D_', filenames[index])]]
    recipient_file <- filenames[index[grepl('_R_', filenames[index])]]
    
    donor_file_filterd <- gsub(".vcf.gz", "_filtered.vcf.gz", donor_file)
    recipient_file_filtered <- gsub(".vcf.gz", "_filtered.vcf.gz", recipient_file)
    
    donor_log <- gsub(".vcf.gz", ".out", donor_file)
    recipient_log <- gsub(".vcf.gz", ".out", recipient_file)
    #cat(paste0(unique_groups[id], ": \n"))
    #cat(paste0(donor_file, " <-> ", recipient_file, "\n"))
    #  RTG_MEM=16G
    system(paste0("echo \"rtg vcffilter -i ", original_file_dir, donor_file, " -o ", destination_dir,
                  donor_file_filterd, " -k PASS > ", destination_dir, donor_log, "\" >> ", destination_dir,"README"), 
           intern = TRUE)
    system(paste0("rtg vcffilter -i ", original_file_dir, donor_file, " -o ", destination_dir,
                  donor_file_filterd, " -k PASS > ", destination_dir ,donor_log), 
           intern = TRUE)
    
    system(paste0("echo \"rtg vcffilter -i ", original_file_dir, recipient_file, " -o ", destination_dir,
                  recipient_file_filtered, " -k PASS > ", destination_dir, recipient_log, "\" >> ", destination_dir, "README"), 
           intern = TRUE)
    system(paste0("rtg vcffilter -i ", original_file_dir, recipient_file, " -o ", destination_dir,
                  recipient_file_filtered, " -k PASS > ", destination_dir, recipient_log), 
           intern = TRUE)

  }
  
}

system(paste0("echo \"Finishing Time: ", Sys.time(), "\" >> ", destination_dir, "README"), 
       intern = TRUE)

## sanity check
#getDoParWorkers() # number of workers doing parallel for-loop
#getDoParName() #  the name and version of the currently registered backend
#getDoParVersion()

#stopCluster(cl)

