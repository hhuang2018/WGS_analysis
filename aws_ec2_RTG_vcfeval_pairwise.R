library(doParallel)
#require(foreach)
# paralle setting
cl <- makeCluster(8)
registerDoParallel(cl)


# list .vcf.gz files
#original_file_dir = "~/data/HLI/"
original_file_dir = "~/data/HLI_filtered/"
destination_dir = "~/data/HLI_output/"
reference_dir = "~/data/references/hg38_ucsc.sdf"
combined_dir = "~/data/HLI_output/combined_cases_filtered/"

filenames <- list.files(original_file_dir, pattern = "\\.vcf.gz$")
groupnames <- sapply(1:length(filenames), function(x) paste0(unlist(strsplit(filenames[x], "_"))[1:2], collapse ="_"))

unique_groups <- unique(groupnames)

num_files <- length(unique_groups)

##### need csv format of summary
foreach(id = 1:num_files) %dopar% {
  
  index <- which(grepl(unique_groups[id], filenames))
  
  if(length(index) == 2){
    
    donor_file <- filenames[index[grepl('_D_', filenames[index])]]
    recipient_file <- filenames[index[grepl('_R_', filenames[index])]]
    
    cat(paste0(unique_groups[id], ": \n"))
    cat(paste0(donor_file, " <-> ", recipient_file, "\n"))
    
    system(paste0("rtg RTG_MEM=4G vcfeval -b ", original_file_dir, donor_file, " -c ",
                  original_file_dir, recipient_file, " -o ", destination_dir,
                  unique_groups[id], " --output-mode combine -t ", reference_dir), 
           intern = TRUE)
    
    system(paste0("cp ", destination_dir, unique_groups[id], "/output.vcf.gz ", combined_dir, unique_groups[id],'_filtered.vcf.gz'))
    system(paste0("cp ", destination_dir, unique_groups[id], "/output.vcf.gz.tbi ", combined_dir, unique_groups[id],'_filtered.vcf.gz.tbi'))
    
  }

}

## sanity check
getDoParWorkers() # number of workers doing parallel for-loop
getDoParName() #  the name and version of the currently registered backend
getDoParVersion()

stopCluster(cl)

