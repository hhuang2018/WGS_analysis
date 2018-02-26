# process HBD/IBD result files from BEAGLE
source('util.R', echo = FALSE)

# load("../Data/ID_table.RData")

beagles_output_dir <- "../Output/beagles_output/"

new_IBD_table <- reformat_beagles(beagles_output_file = paste0(beagles_output_dir, "aGVHD_all_pairs.gt.ibd.ibd"))

save(new_IBD_table, file = paste0(beagles_output_dir, "aGVHD_all_pairs_IBD.RData"))
write.csv(new_IBD_table, file = paste0(beagles_output_dir, "aGVHD_all_pairs_IBD.csv"), row.names = FALSE)
cat("IBD table has been converted and saved!")

new_HBD_table <- reformat_beagles(beagles_output_file = paste0(beagles_output_dir, "aGVHD_all_pairs.gt.ibd.hbd"))

save(new_HBD_table, file = paste0(beagles_output_dir, "aGVHD_all_pairs_HBD.RData"))
write.csv(new_IBD_table, file = paste0(beagles_output_dir, "aGVHD_all_pairs_HBD.csv", row.names = FALSE))
cat("HBD table has been converted and saved!")