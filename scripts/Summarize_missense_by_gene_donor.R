
output_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants_RefSeq_canonical/stats/"

# load("../Output/Missense_variant_stats/RefSeq_canon_0228/Recipient_missesense_stats_RefSeq_canon_0228.RData")
# recipient_missense_stats <- donor_missense_stats
load(paste0(output_dir, "donor_missesense_stats_RefSeq_canon_0228.RData"))

# recipient_missense_stats$NumDiff <- -recipient_missense_stats$NumDiff
# recipient_missense_stats$PertDiff <- recipient_missense_stats$NumDiff/240
donor_missense_stats$PertDiff <- donor_missense_stats$NumDiff/216

# known gene list
load("../Data/hg38_known_gene_symbols_HGNC.RData")

num_gene <- dim(GRCh38_known_gene_list)[1]
donor_gene_stats <- data.frame(GeneName = GRCh38_known_gene_list$name,
                               GeneTxLen = GRCh38_known_gene_list$txEnd - GRCh38_known_gene_list$txStart,
                               GeneMisNum = numeric(num_gene),
                               MisPercent = numeric(num_gene),
                               stringsAsFactors = F)

temp_mat <- donor_missense_stats
counter <- 0
for(chr_id in c(1:22, "X", "Y")){
  
  CH_gene_list <- GRCh38_known_gene_list[which(GRCh38_known_gene_list$Chromosome %in% paste0("chr", chr_id)),]
  temp_mat_list <- temp_mat[which(temp_mat$CHROM %in% paste0("chr", chr_id)), ]
  num_CH_gene <- dim(CH_gene_list)[1]
  
  for(id in 1:num_CH_gene){
    counter <- counter + 1
    gene_index <- which(findInterval(temp_mat_list$POS, c(CH_gene_list$txStart[id], CH_gene_list$txEnd[id]))==1)
    
    if(length(gene_index)>0){
      
      donor_gene_stats$GeneMisNum[counter] <- length(gene_index)
      donor_gene_stats$MisPercent[counter] <- length(gene_index)/donor_gene_stats$GeneTxLen
      
    }
    
    
  }
  
}

for(id in 1:num_gene){
  
  g_index <-  which(donor_gene_stats$GeneName[id] %in% GRCh38_known_gene_list$name)
  donor_gene_stats$GeneSymbol[id] <- GRCh38_known_gene_list$genesymbol[g_index]
  
}

save(donor_gene_stats, file = paste0(output_dir, "Donor_missesense_stats_RefSeq_canon_0228_byGene.RData"))
