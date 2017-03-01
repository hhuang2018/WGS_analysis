GWASH_ID <- read.csv("../ClinVar/GWASH/IDs.csv")

load("../Data/ID_table.RData")

available_IDs <-ID_table[ID_table$GroupID %in% ID_table$GroupID[duplicated(ID_table$GroupID)],]

aa <- intersect(GWASH_ID$did, ID_table$R_D_ID)
bb <- intersect(GWASH_ID$rid, ID_table$R_D_ID)

cc <- intersect(GWASH_ID$did, available_IDs$R_D_ID)
dd <- intersect(GWASH_ID$rid, available_IDs$R_D_ID)

ee <- GWASH_ID[which(GWASH_ID$did %in% aa), ]
ff <- GWASH_ID[which(GWASH_ID$rid %in% bb), ]

# which(ID_table$R_D_ID == 38459665)

write.csv(ee, file = "../CLinVar/GWASH/GWASH_HLI_overlapping_IDs.csv")

# 
HLI_ID <- read.csv("../HLI_ID_pairs.csv")
HLI_ID_pairs <- data.frame(rid = numeric(251),
                           did = numeric(251),
                           bmt = numeric(251), stringsAsFactors = F)

HLI_ID_pairs$rid <- sapply(HLI_ID$RID, function(x) as.numeric(gsub("-", "", x)))
HLI_ID_pairs$did <- sapply(HLI_ID$DID, function(x) as.numeric(gsub("-", "", x)))
HLI_ID_pairs$bmt <- sapply(HLI_ID$BMT, function(x) as.numeric(gsub("-", "", x)))

gg <- intersect(GWASH_ID$did, HLI_ID_pairs$did)
hh <- intersect(GWASH_ID$rid, HLI_ID_pairs$rid)
ii <- intersect(GWASH_ID$bmt_case, HLI_ID_pairs$bmt)

jj <- GWASH_ID[which(GWASH_ID$did %in% gg), ]
kk <- HLI_ID_pairs[which(HLI_ID_pairs$rid %in% hh), ]

###### disease Type
disease_types <- read.csv("../ClinVar/HLI_disease_type.csv")

mm <- intersect(disease_types$rid, available_IDs$R_D_ID)
nn <- intersect(disease_types$did, available_IDs$R_D_ID)

available_groups <- disease_types[which(disease_types$rid %in% mm),]
available_groups2 <- disease_types[which(disease_types$did %in% nn), ]
table(available_groups$disease)
table(available_groups2$disease)

write.csv(available_groups, file = "../ClinVar/HLI_DiseaseTypes.csv")

avail_donors <- available_IDs[which(available_IDs$subjectType == "D"), ]
missing_IDs <- avail_donors[-which(avail_donors$R_D_ID %in% nn), ]

available_groups$Group <- sapply(1:dim(available_groups)[1], function(x) available_IDs$Group[which(available_IDs$R_D_ID %in% available_groups$rid[x])])
available_groups$GroupID <- sapply(1:dim(available_groups)[1], function(x) available_IDs$GroupID[which(available_IDs$R_D_ID %in% available_groups$rid[x])])

available_groups$Group2 <- sapply(1:dim(available_groups)[1], function(x) available_IDs$Group[which(available_IDs$R_D_ID %in% available_groups$did[x])])
available_groups$GroupID2 <- sapply(1:dim(available_groups)[1], function(x) available_IDs$GroupID[which(available_IDs$R_D_ID %in% available_groups$did[x])])

write.csv(available_groups, file = "../ClinVar/HLI_DiseaseTypes_2.csv")


#########
# Overlapping SNPs
#########
GWASH_snps <- read.delim(file = "../ClinVar/GWASH/labcorp0814.map", header = FALSE)
colnames(GWASH_snps) <- c("CHROM", "SNP", "GeneticPos(M)", "POS")

Known_MiHA_SNPs <- read.delim("../WW_MiHA/Known_MiHA_coordinates.txt")

overlapping_SNPs <- intersect(GWASH_snps$SNP, Known_MiHA_SNPs$SNPs)
overlapping_MiHAs_SNPs <- Known_MiHA_SNPs[which(Known_MiHA_SNPs$SNPs %in% c(overlapping_SNPs, "rs10004")), ]
missing_MiHA_snps <- Known_MiHA_SNPs[-which(Known_MiHA_SNPs$SNPs %in% c(overlapping_SNPs, "rs10004")), ]

write.csv(overlapping_MiHAs_SNPs, file = "../ClinVar/GWASH_overlapping_SNPs.csv")
