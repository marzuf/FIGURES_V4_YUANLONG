

mainFolder = "../v2_Yuanlong_Cancer_HiC_data_TAD_DA/"

hicds = "ENCSR489OCU_NCI-H460_40kb"

g2t_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
g2t_dt <- read.delim(g2t_file, col.names=c("entrezID", "chromo", "start", "end", "region"), header=FALSE, stringsAsFactors = FALSE)
g2t_dt$entrezID <- as.character(g2t_dt$entrezID)

setDir="/media/electron"
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))


require(dplyr)

g2t_dt_symb <- merge(g2t_dt[,c("entrezID", "region")], gff_dt[,c("entrezID", "symbol")], by="entrezID", all.x=T, all.y=F)

g2t_dt_symb[grepl("^HOXB", g2t_dt_symb$symbol),]
unique(g2t_dt_symb$region[grepl("^HOXB", g2t_dt_symb$symbol)])
# "chr17_TAD162" "chr17_TAD163"

g2t_dt_symb[grepl("^AKR1C", g2t_dt_symb$symbol),]
unique(g2t_dt_symb$region[grepl("^AKR1C", g2t_dt_symb$symbol)])
# "chr10_TAD17" "chr10_TAD16"

g2t_dt_symb[grepl("^MMP", g2t_dt_symb$symbol),]
unique(g2t_dt_symb$region[grepl("^MMP", g2t_dt_symb$symbol)])
# [1] "chr20_TAD115" "chr10_TAD446" "chr11_TAD390" "chr16_TAD150" "chr11_TAD389" "chr20_TAD152" "chr22_TAD29"  "chr14_TAD15" 
# [9] "chr16_TAD160" "chr8_TAD347"  "chr12_TAD501" "chr12_TAD210" "chr11_TAD18"  "chr16_TAD14"  "chr17_TAD116" "chr1_TAD4"   

g2t_dt_symb[grepl("^SFTP", g2t_dt_symb$symbol),]
unique(g2t_dt_symb$region[grepl("^SFTP", g2t_dt_symb$symbol)])
# "chr10_TAD268" "chr2_TAD331"  "chr8_TAD85"   "chr10_TAD269"

Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr17_TAD162
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr17_TAD163
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD17
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD16

Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr20_TAD115
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr10_TAD446
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD390
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr16_TAD150
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD389
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr20_TAD152
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr22_TAD29
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr14_TAD15
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr16_TAD160
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr8_TAD347
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr12_TAD501
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr12_TAD210
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD18
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr16_TAD14
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr17_TAD116
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr1_TAD4

Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr10_TAD268
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr2_TAD331
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr8_TAD85
Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr10_TAD269


Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr17_TAD162
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr17_TAD163
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD17
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD16

Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr20_TAD115
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr10_TAD446
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD390
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr16_TAD150
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD389
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr20_TAD152
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr22_TAD29
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr14_TAD15
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr16_TAD160
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr8_TAD347
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr12_TAD501
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr12_TAD210
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD18
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr16_TAD14
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr17_TAD116
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr1_TAD4

Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr10_TAD268
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr2_TAD331
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr8_TAD85
Rscript look_TAD_expression_withRank.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr10_TAD269




