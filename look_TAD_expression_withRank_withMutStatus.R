########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript look_TAD_expression_withRank_withMutStatus.R\n"))

script_name <- "look_TAD_expression_withRank_withMutStatus.R"

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
require(ggsci)
ggsci_pal <- "d3"
ggsci_subpal <- ""
require(reshape2)
log10_offset <- 0.01


# Rscript look_TAD_expression_withRank_withMutStatus.R <hicds> <exprds> <gene_symbol>
# Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr11_TAD390 # MMPP
# Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr10_TAD268 # SFTPA


# ENCSR489OCU_NCI-H460_40kb	TCGAluad_mutKRAS_mutEGFR	chr10_TAD16	AKR1C1,AKR1C2,AKR1C3
# ENCSR489OCU_NCI-H460_40kb	TCGAluad_mutKRAS_mutEGFR	chr17_TAD162	HOXB2,HOXB3,HOXB4,HOXB5,HOXB6,HOXB7
# ENCSR489OCU_NCI-H460_40kb	TCGAluad_nonsmoker_smoker	chr10_TAD16	AKR1C1,AKR1C2,AKR1C3
# ENCSR489OCU_NCI-H460_40kb	TCGAluad_nonsmoker_smoker	chr17_TAD162	HOXB2,HOXB3,HOXB4,HOXB5,HOXB6,HOXB7
# 
# ENCSR489OCU_NCI-H460_40kb	TCGAluad_norm_luad	chr11_TAD390	MMP1,MMP12,MMP13
# ENCSR489OCU_NCI-H460_40kb	TCGAluad_norm_luad	chr10_TAD268	BEND3P3,SFTPA1,SFTPA2
# ENCSR489OCU_NCI-H460_40kb	TCGAlusc_norm_lusc	chr11_TAD390	MMP1,MMP10,MMP12,MMP13,MMP3
# ENCSR489OCU_NCI-H460_40kb	TCGAlusc_norm_lusc	chr10_TAD268	BEND3P3,SFTPA1,SFTPA2


# Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD16
# Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr17_TAD162
# Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker chr10_TAD16
# Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker chr17_TAD162
# Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD390
# Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr10_TAD268
# Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr11_TAD390
# Rscript look_TAD_expression_withRank_withMutStatus.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr10_TAD268
# 




hicds="ENCSR489OCU_NCI-H460_40kb"
exprds="TCGAlusc_norm_lusc"
tad_to_plot="chr11_TAD390"

col1 <- pal_futurama()(5)[1]
col2 <- pal_futurama()(5)[5]
col1 <- pal_aaas()(5)[4]
col2 <- pal_npg()(5)[5]

mutCol <- "chartreuse3"
notmutCol <- "darkgrey"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
hicds <- args[1]
exprds <- args[2]
tad_to_plot <- args[3]

cat("load inDT \n")
inDT <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
inDT <- inDT[inDT$hicds == hicds & inDT$exprds == exprds,]

tad_plot_rank <- unique(inDT$tad_rank[inDT$hicds == hicds & inDT$exprds == exprds & inDT$region == tad_to_plot])
stopifnot(!is.na(tad_plot_rank))
stopifnot(length(tad_plot_rank) == 1)

plotSub <- paste0(tad_to_plot, " - rank: ", tad_plot_rank)

my_xlab <- "TAD genes (ordered by start positions)"
my_ylab <- "RNA-seq expression count [log10]"

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4
myHeightGG <- 6
myWidthGG <- 7.5

source("../FIGURES_V2_YUANLONG/settings.R")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

mainFolder <- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA")
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path(mainFolder, "PIPELINE", "INPUT_FILES")

outFolder <- file.path("LOOK_TAD_EXPRESSION_WITHRANK_WITHMUTSTATUS")
dir.create(outFolder, recursive = TRUE)

mutSamples <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/NFE2L2_KEAP1_MUTSAMPLES/mut_samples.Rdata"))

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
gff_dt <- gff_dt[order(gff_dt$chromo, gff_dt$start, gff_dt$end, gff_dt$entrezID),]

g2t_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
g2t_dt <- read.delim(g2t_file, col.names=c("entrezID", "chromo", "start", "end", "region"), header=FALSE, stringsAsFactors = FALSE)

stopifnot(tad_to_plot %in% g2t_dt$region)

settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
stopifnot(file.exists(settingFile))
source(settingFile)

samp1 <- get(load(file.path(setDir, sample1_file)))
samp2 <- get(load(file.path(setDir, sample2_file)))

geneList <- get(load(file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
stopifnot(geneList %in% g2t_dt$entrezID)

toplot_dt <- g2t_dt[g2t_dt$region == tad_to_plot & g2t_dt$entrezID %in% geneList,]
toplot_dt <- toplot_dt[order(toplot_dt$start, toplot_dt$end),]

stopifnot(toplot_dt$entrezID %in% geneList)
plot_geneList <- geneList[geneList %in% toplot_dt$entrezID]

count_file <- file.path(pipFolder, hicds, exprds, "1_runGeneDE", "DE_rnaseqDT.Rdata")
stopifnot(file.exists(count_file))
fpkm_dt <- get(load(count_file))
stopifnot(names(geneList) %in% rownames(fpkm_dt))

fpkm_plot_dt <- fpkm_dt[rownames(fpkm_dt) %in% names(plot_geneList),]
fpkm_plot_dt$entrezID <- plot_geneList[paste0(rownames(fpkm_plot_dt))]
stopifnot(!is.na(fpkm_plot_dt$entrez))

stopifnot(samp1 %in% colnames(fpkm_plot_dt))
stopifnot(samp2 %in% colnames(fpkm_plot_dt))

fpkm_plot_dt <- fpkm_plot_dt[,c("entrezID", samp1, samp2)]
m_fpkm_dt <- melt(fpkm_plot_dt, id="entrezID")

stopifnot(m_fpkm_dt$entrezID %in% toplot_dt$entrezID)

toplot_dt <- merge(m_fpkm_dt, toplot_dt, merge="entrezID", all.x=TRUE, all.y=FALSE)
stopifnot(!is.na(toplot_dt))
toplot_dt$cond <- ifelse(toplot_dt$variable %in% samp1, cond1, 
                         ifelse(toplot_dt$variable %in% samp2, cond2, NA ))
stopifnot(!is.na(toplot_dt$cond))
toplot_dt$symbol <- entrez2symb[paste0(toplot_dt$entrezID)]
stopifnot(!is.na(toplot_dt$symbol))

toplot_dt <- toplot_dt[order(toplot_dt$chromo, toplot_dt$start, toplot_dt$end, toplot_dt$value),]
toplot_dt$value_log10 <- log10(toplot_dt$value + log10_offset)

# withRank_toplot_dt <- do.call(rbind, by(toplot_dt, list(toplot_dt$symbol, toplot_dt$cond), function(x) {
#   dt <- x[order(x$value, decreasing = TRUE),]
#   dt$samp_rank <- 1:nrow(dt)
#   dt
# }))


withRank_toplot_dt2 <- do.call(rbind, by(toplot_dt, list(toplot_dt$symbol), function(x) {
  x$cond <- factor(x$cond, levels=c(cond1,cond2))
  dt <- x[order(-as.numeric(x$cond), x$value, decreasing = TRUE),]
  dt$samp_rank <- 1:nrow(dt)
  dt
}))
withRank_toplot_dt2$hicds <- hicds
withRank_toplot_dt2$exprds <- exprds

cat("merge withRank and inDT \n")

withRank_toplot_dt2 <- merge(withRank_toplot_dt2, inDT, all.x=TRUE, all.y=FALSE, by=intersect(colnames(inDT), colnames(withRank_toplot_dt2)))


tmp <- withRank_toplot_dt2[,c("symbol", "start", "end", "gene_rank")]
tmp <- unique(tmp)
tmp <- tmp[order(tmp$start, tmp$end),]
tmp$lab <- paste0(tmp$symbol, "\n(rank: ", tmp$gene_rank, ")")

withRank_toplot_dt2$symbol <- factor(withRank_toplot_dt2$symbol, levels=tmp$symbol)

withRank_toplot_dt2$cond_dots <- withRank_toplot_dt2$cond
withRank_toplot_dt2$cond_dots <- ifelse(as.character(withRank_toplot_dt2$variable) %in% mutSamples, "withMut", as.character(withRank_toplot_dt2$cond_dots))
withRank_toplot_dt2$cond_dots <- factor(withRank_toplot_dt2$cond_dots, levels = c(cond1,cond2, "withMut"))
stopifnot(!is.na(withRank_toplot_dt2$cond_dots))

withRank_toplot_dt2$cond_sh <- withRank_toplot_dt2$cond
withRank_toplot_dt2$cond_sh <- ifelse(as.character(withRank_toplot_dt2$variable) %in% mutSamples, "withMut", "noMut")
withRank_toplot_dt2$cond_sh <- factor(withRank_toplot_dt2$cond_sh, levels = c("noMut", "withMut"))
stopifnot(!is.na(withRank_toplot_dt2$cond_sh))


withRank_toplot_dt2$cond <- factor(withRank_toplot_dt2$cond, levels = c(cond1,cond2))
stopifnot(!is.na(withRank_toplot_dt2$cond))

withRank_toplot_dt2$symbol_lab <- paste0(withRank_toplot_dt2$symbol, "\n(rank: ", withRank_toplot_dt2$gene_rank, ")")
withRank_toplot_dt2$symbol_lab <- factor(withRank_toplot_dt2$symbol_lab, levels=tmp$lab)
stopifnot(!is.na(withRank_toplot_dt2$symbol_lab))

save(withRank_toplot_dt2, file ="withRank_toplot_dt2.Rdata")

subTit <- paste0(tad_to_plot, " (rank: ", tad_plot_rank, ")")

p_var_boxplot <-  ggplot(withRank_toplot_dt2, aes(x = symbol_lab, y = value_log10, fill = cond)) + 
  
  geom_boxplot(notch = TRUE, outlier.shape=NA)+
  geom_jitter(aes(colour = cond, fill =cond, shape=cond_sh), position=position_jitterdodge(), alpha=0.7)+
  # geom_jitter(aes(colour = cond, shape=cond_sh), position=position_jitterdodge())+
  
  ggtitle(paste0(hicds, " - ", exprds), subtitle = paste0(subTit))+
  scale_x_discrete(name=my_xlab)+
  scale_y_continuous(name=paste0(my_ylab),
                     breaks = scales::pretty_breaks(n = 20))+
  
  scale_shape_manual(
    values = c(15,8),
    breaks = c("noMut", "withMut"),
    labels = c("not mut.", "mut.")
  )+
  
  scale_color_manual(values=c(col1, col2))+
  scale_fill_manual(values=c(col1, col2))+
  
  labs(fill  = paste0("Cond."), color=paste0("Cond."), shape=paste0("KEAP1|NFE2L2")) +
  theme( 
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    axis.line.x= element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
    axis.text.x =element_text(color="black", hjust=0.5,vjust = 0.5, size=12, face="bold"),
    # axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.text = element_text(size=12),
    legend.key = element_blank(),
    legend.title = element_text(face="bold", size=12)
  )

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_plot, "_allSamples_exprValues_boxplot_vShape.", plotType))
ggsave(plot = p_var_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))


withRank_toplot_dt2$cond_border <- paste0(withRank_toplot_dt2$cond , "_", withRank_toplot_dt2$cond_sh)
withRank_toplot_dt2$cond_border <- factor(withRank_toplot_dt2$cond_border, levels = c(
  paste0(cond1, "_noMut"),  paste0(cond1, "_withMut"), paste0(cond2, "_noMut"), paste0(cond2, "_withMut")
))
stopifnot(!is.na(withRank_toplot_dt2$cond_border))


save(withRank_toplot_dt2, file="withRank_toplot_dt2.Rdata", version=2)

p_var_boxplot <- ggplot(withRank_toplot_dt2, aes(x = symbol_lab, y = value_log10, fill = cond)) +
  geom_boxplot(notch = TRUE, outlier.shape=NA)+
  geom_point(aes(color=cond_border), position=position_jitterdodge(), stroke=0.8, shape=21, alpha=0.8) +
  ggtitle(paste0(hicds, " - ", exprds), subtitle = paste0(subTit))+
  scale_x_discrete(name=my_xlab)+
  scale_y_continuous(name=paste0(my_ylab),
                     breaks = scales::pretty_breaks(n = 20))+
  
  ### NEED TO USE setNames => TO HAVE THE MATCH WHEN FOR ONE CATEGORY I HAVE NO DATA !!!
  scale_color_manual( values=c( setNames(c(col1, mutCol,  col2, mutCol),   c(paste0(cond1, "_noMut"),  paste0(cond1, "_withMut"), paste0(cond2, "_noMut"), paste0(cond2, "_withMut")))),
                      labels = c(setNames(c(paste0(cond1," -\nnot mut."), paste0(cond1," -\nmut."), paste0(cond2," -\nnot mut."), paste0(cond2," -\nmut.")), 
                                          c(paste0(cond1, "_noMut"),  paste0(cond1, "_withMut"), paste0(cond2, "_noMut"), paste0(cond2, "_withMut")))))+
  
  scale_fill_manual( values=c(col1, col2))+
  labs(fill  = paste0("Cond."), fill=paste0("Cond."), color=paste0("KEAP1|NEF2L2")) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    axis.line.x= element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
    axis.text.x =element_text(color="black", hjust=0.5,vjust = 0.5, size=12, face="bold"),
    # axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.text = element_text(size=12),
    legend.key = element_blank(),
    
    legend.key.size = unit(1.2, 'cm'),
    
    legend.title = element_text(face="bold", size=12)
  )

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_plot, "_allSamples_exprValues_boxplot_vCol.", plotType))
ggsave(plot = p_var_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))





save(withRank_toplot_dt2, file="withRank_toplot_dt2.Rdata", version=2)

p_var_boxplot <- ggplot(withRank_toplot_dt2, aes(x = symbol_lab, y = value_log10, fill = cond)) +
  
  geom_boxplot(notch = TRUE, outlier.shape=NA)+
  geom_point(aes(color=cond_sh, shape=cond), position=position_jitterdodge(), stroke=0.8, size=1,alpha=0.8) +
  ggtitle(paste0(hicds, " - ", exprds), subtitle = paste0(subTit))+
  
  scale_x_discrete(name=my_xlab)+
  scale_y_continuous(name=paste0(my_ylab),
                     breaks = scales::pretty_breaks(n = 20))+
  
  ### NEED TO USE setNames => TO HAVE THE MATCH WHEN FOR ONE CATEGORY I HAVE NO DATA !!!
  scale_color_manual( values=c( setNames(c(notmutCol, mutCol),   c(paste0("noMut"), paste0( "withMut")))),
                      labels = c(setNames(c(paste0("not mut."), paste0("mut.")), 
                                          c( paste0("noMut"), paste0( "withMut")))))+
  
  scale_shape_manual(values=c(1,2))+
  
  scale_fill_manual( values=c(col1, col2))+
  labs(shape  = paste0("Cond."), fill=paste0("Cond."), color=paste0("KEAP1|NEF2L2")) +
  
  
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    axis.line.x= element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
    axis.text.x =element_text(color="black", hjust=0.5,vjust = 0.5, size=12, face="bold"),
    # axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.text = element_text(size=12),
    legend.key = element_blank(),
    
    legend.key.size = unit(1.2, 'cm'),
    
    legend.title = element_text(face="bold", size=12)
  )

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_plot, "_allSamples_exprValues_boxplot_vShapeCol.", plotType))
ggsave(plot = p_var_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))








##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))

