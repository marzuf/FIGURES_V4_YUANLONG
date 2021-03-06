########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript look_fam_expression_withRank_withMutStatus_facet.R\n"))

script_name <- "look_fam_expression_withRank_withMutStatus_facet.R"

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

widthExp <- 1.8

xtxtSize <- 10

# Rscript look_fam_expression_withRank_withMutStatus_facet.R <hicds> <exprds> <gene_symbol>

# Rscript look_fam_expression_withRank_withMutStatus_facet.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad MMP
# Rscript look_fam_expression_withRank_withMutStatus_facet.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad SFTPA
# Rscript look_fam_expression_withRank_withMutStatus_facet.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR HOXB
# Rscript look_fam_expression_withRank_withMutStatus_facet.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR AKR1C

hicds="ENCSR489OCU_NCI-H460_40kb"
exprds="TCGAlusc_norm_lusc"
fam_to_plot="MMP"

col1 <- pal_futurama()(5)[1]
col2 <- pal_futurama()(5)[5]
col1 <- pal_aaas()(5)[4]
col2 <- pal_npg()(5)[5]
col1_mut <- "#40104e"
col2_mut <- "#b23811"


mutCol <- "chartreuse3"
notmutCol <- "darkgrey"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
hicds <- args[1]
exprds <- args[2]
fam_to_plot <- args[3]

cat("load inDT \n")
inDT <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
inDT <- inDT[inDT$hicds == hicds & inDT$exprds == exprds,]

my_xlab <- paste0(fam_to_plot, "* genes (ordered by start positions)")
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

outFolder <- file.path("LOOK_FAM_EXPRESSION_WITHRANK_WITHMUTSTATUS_FACET")
dir.create(outFolder, recursive = TRUE)

mutSamples <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/NFE2L2_KEAP1_MUTSAMPLES/mut_samples.Rdata"))

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
gff_dt <- gff_dt[order(gff_dt$chromo, gff_dt$start, gff_dt$end, gff_dt$entrezID),]

fam_entrez <- gff_dt$entrezID[grepl(paste0("^", fam_to_plot), gff_dt$symbol)]

g2t_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
g2t_dt <- read.delim(g2t_file, col.names=c("entrezID", "chromo", "start", "end", "region"), header=FALSE, stringsAsFactors = FALSE)


settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
stopifnot(file.exists(settingFile))
source(settingFile)

samp1 <- get(load(file.path(setDir, sample1_file)))
samp2 <- get(load(file.path(setDir, sample2_file)))

geneList <- get(load(file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
stopifnot(geneList %in% g2t_dt$entrezID)

toplot_dt <- g2t_dt[g2t_dt$entrezID %in% geneList & g2t_dt$entrezID %in% fam_entrez,]
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

tmp <- withRank_toplot_dt2[,c("symbol", "start", "end", "gene_rank", "region", "chromo", "tad_rank")]
tmp <- unique(tmp)
tmp <- tmp[order(tmp$chromo, tmp$start, tmp$end),]
tmp$lab <- paste0(tmp$symbol, " (", tmp$gene_rank,")\n", tmp$region, " (", tmp$tad_rank ,")")

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

save(withRank_toplot_dt2, file ="withRank_toplot_dt2.Rdata", version=2)

withRank_toplot_dt2$symbol_lab <- paste0(withRank_toplot_dt2$symbol, " (", withRank_toplot_dt2$gene_rank,")\n",  withRank_toplot_dt2$region, " (", withRank_toplot_dt2$tad_rank, ")")
withRank_toplot_dt2$symbol_lab <- factor(withRank_toplot_dt2$symbol_lab, levels=tmp$lab)
stopifnot(!is.na(withRank_toplot_dt2$symbol_lab))

save(withRank_toplot_dt2, file ="withRank_toplot_dt2.Rdata")

# > unique(withRank_toplot_dt2$cond_dots)
# [1] norm    luad    withMut
# Levels: norm luad withMut
# > unique(withRank_toplot_dt2$cond_sh)
# [1] noMut   withMut
# Levels: noMut withMut
# > unique(withRank_toplot_dt2$cond)
# [1] norm luad


withRank_toplot_dt2$all_cond <- interaction(withRank_toplot_dt2$cond, withRank_toplot_dt2$cond_sh)
withRank_toplot_dt2$all_cond <- factor(withRank_toplot_dt2$all_cond, levels = 
                                         c(paste0(cond1, ".noMut"),paste0(cond1, ".withMut"),  paste0(cond2, ".noMut"), paste0(cond2, ".withMut")))

withRank_toplot_dt2$group <- interaction(withRank_toplot_dt2$symbol_lab, withRank_toplot_dt2$cond)


subTit <- paste0(fam_to_plot, " genes")

p_var_boxplot <-  ggplot(withRank_toplot_dt2, aes(x = symbol_lab, y = value_log10, fill = cond, group=group)) + 
  geom_boxplot(notch = TRUE, outlier.shape=NA)+
  geom_jitter(aes(colour = all_cond), position=position_jitterdodge(), alpha=0.7)+
  # geom_jitter(aes(colour = cond, shape=cond_sh), position=position_jitterdodge())+
  ggtitle(paste0(hicds, " - ", exprds), subtitle = paste0(subTit))+
  scale_x_discrete(name=my_xlab,drop=FALSE)+
  scale_y_continuous(name=paste0(my_ylab),
                     breaks = scales::pretty_breaks(n = 20))+

  
  scale_color_manual(values=setNames(c(col1, col1_mut, col2, col2_mut), 
                                     c(paste0(cond1, ".noMut"),paste0(cond1, ".withMut"),  paste0(cond2, ".noMut"), paste0(cond2, ".withMut"))),
                     labels=setNames(c(paste0(cond1, "-wt"),paste0(cond1, "-mut"),  paste0(cond2, "-wt"), paste0(cond2, "-mut")),
                                     c(paste0(cond1, ".noMut"),paste0(cond1, ".withMut"),  paste0(cond2, ".noMut"), paste0(cond2, ".withMut")))) + 
  
  scale_fill_manual(values=setNames(c(col1, col2),
                                    c(cond1, cond2)))+
  
  
  labs(fill  = paste0(""), color=paste0("Cond.-\nKEAP1|NEF2L2 stat."))+
  theme( 
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    axis.line.x= element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
    axis.text.x =element_text(color="black", hjust=1,vjust = 0.5, size=xtxtSize, face="bold", angle=90),
    # axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.text = element_text(size=12),
    legend.key = element_blank(),
    legend.title = element_text(face="bold", size=12, hjust=0.5)
  )

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", fam_to_plot, "_allSamples_exprValues_boxplot_condBox_condmutJitter.", plotType))
ggsave(plot = p_var_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*widthExp)
cat(paste0("... written: ", outFile, "\n"))


p_var_boxplot <-  ggplot(withRank_toplot_dt2, aes(x = symbol_lab, y = value_log10, fill = all_cond)) + 
  geom_boxplot(notch = TRUE, outlier.shape=NA)+
  geom_jitter(aes(colour = all_cond, fill = all_cond), position=position_jitterdodge(), alpha=0.7)+
  # geom_jitter(aes(colour = cond, shape=cond_sh), position=position_jitterdodge())+
  ggtitle(paste0(hicds, " - ", exprds), subtitle = paste0(subTit))+
  scale_x_discrete(name=my_xlab)+
  scale_y_continuous(name=paste0(my_ylab),
                     breaks = scales::pretty_breaks(n = 20))+

  scale_color_manual(values=setNames(c(col1, col1_mut, col2, col2_mut), 
                                     c(paste0(cond1, ".noMut"),paste0(cond1, ".withMut"),  paste0(cond2, ".noMut"), paste0(cond2, ".withMut"))),
                     labels=setNames(c(paste0(cond1, "-wt"),paste0(cond1, "-mut"),  paste0(cond2, "-wt"), paste0(cond2, "-mut")),
                              c(paste0(cond1, ".noMut"),paste0(cond1, ".withMut"),  paste0(cond2, ".noMut"), paste0(cond2, ".withMut")))) + 
  
  scale_fill_manual(values=setNames(c(col1, col1_mut, col2, col2_mut), 
                                     c(paste0(cond1, ".noMut"),paste0(cond1, ".withMut"),  paste0(cond2, ".noMut"), paste0(cond2, ".withMut")))) + 
  
  guides(fill = FALSE) +
  
  labs( color=paste0("Cond.-\nKEAP1|NEF2L2 stat.")) +
  
  theme( 
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    axis.line.x= element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
    axis.text.x =element_text(color="black", hjust=1,vjust = 0.5, size=xtxtSize, face="bold", angle=90),
    # axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.text = element_text(size=12),
    legend.key = element_blank(),
    legend.title = element_text(face="bold", size=12, hjust=0.5)
  )

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", fam_to_plot, "_allSamples_exprValues_boxplot_condmutBox_condmutJitter.", plotType))
ggsave(plot = p_var_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*widthExp)
cat(paste0("... written: ", outFile, "\n"))


##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))

