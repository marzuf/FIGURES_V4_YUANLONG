
startTime <- Sys.time()
script_name <- "corrExpr_geneWithTAD_condCol.R"

cat(">START ", script_name, "\n")


# Rscript corrExpr_geneWithTAD_condCol.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD268 DNMT1

require(reshape2)
require(ggplot2)
require(ggsci)
require(foreach)
require(doMC)
registerDoMC(40)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 4) {
  hicds <- args[1] 
  exprds <- args[2] 
  tad_to_plot <- args[3]
  symbol_to_plot <- args[4]
} else {
  hicds <-  "ENCSR489OCU_NCI-H460_40kb"
  exprds <- "TCGAluad_norm_luad"
  tad_to_plot <- "chr10_TAD268"
  symbol_to_plot <- "DNMT1"
}
col1 <- pal_futurama()(5)[1]
col2 <- pal_futurama()(5)[5]
col1 <- pal_aaas()(5)[4]
col2 <- pal_npg()(5)[5]


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
log10_offset <- 0.01


plotType <- "svg"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4
myHeightGG <- 6
myWidthGG <- 7.5

log10_offset <- 0.01

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)

stopifnot(symbol_to_plot %in% gff_dt$symbol)

entrez_to_plot <- gff_dt$entrezID[gff_dt$symbol == symbol_to_plot]
stopifnot(length(entrez_to_plot) == 1)

stopifnot(!duplicated(gff_dt$entrezID))
stopifnot(!duplicated(gff_dt$symbol))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)

outFolder <- "CORREXPR_GENEWITHTAD_CONDCOL"
dir.create(outFolder, recursive = TRUE)

runFolder <- "../v2_Yuanlong_Cancer_HiC_data_TAD_DA"

pipFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path(runFolder, "PIPELINE", "INPUT_FILES")

settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
stopifnot(file.exists(settingFile))
source(settingFile)

g2t_file <- file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt")
g2t_dt <- read.delim(g2t_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2t_dt$entrezID <- as.character(g2t_dt$entrezID)

gene_list <- get(load(file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata") ))


tad_g2t_dt <- g2t_dt[g2t_dt$region == tad_to_plot & g2t_dt$entrezID %in% gene_list,]

stopifnot(tad_g2t_dt$entrezID %in% gene_list)
geneList_to_plot <- gene_list[gene_list %in% tad_g2t_dt$entrezID]

## => change here, take only the expressed TFs
de_dt <- get(load(file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "1_runGeneDE", "DE_topTable.Rdata") ))

stopifnot(names(geneList_to_plot) %in% de_dt$genes)

stopifnot(entrez_to_plot %in% de_dt$genes)

plot_de_dt <- de_dt[de_dt$genes %in% names(geneList_to_plot),]

plot_de_dt$entrezID <- sapply(plot_de_dt$genes, function(x) geneList_to_plot[names(geneList_to_plot) == x ])
stopifnot(!is.na(plot_de_dt$entrezID))
stopifnot(setequal(plot_de_dt$entrezID, geneList_to_plot))

plot_de_dt$symbol <- sapply(plot_de_dt$entrezID, function(x) gff_dt$symbol[gff_dt$entrezID == x])
stopifnot(!is.na(plot_de_dt$symbol))

count_file <- file.path(pipFolder, hicds, exprds, "1_runGeneDE", "DE_rnaseqDT.Rdata")
stopifnot(file.exists(count_file))
fpkm_dt <- get(load(count_file))

stopifnot(plot_de_dt$genes %in% rownames(fpkm_dt))

stopifnot(nrow(plot_de_dt) > 0)

stopifnot(length(entrez_to_plot) == 1)
stopifnot(entrez_to_plot %in% rownames(fpkm_dt))

tf_count <- as.numeric(fpkm_dt[paste0(entrez_to_plot),])

samp1_id <- get(load(file.path(setDir, sample1_file)))
samp2_id <- get(load(file.path(setDir, sample2_file)))

stopifnot(colnames(fpkm_dt) %in% samp1_id | colnames(fpkm_dt) %in% samp2_id)

dotCols <- ifelse(colnames(fpkm_dt) %in% samp1_id, col1, 
                  ifelse(colnames(fpkm_dt) %in% samp2_id, col2, NA))
stopifnot(!is.na(dotCols))

mutSamples <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/NFE2L2_KEAP1_MUTSAMPLES/mut_samples.Rdata"))

shapeMut <- 15
shapeWt <- 8

dotShape <- ifelse(colnames(fpkm_dt) %in% mutSamples, shapeMut, shapeWt)
stopifnot(!is.na(dotShape))


for(i in 1:nrow(plot_de_dt)) {
  
  curr_entrez <- plot_de_dt$genes[i]
  stopifnot(curr_entrez %in% rownames(fpkm_dt))
  gene_count <- as.numeric(fpkm_dt[paste0(curr_entrez),])
  
  curr_symbol <- plot_de_dt$symbol[i]
  
  stopifnot(length(tf_count) == length(gene_count))

  myx <- log10(tf_count+log10_offset)
  myy <- log10(gene_count+log10_offset)
  
  
  
  outFile <- file.path(outFolder, paste0(symbol_to_plot, "_vs_", curr_symbol, "_", tad_to_plot, "_fpkmCount_colplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
  
  
  par(bty="L")
  par(mar=c(5,4,2,6))
  plot(
    x = myx,
    y = myy,
    main = "RNA-seq FPKM [log10]",
    xlab = paste0(symbol_to_plot, " (reg. element)"),
    ylab = paste0(curr_symbol, " (", tad_to_plot, " gene)"),
    cex = 0.7,
    pch = dotShape,
    col= dotCols,
    cex.lab = plotCex,
    cex.axis = plotCex,
    cex.main = plotCex
  )
  mtext(side = 3, text = paste0(hicds, " - ", exprds))
  addCorr(x = myx, y=myy, bty="n", legPos="topleft"  )
  abline(lm(myy~myx), lty=2, col="darkgrey")  
  
  par(xpd=TRUE)
  
  legend(
    # "bottomleft",
    c(max(myx), max(myy)),
    legend= c(cond1, cond2, "KEAP1|NFE2L mut.", "KEAP1|NFE2L wt."),
    pch = c(NA, NA, shapeMut, shapeWt),
    text.col = c(col1, col2, "black", "black"),
    bty="n"
  )
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}
stopifnot(plot_de_dt$genes %in% rownames(fpkm_dt))
tad_count_dt <- fpkm_dt[rownames(fpkm_dt) %in% plot_de_dt$genes,]

tad_mean_count <- as.numeric(colMeans(tad_count_dt))
stopifnot(length(tf_count) == length(tad_mean_count))


myx <- log10(tf_count+log10_offset)
myy <- log10(tad_mean_count+log10_offset)

outFile <- file.path(outFolder, paste0(symbol_to_plot, "_vs_", tad_to_plot, "_mean_fpkmCount_colplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))

par(bty="L")
par(mar=c(5,4,2,6))
plot(
  x = myx,
  y = myy,
  main = "RNA-seq FPKM [log10]",
  xlab = paste0(symbol_to_plot, " (reg. element)"),
  ylab = paste0(tad_to_plot, " mean counts"),
  cex = 0.7,
  col = dotCols,
  pch = dotShape,
  cex.lab = plotCex,
  cex.axis = plotCex,
  cex.main = plotCex
)
mtext(side = 3, text = paste0(hicds, " - ", exprds))
addCorr(x = myx, y=myy, bty="n", legPos="topleft"  )
abline(lm(myy~myx), lty=2, col="darkgrey")


foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))









##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))
