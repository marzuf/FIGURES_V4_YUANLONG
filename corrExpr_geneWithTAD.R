
startTime <- Sys.time()
script_name <- "corrExpr_geneWithTAD.R"

cat(">START ", script_name, "\n")


# Rscript corrExpr_geneWithTAD.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr11_TAD390 SOX2

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
  tad_to_plot <- "chr11_TAD390"
  symbol_to_plot <- "SOX2"
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

outFolder <- "CORREXPR_GENEWITHTAD"
dir.create(outFolder, recursive = TRUE)

runFolder <- "../v2_Yuanlong_Cancer_HiC_data_TAD_DA"

pipFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path(runFolder, "PIPELINE", "INPUT_FILES")

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

for(i in 1:nrow(plot_de_dt)) {
  
  curr_entrez <- plot_de_dt$genes[i]
  stopifnot(curr_entrez %in% rownames(fpkm_dt))
  gene_count <- as.numeric(fpkm_dt[paste0(curr_entrez),])
  
  curr_symbol <- plot_de_dt$symbol[i]
  
  stopifnot(length(tf_count) == length(gene_count))

  myx <- log10(tf_count+log10_offset)
  myy <- log10(gene_count+log10_offset)
  
  outFile <- file.path(outFolder, paste0(symbol_to_plot, "_vs_", curr_symbol, "_", tad_to_plot, "_fpkmCount_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
  
  
  par(bty="L")
  
  densplot(
    x = myx,
    y = myy,
    main = "RNA-seq FPKM [log10]",
    xlab = paste0(symbol_to_plot, " (reg. element)"),
    ylab = paste0(curr_symbol, " (", tad_to_plot, " gene)"),
    cex = 0.7,
    cex.lab = plotCex,
    cex.axis = plotCex,
    cex.main = plotCex
  )
  mtext(side = 3, text = paste0(hicds, " - ", exprds))
  addCorr(x = myx, y=myy, bty="n", legPos="topleft"  )
  abline(lm(myy~myx), lty=2, col="darkgrey")  
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}
stopifnot(plot_de_dt$genes %in% rownames(fpkm_dt))
tad_count_dt <- fpkm_dt[rownames(fpkm_dt) %in% plot_de_dt$genes,]

tad_mean_count <- as.numeric(colMeans(tad_count_dt))
stopifnot(length(tf_count) == length(tad_mean_count))


myx <- log10(tf_count+log10_offset)
myy <- log10(tad_mean_count+log10_offset)

outFile <- file.path(outFolder, paste0(symbol_to_plot, "_vs_", tad_to_plot, "_mean_fpkmCount_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))

par(bty="L")

densplot(
  x = myx,
  y = myy,
  main = "RNA-seq FPKM [log10]",
  xlab = paste0(symbol_to_plot, " (reg. element)"),
  ylab = paste0(tad_to_plot, " mean counts"),
  cex = 0.7,
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
