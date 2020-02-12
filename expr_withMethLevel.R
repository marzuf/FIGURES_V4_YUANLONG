
startTime <- Sys.time()
script_name <- "corrExpr_geneWithTAD_condCol.R"

cat(">START ", script_name, "\n")


# Rscript expr_withMethLevel.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR

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
}
col1 <- pal_futurama()(5)[1]
col2 <- pal_futurama()(5)[5]
col1 <- pal_aaas()(5)[4]
col2 <- pal_npg()(5)[5]

aggregFunc <- "max"

dotShape <- 16
shapeMut <- 15
shapeWt <- 8


withMut <- TRUE
withMut_suffix <- ifelse(withMut, "_withMut", "")

require(reshape2)

met_dt <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/luad_meth_sftpa12.RData"))
met_dt_m <- melt(met_dt, id=colnames(met_dt)[1:2])
colnames(met_dt_m)[colnames(met_dt_m) == "variable"] <- "sample"
colnames(met_dt_m)[colnames(met_dt_m) == "value"] <- "beta_met"

agg_met_dt <- aggregate(beta_met ~ FantomPromAnnot_Gene_Name + sample, FUN=aggregFunc, data = met_dt_m[,c("FantomPromAnnot_Gene_Name", "sample", "beta_met")])

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

stopifnot(!duplicated(gff_dt$entrezID))
stopifnot(!duplicated(gff_dt$symbol))

stopifnot(agg_met_dt$FantomPromAnnot_Gene_Name %in% gff_dt$symbol)

entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)

agg_met_dt$entrezID <- symb2entrez[paste0(agg_met_dt$FantomPromAnnot_Gene_Name)]

outFolder <- "EXPR_WITH_METLEVEL"
dir.create(outFolder, recursive = TRUE)

runFolder <- "../v2_Yuanlong_Cancer_HiC_data_TAD_DA"

pipFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path(runFolder, "PIPELINE", "INPUT_FILES")

settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
stopifnot(file.exists(settingFile))
source(settingFile)

samp1_id <- get(load(file.path(setDir, sample1_file)))
samp2_id <- get(load(file.path(setDir, sample2_file)))

g2t_file <- file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt")
g2t_dt <- read.delim(g2t_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
stopifnot(!duplicated(g2t_dt$entrezID))
g2t_vect <- setNames(g2t_dt$region, g2t_dt$entrezID)


gene_list <- get(load(file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata") ))

stopifnot(agg_met_dt$entrezID %in% gene_list)



count_file <- file.path(pipFolder, hicds, exprds, "1_runGeneDE", "DE_rnaseqDT.Rdata")
stopifnot(file.exists(count_file))
fpkm_dt <- get(load(count_file))

stopifnot(agg_met_dt$entrezID %in% rownames(fpkm_dt))

entrez_to_plot <- unique(agg_met_dt$entrezID)

fpkm_plot_dt <- fpkm_dt[paste0(entrez_to_plot),]
stopifnot(nrow(fpkm_plot_dt) == length(entrez_to_plot))
count_m <- melt(t(fpkm_plot_dt))
colnames(count_m)[colnames(count_m) == "Var1" ] <- "sample"
colnames(count_m)[colnames(count_m) == "Var2" ] <- "entrezID"
colnames(count_m)[colnames(count_m) == "value" ] <- "count"

stopifnot(count_m$entrezID %in% entrez_to_plot)

plot_dt <- merge(count_m, agg_met_dt,  by = c("sample", "entrezID"))

plot_dt$cond <- ifelse(plot_dt$sample %in% samp1_id, cond1, 
                          ifelse(plot_dt$sample %in% samp2_id, cond2, NA))
stopifnot(!is.na(plot_dt$cond))

plot_dt$dotCols <- ifelse(plot_dt$sample %in% samp1_id, col1, 
                  ifelse(plot_dt$sample %in% samp2_id, col2, NA))
stopifnot(!is.na(plot_dt$dotCols))

mutSamples <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/NFE2L2_KEAP1_MUTSAMPLES/mut_samples.Rdata"))


plot_dt$dotShape <- ifelse(plot_dt$sample %in% mutSamples, shapeMut, shapeWt)
stopifnot(!is.na(plot_dt$dotShape))

if(!withMut) plot_dt$dotShape <- dotShape


plot_entrez =entrez_to_plot[1]
for(plot_entrez in entrez_to_plot) {
  
  sub_dt <- plot_dt[plot_dt$entrezID == plot_entrez,]
  stopifnot(!duplicated(sub_dt$sample))
  
  plot_region <- as.character(g2t_vect[paste0(plot_entrez)] )
  stopifnot(!is.na(plot_region))
  plot_symbol <- as.character(entrez2symb[paste0(plot_entrez)])
  stopifnot(!is.na(plot_symbol))
  
  myx <- log10(sub_dt$count+log10_offset)
  # myx <- sub_dt$count
  myy <- sub_dt$beta_met
  
  outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", plot_symbol, "_", plot_region, "_exprFPKM_vs_betaMet_agg", aggregFunc, "_colplot", withMut_suffix, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
  
  par(bty="L")
  # par(mar=c(5,4,2,6))
  plot(
    x = myx,
    y = myy,
    main = "Expression and methylation",
    xlab = paste0(plot_symbol, " - RNA-seq FPKM [log10]"),
    ylab = paste0("\u0392 value methylation (", aggregFunc, ")"),
    cex = 0.7,
    pch = sub_dt$dotShape,
    col= sub_dt$dotCols,
    cex.lab = plotCex,
    cex.axis = plotCex,
    cex.main = plotCex
  )
  mtext(side = 3, text = paste0(hicds, " - ", exprds, " - ", cond1, ": ", sum(sub_dt$cond == cond1), "/", length(cond1_ID),"; ", cond2, ": ", sum(sub_dt$cond == cond2), "/", length(cond2_ID) ))
  addCorr(x = myx, y=myy, bty="n", legPos="topleft"  )
  abline(lm(myy~myx), lty=2, col="darkgrey")  
  
  par(xpd=TRUE)
  
  if(withMut) {
    legend(
      "bottomleft",
      # c(max(myx), max(myy)),
      legend= c(cond1, cond2, "KEAP1|NFE2L mut.", "KEAP1|NFE2L wt."),
      pch = c(NA, NA, shapeMut, shapeWt),
      text.col = c(col1, col2, "black", "black"),
      bty="n"
    )
    
  } else {
    legend(
      "bottomleft",
      # c(max(myx), max(myy)),
      legend= c(cond1, cond2),
      pch = c(dotShape),
      text.col = c(col1, col2),
      bty="n"
    )
    
  }
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}


##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))
