

options(scipen=100)


# Rscript meanFC_meanCorr_dotSize.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr11_TAD390 chr10_TAD268
# Rscript meanFC_meanCorr_dotSize.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD390 chr10_TAD268
# Rscript meanFC_meanCorr_dotSize.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD16 chr17_TAD162
# Rscript meanFC_meanCorr_dotSize.R ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker chr10_TAD16 chr17_TAD162
# 




script_name <- "meanFC_meanCorr_dotSize.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(ggplot2)
require(reshape2)

registerDoMC(4)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

all_cols[all_cols == "red"] <- "brown3"
all_cols[all_cols == "blue"] <- "darkblue"
all_cols[all_cols == "green"] <- "forestgreen"

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

setDir <- "/media/electron"
setDir <- ""

mainFolder <- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

outFolder <- "MEANFC_MEANCORR_DOTSIZE"
dir.create(outFolder, recursive = TRUE)

ex_hicds <- "ENCSR489OCU_NCI-H460_40kb"
ex_exprds <- "TCGAlusc_norm_lusc"
hicds_tit <- "ENCSR489OCU_NCI-H460_40kb"
exprds_tit <- "TCGAlusc_norm_lusc"
tads_to_annot <- c("chr11_TAD390", "chr10_TAD268")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)
ex_hicds <- args[1]
ex_exprds <- args[2]
hicds_tit <- ex_hicds
exprds_tit <- ex_exprds
if(length(args) > 2) {
  tads_to_annot <- args[3:length(args)]  
}else{
  tads_to_annot <- NULL
}



plotCex <- 1.2

signifThresh <- 0.01


final_DT <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/CREATE_FINAL_TABLE/all_result_dt.Rdata"))
final_DT_signif <- final_DT[final_DT$adjPvalComb <= signifThresh,]

cat(paste0(ex_hicds, "\n", ex_exprds, "\n"))

ex_DT <- final_DT[final_DT$hicds == ex_hicds &
                    final_DT$exprds == ex_exprds,
                  ]

stopifnot(nrow(ex_DT) > 0)

ex_DT <- ex_DT[order(ex_DT$adjPvalComb),]



ex_DT_signif <- ex_DT[ex_DT$adjPvalComb <= signifThresh,]

stopifnot(nrow(ex_DT_signif) > 0)



# if <= 0.01 => log10 <= 2 => cex = 0.7
ex_DT$adjPvalComb_log10 <- -log10(ex_DT$adjPvalComb)

outFile <- file.path(outFolder, paste0(ex_hicds, "_", ex_exprds, "_dt.Rdata"))
save(ex_DT, file=outFile, version=2)



notSignifSize <- 0.6
new_min <-0.7
new_max <- 3

# ex_DT$dotSize <- (new_max-new_min) * (ex_DT$adjPvalComb_log10- min(ex_DT$adjPvalComb_log10))/(max(ex_DT$adjPvalComb_log10)-min(ex_DT$adjPvalComb_log10)) + new_min

# signif_min <- min(ex_DT$adjPvalComb_log10[ex_DT$adjPvalComb_log10 >= 2])
signif_min <- 2
signif_max <- max(ex_DT$adjPvalComb_log10[ex_DT$adjPvalComb_log10 >= 2])

ex_DT$dotSize <- (new_max-new_min) * (ex_DT$adjPvalComb_log10 - signif_min)/(signif_max-signif_min) + new_min

ex_DT$dotSize[ex_DT$dotSize <= new_min ] <- notSignifSize

# myvect <- c(1, 2.2,2.3,3,6)
# (new_max-new_min) * (myvect- min(myvect))/(max(myvect)-min(myvect)) + new_min

ex_DT$dotcol <- ifelse(ex_DT$adjPvalComb <= signifThresh,"red", "darkgrey")

save(ex_DT, file="ex_DT.Rdata", version=2)

xoffset <- 1
yoffset <- 0.2

outFile <- file.path(outFolder, paste0(ex_hicds, "_", ex_exprds, "_meanIntraCorr_meanLogFC_dotplot_with_signif.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

par(bty="l")
plot(
  main = paste0(hicds_tit, " - ", exprds_tit),
  x = ex_DT$meanLogFC,
  y = ex_DT$meanCorr,
  xlim = c(min(ex_DT$meanLogFC) - xoffset , max(ex_DT$meanLogFC) + xoffset),
  ylim = c(min(ex_DT$meanCorr) -yoffset , max(ex_DT$meanCorr) + yoffset),
  xlab = "meanLogFC",
  ylab = "meanIntraCorr",
  pch=16,
  col = ex_DT$dotcol,
  cex=ex_DT$dotSize,
  cex.axis=plotCex,
  cex.lab=plotCex
)
points(x = ex_DT$meanLogFC[1:10],
     y = ex_DT$meanCorr[1:10],
     col="blue", 
     cex=ex_DT$dotSize[1:10], 
     lwd=2)
mtext(side=3, text = paste0("# TADs = ", nrow(ex_DT), "; # signif. TADs = ", nrow(ex_DT_signif)))

legend("bottomright", legend=c(paste0("adj. p-val <= ", signifThresh), "top 10"), pch=c(16, 1), col=c("red", "blue"), bty="n")

pval_001 <- (new_max-new_min) * (-log10(0.01)- signif_min)/(signif_max-signif_min) + new_min
pval_0001 <- (new_max-new_min) * (-log10(0.001)- signif_min)/(signif_max-signif_min) + new_min
pval_00001 <- (new_max-new_min) * (-log10(0.0001)- signif_min)/(signif_max-signif_min) + new_min

legend("bottomleft", 
       legend = c( as.expression(bquote("p-val = " ~ 10^-2)),
                  as.expression(bquote("p-val = " ~ 10^-3))
                                           ),
       pch = 16, 
       pt.cex = c(pval_001, pval_0001),
       col = "red",
       bty="n"
       )

if(!is.null(tads_to_annot)) {
  
  tads_to_annot <- tads_to_annot[tads_to_annot %in% ex_DT$region]
  
  text(
    x = ex_DT$meanLogFC[ex_DT$region %in% tads_to_annot],
    y = ex_DT$meanCorr[ex_DT$region %in% tads_to_annot],
    labels = ex_DT$region[ex_DT$region %in% tads_to_annot],
    pos=3
  )
  
}

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

nAxBreaks <- 10

x_breaks <- scales::pretty_breaks(n = nAxBreaks)(range(ex_DT$meanLogFC) + c(-xoffset, xoffset))
y_breaks <- scales::pretty_breaks(n = nAxBreaks)(range(ex_DT$meanCorr) + c(-yoffset, yoffset))

outFile <- file.path(outFolder, paste0(ex_hicds, "_", ex_exprds, "_meanIntraCorr_meanLogFC_dotplot_with_signif_axBreaks.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

par(bty="l")
plot(
  main = paste0(hicds_tit, " - ", exprds_tit),
  x = ex_DT$meanLogFC,
  y = ex_DT$meanCorr,
  xlim = c(min(ex_DT$meanLogFC) - xoffset , max(ex_DT$meanLogFC) + xoffset),
  ylim = c(min(ex_DT$meanCorr) -yoffset , max(ex_DT$meanCorr) + yoffset),
  xlab = "meanLogFC",
  ylab = "meanIntraCorr",
  pch=16,
  col = ex_DT$dotcol,
  cex=ex_DT$dotSize,
  cex.axis=plotCex,
  cex.lab=plotCex,
  axes=FALSE
)
axis(1, at = x_breaks, lwd=0, lwd.ticks = 1, cex = plotCex, cex.axis=plotCex, cex.lab=plotCex)
axis(2, at = y_breaks, lwd=0, lwd.ticks = 1, cex = plotCex, cex.axis=plotCex, cex.lab=plotCex)
box(bty="L")
points(x = ex_DT$meanLogFC[1:10],
       y = ex_DT$meanCorr[1:10],
       col="blue", 
       cex=ex_DT$dotSize[1:10], 
       lwd=2)
mtext(side=3, text = paste0("# TADs = ", nrow(ex_DT), "; # signif. TADs = ", nrow(ex_DT_signif)))

legend("bottomright", legend=c(paste0("adj. p-val <= ", signifThresh), "top 10"), pch=c(16, 1), col=c("red", "blue"), bty="n")

pval_001 <- (new_max-new_min) * (-log10(0.01)- signif_min)/(signif_max-signif_min) + new_min
pval_0001 <- (new_max-new_min) * (-log10(0.001)- signif_min)/(signif_max-signif_min) + new_min
pval_00001 <- (new_max-new_min) * (-log10(0.0001)- signif_min)/(signif_max-signif_min) + new_min

legend("bottomleft", 
       legend = c( as.expression(bquote("p-val = " ~ 10^-2)),
                   as.expression(bquote("p-val = " ~ 10^-3))
       ),
       pch = 16, 
       pt.cex = c(pval_001, pval_0001),
       col = "red",
       bty="n"
)

if(!is.null(tads_to_annot)) {
  
  tads_to_annot <- tads_to_annot[tads_to_annot %in% ex_DT$region]
  
  text(
    x = ex_DT$meanLogFC[ex_DT$region %in% tads_to_annot],
    y = ex_DT$meanCorr[ex_DT$region %in% tads_to_annot],
    labels = ex_DT$region[ex_DT$region %in% tads_to_annot],
    pos=3
  )
  
}

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




