hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"
tad_to_show <- "chr11_TAD390"

plotType <- "svg"

source("../FIGURES_V2_YUANLONG/settings.R")
require(ggsci)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr11_TAD390
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr10_TAD268
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD390
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr10_TAD268
# 
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD16
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr17_TAD162
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker chr10_TAD16
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker chr17_TAD162
# 


outFolder <- file.path("GENERANK_TADRANK_SELECTTAD")
dir.create(outFolder, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
hicds <- args[1]
exprds <- args[2]
tad_to_show <- args[3]

inDT <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))

ds_dt <- inDT[inDT$hicds == hicds & inDT$exprds == exprds,]

tad_dt <- ds_dt[,c("hicds", "exprds", "region", "tad_adjCombPval")]
tad_dt <- unique(tad_dt)
tad_dt$tad_rank <- rank(tad_dt$tad_adjCombPval, ties="min")
tad_dt$rev_tad_rank <- rank(-tad_dt$tad_adjCombPval, ties="min")

ds_dt <- merge(ds_dt, tad_dt[, c("hicds", "exprds", "region", "rev_tad_rank")], all=TRUE, by=c("hicds", "exprds", "region"))

nTADs <- length(unique(ds_dt$region))

head(ds_dt[order(ds_dt$tad_adjCombPval),])


geneSignifPval <- 0.05
tadSignifPval <- 0.01


tad_to_show_col <-  "red"
hide_col <- "grey"
showCex <- 1.2

nTop <- 10

top_reg <- unique(ds_dt$region[ds_dt$tad_rank <= nTop])
last_reg <- unique(ds_dt$region[ds_dt$rev_tad_rank <= nTop])



ds_dt$showCol <- ifelse(ds_dt$region == tad_to_show, tad_to_show_col, hide_col)

stopifnot(!is.na(ds_dt$dotCol))
plot_dt <- ds_dt[ds_dt$region %in% top_reg | ds_dt$region %in% last_reg ,]


#########################################################################################################################
#### BARPLOTS - all gene ranks
#########################################################################################################################
maxGeneRank <- max(ds_dt$gene_rank)
maxTADrank <- max(ds_dt$tad_rank)

ds_dt$gene_rank_rel <- ds_dt$gene_rank/maxGeneRank
ds_dt$tad_rank_rel <- ds_dt$tad_rank/maxTADrank


geneBar_pos <- 1
tadBar_pos <- 2
axisOffset <- 0.5
labOffset <- 0.02

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_show, "_geneRank_tadRank_showBars_segments.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))

par(bty="n", family=fontFamily)
plot(NULL,
     xlab="",
     ylab="",
     axes=F,
     # main = paste0(hicds, " - " , exprds),
     main = paste0(hicds, "\n" , exprds),
     # sub=paste0(),
     cex.axis=plotCex,
     cex.lab = plotCex,
     cex.main = plotCex,
  xlim = c(geneBar_pos-axisOffset, tadBar_pos+axisOffset),
  # ylim = c(-axisOffset, 1+axisOffset)
  ylim = c(1+axisOffset, -axisOffset )
)
# mtext(side=3, paste0(hicds, " - ", exprds), cex=plotCex)
mtext(side=3, paste0(tad_to_show), cex=plotCex-0.2, line=-1)

segments(x0=geneBar_pos,
         x1=tadBar_pos,
         y0=ds_dt$gene_rank_rel[ds_dt$showCol == tad_to_show_col],
         y1=ds_dt$tad_rank_rel[ds_dt$showCol == tad_to_show_col],
         col = tad_to_show_col
         )
# add bar for the genes
segments(x0=geneBar_pos, y0=0, x1=geneBar_pos, y1=1, lwd=5)
# add bar for the TADs
segments(x0=tadBar_pos, y0=0, x1=tadBar_pos, y1=1, lwd=5)

text(
  x = c(geneBar_pos, tadBar_pos),
  # y = 1 + 0.2,
  y = 0- 0.2,
  labels=c("Gene ranks", "TAD ranks"),
  cex = plotCex
)
show_dt <- ds_dt[ds_dt$region == tad_to_show,]

text(
  y = c(unique(show_dt$tad_rank_rel), unique(show_dt$gene_rank_rel)),
  x = c(rep(tadBar_pos, length(unique(show_dt$tad_rank_rel)))+labOffset, rep(geneBar_pos, length(unique(show_dt$gene_rank_rel)))-labOffset),
  labels=c(unique(show_dt$tad_rank), unique(show_dt$gene_rank)),
  cex = 1,
  pos=c(rep(4, length(unique(show_dt$tad_rank_rel))), rep(2, length(unique(show_dt$gene_rank_rel))))
)
segments(
  y0 = c(unique(show_dt$tad_rank_rel), unique(show_dt$gene_rank_rel)),
  y1 = c(unique(show_dt$tad_rank_rel), unique(show_dt$gene_rank_rel)),
  x0 = c(rep(tadBar_pos, length(unique(show_dt$tad_rank_rel))), rep(geneBar_pos, length(unique(show_dt$gene_rank_rel)))),
  x1 = c(rep(tadBar_pos, length(unique(show_dt$tad_rank_rel)))+labOffset, rep(geneBar_pos, length(unique(show_dt$gene_rank_rel)))-labOffset),
  lwd=1
)

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


#########################################################################################################################
#### BARPLOTS - all gene med ranks
#########################################################################################################################

med_ds_dt <- aggregate(gene_rank ~ hicds + exprds+ region, data=ds_dt, FUN=median)
colnames(med_ds_dt)[colnames(med_ds_dt) == "gene_rank"] <- "med_gene_rank"

med_plot_dt <- merge(med_ds_dt, tad_dt, all=TRUE, by=c("hicds", "exprds","region"))
stopifnot(!duplicated(med_plot_dt$region))


maxGeneRank <- max(med_plot_dt$med_gene_rank)
maxTADrank <- max(med_plot_dt$tad_rank)

med_plot_dt$gene_rank_rel <- med_plot_dt$med_gene_rank/maxGeneRank
med_plot_dt$tad_rank_rel <- med_plot_dt$tad_rank/maxTADrank


med_plot_dt$showCol <- ifelse(med_plot_dt$region == tad_to_show, tad_to_show_col, hide_col)
stopifnot(!is.na(med_plot_dt$showCol))

geneBar_pos <- 1
tadBar_pos <- 2
axisOffset <- 0.5

outFile <- file.path(outFolder, paste0(hicds,"_", exprds, "_", tad_to_show, "_medGeneRank_tadRank_showBars_segments.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))

par(bty="n", family=fontFamily)
plot(NULL,
     xlab="",
     ylab="",
     axes=F,
     # main=paste0("Top and last ", nTop, " TADs"),
     # main=paste0(hicds, " - " , exprds),
     main=paste0(hicds, "\n" , exprds),
     # sub=paste0(),
     cex.axis=plotCex,
     cex.lab = plotCex,
     cex.main = plotCex,
     xlim = c(geneBar_pos-axisOffset, tadBar_pos+axisOffset),
     # ylim = c(-axisOffset, 1+axisOffset)
     ylim = c(1+axisOffset,-axisOffset)
)
# mtext(side=3, paste0(hicds, " - ", exprds), cex=plotCex)
mtext(side=3, paste0(tad_to_show), cex=plotCex-0.2, line=-1)

segments(x0=geneBar_pos,
         x1=tadBar_pos,
         y0=med_plot_dt$gene_rank_rel[med_plot_dt$showCol == hide_col],
         y1=med_plot_dt$tad_rank_rel[med_plot_dt$showCol == hide_col],
         lty=3,
         tcl=-.2,
         col = hide_col
)
segments(x0=geneBar_pos,
         x1=tadBar_pos,
         y0=med_plot_dt$gene_rank_rel[med_plot_dt$showCol ==  tad_to_show_col],
         y1=med_plot_dt$tad_rank_rel[med_plot_dt$showCol == tad_to_show_col],
         col = tad_to_show_col
)



# add bar for the genes
segments(x0=geneBar_pos, y0=0, x1=geneBar_pos, y1=1, lwd=5)

# add bar for the TADs
segments(x0=tadBar_pos, y0=0, x1=tadBar_pos, y1=1, lwd=5)

text(
  x = c(geneBar_pos, tadBar_pos),
  # y = 1 + 0.2,
  y = 0- 0.2,
  labels=c("Median gene rank", "TAD ranks"),
  cex = plotCex
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




