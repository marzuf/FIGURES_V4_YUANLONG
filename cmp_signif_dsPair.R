


# EGFR mutations are most common in people with lung adenocarcinoma (a form of non-small cell lung cancer), are more common with lung cancer in non-smokers, and are more common in women than in men.
#KRAS mutations are found in ~ 25% of lung adenocarcinomas in Western countries and, as a group, have been strongly associated with cigarette smoking. These mutations are predictive of poor prognosis in resected disease as well as resistance to treatment with erlotinib or gefitinib.

# KRAS ~ smoker
# EGFR ~ non smoker

# Rscript cmp_signif_dsPair.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad 0.01
# Rscript cmp_signif_dsPair.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker 0.01
# Rscript cmp_signif_dsPair.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc 0.01

require(VennDiagram)


plotType <- "svg"
source("../FIGURES_V2_YUANLONG/settings.R")

myWidth <- myWidth * 1.2

outFolder <- "CMP_SIGNIF_DSPAIR"
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

ggsci_pal <- "d3"
ggsci_subpal <- ""

plotMargin <- c(1,2,1,1)

hicds1 <- "LG1_40kb"
exprds1 <- "TCGAlusc_norm_lusc"
hicds2 <- "LG1_40kb"
exprds2 <- "TCGAluad_norm_luad"
tadPval_thresh <- 0.01

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  
  
} else {
  stopifnot(length(args) == 5)
  hicds1 <- args[1]
  exprds1 <- args[2]
  hicds2 <- args[3]
  exprds2 <- args[4]
  tadPval_thresh <- as.numeric(args[5])
  stopifnot(!is.na(tadPval_thresh))
  
}

# to compare the TADs -> they must have same hicds
stopifnot(hicds1 == hicds2)

inDT <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
inDT_check <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/CREATE_FINAL_TABLE/all_result_dt.Rdata"))


signif_ds1 <- inDT[inDT$hicds == hicds1 & inDT$exprds == exprds1 & inDT$tad_adjCombPval <= tadPval_thresh,]
signif_ds1_check <- inDT_check[inDT_check$hicds == hicds1 & inDT_check$exprds == exprds1 & inDT_check$adjPvalComb <= tadPval_thresh,]
stopifnot(setequal(signif_ds1$region, signif_ds1_check$region))

signif_ds2 <- inDT[inDT$hicds == hicds2 & inDT$exprds == exprds2 & inDT$tad_adjCombPval <= tadPval_thresh,]
signif_ds2_check <- inDT_check[inDT_check$hicds == hicds2 & inDT_check$exprds == exprds2 & inDT_check$adjPvalComb <= tadPval_thresh,]
stopifnot(setequal(signif_ds2$region, signif_ds2_check$region))

ds1_dt <- inDT_check[inDT_check$hicds == hicds1 & inDT_check$exprds == exprds1 ,]
ds2_dt <-  inDT_check[inDT_check$hicds == hicds2 & inDT_check$exprds == exprds2 ,]

######################################################################
#### CMP_PVAL
######################################################################

merge_dt <- merge(ds1_dt[,c("hicds",  "region", "adjPvalComb")], 
                  ds2_dt[,c("hicds", "region", "adjPvalComb")], 
                  by=c("hicds",  "region"), 
                  # suffixes = c(paste0("_", exprds1), paste0("_", exprds2)),
                  suffixes = c(paste0("_", "ds1"), paste0("_", "ds2")),
                  # all.x=T, all.y=T
                  all.x=F, all.y=F
                  )


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


ds1 <- paste0(hicds1, "-", exprds1)
ds2 <- paste0(hicds2, "-", exprds2)

myxlab <- paste0(exprds1)
myylab <- paste0(exprds2)

mysub <- paste0(hicds1)

myx <- -log10(merge_dt$adjPvalComb_ds1)
myy <- -log10(merge_dt$adjPvalComb_ds2)

outFile <- file.path(outFolder, paste0(hicds1, "_", exprds1, "_vs_", hicds2, "_", exprds2, "_adjPvalComb_densplot.", "png"))
do.call("png", list(outFile, height=400, width=400))

par(bty="L", family=fontFamily)

densplot(
  x=myx,
  y=myy,
  xlab = myxlab,
  ylab = myylab,
  main = paste0("TAD adj. comb. p-val. [log10]"),
  pch =16,
  cex.axis = plotCex,
  cex.lab = plotCex,
  cex.main = plotCex
         )
addCorr(x=myx,
        y=myy,
        legPos = "bottomright",
        bty="n"
        )
mtext(side=3, text = paste0(mysub), cex=plotCex)
legend("topleft",
       legend=paste0("# common TADs = ", nrow(merge_dt)),
       bty="n")
abline(
  lm(myy~myx),
  lty=2,
  col = "darkgrey"
)

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

######################################################################
#### VENN DIAGRAM
######################################################################
stopifnot(!duplicated(signif_ds1$entrezID))
stopifnot(!duplicated(signif_ds2$entrezID))

signif_ds1 <- unique(signif_ds1$entrezID)
signif_ds2 <- unique(signif_ds2$entrezID)

signif_ds1Only <- setdiff(signif_ds1, signif_ds2)
stopifnot(signif_ds1Only %in% signif_ds1 & ! signif_ds1Only %in% signif_ds2)
nsignif_ds1Only <- length(signif_ds1Only)

signif_ds2Only <- setdiff(signif_ds2, signif_ds1)
stopifnot(signif_ds2Only %in% signif_ds2 & ! signif_ds2Only %in% signif_ds1)
nsignif_ds2Only <- length(signif_ds2Only)

signif_ds1And2 <- intersect(signif_ds1, signif_ds2)
nsignif_ds1And2 <- length(signif_ds1And2)

stopifnot(nsignif_ds1And2+nsignif_ds2Only+nsignif_ds1Only == length(unique(c(signif_ds1, signif_ds2))))

outFile <- file.path(outFolder, paste0(hicds1, "_", exprds1, "_vs_", hicds2, "_", exprds2, "_nbrSignifGenes_vennDiagram.", plotType))

require(ggplot2)
require(ggsci)

ds1 <- paste0(hicds1, "-", exprds1)
ds2 <- paste0(hicds2, "-", exprds2)

myTit <- paste0("# signif. genes\nTAD p-val <= ", tadPval_thresh )
myTit <- paste0("# signif. genes")
mySub <- paste0(ds1, "\n", ds2)
mySub <- paste0("TAD p-val <= ", tadPval_thresh )


vd <- venn.diagram(
  x = list(signif_ds1, signif_ds2),
  main = myTit,
  sub = mySub,
  category.names = c(paste0("signif. ", gsub("-", "\n", ds1), "\n(", nsignif_ds1Only+nsignif_ds1And2, ")") , paste0("signif. ", gsub("-", "\n", ds2), "\n(", nsignif_ds2Only+nsignif_ds1And2, ")")),
  
  fontfamily="Hershey",
  cat.fontfamily="Hershey",
  main.fontfamily="Hershey",
  sub.fontfamily="Hershey",
  
  sub.fontface="plain",
  main.fontface="bold",
  
  margin=c(0,0,0,0),
  
  main.cex = 2,
  sub.cex=1.6,
  cat.cex=1.4,
  cex = 2,
  
  cat.default.pos="outer",
  cat.pos=c(-25, 25),
  
  cat.col = pal_lancet()(2),
  fill = pal_lancet()(2),
  col = pal_lancet()(2),
  alpha=0.2,
  scaled=FALSE,
  filename = NULL
)
ggsave(vd, file=outFile,height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

system(paste0("rm -f VennDiagram2020*.log"))



