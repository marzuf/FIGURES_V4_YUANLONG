
startTime <- Sys.time()
script_name <- "TF_expr_FC.R"

cat(">START ", script_name, "\n")

# Rscript TF_expr_symbol_withMutStat.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad NKX2-1
# Rscript TF_expr_symbol_withMutStat.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR

require(reshape2)
require(ggplot2)
require(ggsci)
require(foreach)
require(doMC)
registerDoMC(40)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 3) {
  hicds <- args[1] 
  exprds <- args[2] 
  symbol <- args[3]
} else {
  hicds <-  "ENCSR489OCU_NCI-H460_40kb"
  exprds <- "TCGAluad_norm_luad"
  exprds <- "TCGAlusc_norm_lusc"
  symbol <- "NKX2-1"
}
col1 <- pal_futurama()(5)[1]
col2 <- pal_futurama()(5)[5]
col1 <- pal_aaas()(5)[4]
col2 <- pal_npg()(5)[5]


# my_xlab <- "TAD genes (ordered by start positions)"
my_xlab <- ""
my_ylab <- "RNA-seq expression count [log10]"

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
stopifnot(symbol %in% gff_dt$symbol)

mutSamples <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/NFE2L2_KEAP1_MUTSAMPLES/mut_samples.Rdata"))
mutCol <- "chartreuse3"
notmutCol <- "darkgrey"


stopifnot(!duplicated(gff_dt$entrezID))
stopifnot(!duplicated(gff_dt$symbol))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)

outFolder <- "TF_EXPR_SYMBOL_WITHMUTSTAT"
dir.create(outFolder, recursive = TRUE)

runFolder <- "../v2_Yuanlong_Cancer_HiC_data_TAD_DA"

pipFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path(runFolder, "PIPELINE", "INPUT_FILES")

reg_file <- file.path(runFolder, "chea3_lung_TFs_processed.txt")
reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)

cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
cat(paste0("with Target Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)

reg_dt <- reg_dt[reg_dt$regSymbol %in% names(symb2entrez),]
cat(paste0("with Reg Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
reg_dt$regEntrezID <- symb2entrez[reg_dt$regSymbol]
reg_dt$regEntrezID <- as.character(reg_dt$regEntrezID)

all_reg_entrez <- unique(reg_dt$regEntrezID)

g2t_file <- file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt")
g2t_dt <- read.delim(g2t_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2t_dt$entrezID <- as.character(g2t_dt$entrezID)


## => change here, take only the expressed TFs
de_dt <- get(load(file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "1_runGeneDE", "DE_topTable.Rdata") ))

reg_dt <- reg_dt[reg_dt$regEntrezID %in% de_dt$genes,]
cat(paste0("with expressed Reg Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))

all_tf_dt <- data.frame(tf_symbol = as.character(symbol), tf_entrez = gff_dt$entrezID[gff_dt$symbol == symbol], stringsAsFactors = FALSE)
stopifnot(nrow(all_tf_dt) == 1)
all_tf_dt <- unique(all_tf_dt)
stopifnot(!duplicated(all_tf_dt$tf_symbol))
stopifnot(!duplicated(all_tf_dt$tf_entrez))

all_tf_unique <- setNames(all_tf_dt$tf_entrez, all_tf_dt$tf_symbol)


stopifnot(all_tf_dt$tf_entrez %in% de_dt$genes)
tmp_dt <- de_dt
tmp_dt$tf_entrez <- tmp_dt$genes
tmp_dt <- tmp_dt[order(abs(tmp_dt$logFC), decreasing = TRUE),]
tmp_dt$absLogFC_rank <- rank(-abs(tmp_dt$logFC))
rank_dt <- merge(all_tf_dt, tmp_dt[,c("tf_entrez", "logFC", "absLogFC_rank")], by=c("tf_entrez"), all.x=TRUE, all.y=FALSE)
stopifnot(!is.na(rank_dt))

all_tf_rank <- setNames(rank_dt$absLogFC_rank, rank_dt$tf_entrez)
all_tf_rank <- sort(all_tf_rank)

count_file <- file.path(pipFolder, hicds, exprds, "1_runGeneDE", "DE_rnaseqDT.Rdata")
stopifnot(file.exists(count_file))
fpkm_dt <- get(load(count_file))

stopifnot(all_tf_unique %in% rownames(fpkm_dt))

settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
stopifnot(file.exists(settingFile))
source(settingFile)

samp1 <- get(load(file.path(setDir, sample1_file)))
samp2 <- get(load(file.path(setDir, sample2_file)))

## => boxplot for selected TF elements here !


stopifnot(setequal(all_tf_unique, names(all_tf_rank)))

tf = all_tf_unique[1]
out_dt <- foreach(tf = names(all_tf_rank), .combine='rbind')%do%{

  
  tf_rank <- all_tf_rank[paste0(tf)]
  
  stopifnot(tf %in% rownames(fpkm_dt))
  
  stopifnot(tf %in% names(all_tf_rank))  
  
  subTit <- paste0("TF rank: ", tf_rank)
  
  fpkm_plot_dt <- fpkm_dt[rownames(fpkm_dt) == tf,,drop=FALSE]
  fpkm_plot_dt$entrezID <- rownames(fpkm_plot_dt)
  
  stopifnot(samp1 %in% colnames(fpkm_plot_dt))
  stopifnot(samp2 %in% colnames(fpkm_plot_dt))
  
  fpkm_plot_dt <- fpkm_plot_dt[,c("entrezID", samp1, samp2)]
  m_fpkm_dt <- melt(fpkm_plot_dt, id="entrezID")
  
  stopifnot(m_fpkm_dt$entrezID == tf)
  
  toplot_dt <-  m_fpkm_dt
  stopifnot(!is.na(toplot_dt))
  toplot_dt$cond <- ifelse(toplot_dt$variable %in% samp1, cond1, 
                           ifelse(toplot_dt$variable %in% samp2, cond2, NA ))
  stopifnot(!is.na(toplot_dt$cond))
  toplot_dt$symbol <- names(all_tf_unique)[all_tf_unique==tf]
  stopifnot(length(unique(toplot_dt$symbol)) == 1)
  
  toplot_dt <- toplot_dt[order(toplot_dt$value),]
  toplot_dt$value_log10 <- log10(toplot_dt$value + log10_offset)
  

    
  withRank_toplot_dt2 <- do.call(rbind, by(toplot_dt, list(toplot_dt$symbol), function(x) {
    x$cond <- factor(x$cond, levels=c(cond1,cond2))
    dt <- x[order(-as.numeric(x$cond), x$value, decreasing = TRUE),,drop=FALSE]
    dt$samp_rank <- 1:nrow(dt)
    dt
  }))
  withRank_toplot_dt2$hicds <- hicds
  withRank_toplot_dt2$exprds <- exprds
  
  cat("merge withRank and inDT \n")
  
  withRank_toplot_dt2$cond <- factor(withRank_toplot_dt2$cond, levels = c(cond1,cond2))
  
  tf_symbol <- unique(withRank_toplot_dt2$symbol)
  stopifnot(length(tf_symbol) == 1)
  

  withRank_toplot_dt2$symbol_lab <- paste0(withRank_toplot_dt2$symbol)
  
  save(withRank_toplot_dt2, file ="withRank_toplot_dt2.Rdata")
  
  
  p_var_boxplot <-  ggplot(withRank_toplot_dt2, aes(x = symbol_lab, y = value_log10, fill = cond)) + 
    # geom_boxplot(notch = TRUE, outlier.shape=NA)+
    # geom_jitter(aes(colour = cond), position=position_jitterdodge())+
    
    geom_point(aes(color=cond), position=position_jitterdodge(),  alpha=0.5) +
    
    
    geom_boxplot(notch = TRUE, outlier.shape=NA)+
    ggtitle(paste0(hicds, " - ", exprds), subtitle = paste0(subTit))+
    scale_x_discrete(name=my_xlab)+
    scale_y_continuous(name=paste0(my_ylab),
                       breaks = scales::pretty_breaks(n = 20))+
    
    scale_color_manual(values=c(col1, col2))+
    scale_fill_manual(values=c(col1, col2))+
    
    labs(fill  = paste0("Cond."), color=paste0("Cond.")) +
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
      axis.title.y = element_text(color="black", size=13),
      axis.title.x = element_text(color="black", size=13),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.text = element_text(size=12),
      legend.key = element_blank(),
      legend.title = element_text(face="bold")
    )
  
  outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tf_symbol, "_rank", tf_rank , "_allSamples_exprValues_boxplot.", plotType))
  ggsave(plot = p_var_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
  cat(paste0("... written: ", outFile, "\n"))
  
  withRank_toplot_dt2$cond_dots <- withRank_toplot_dt2$cond
  withRank_toplot_dt2$cond_dots <- ifelse(as.character(withRank_toplot_dt2$variable) %in% mutSamples, "withMut", as.character(withRank_toplot_dt2$cond_dots))
  withRank_toplot_dt2$cond_dots <- factor(withRank_toplot_dt2$cond_dots, levels = c(cond1,cond2, "withMut"))
  stopifnot(!is.na(withRank_toplot_dt2$cond_dots))
  
  
  
  withRank_toplot_dt2$cond_sh <- withRank_toplot_dt2$cond
  withRank_toplot_dt2$cond_sh <- ifelse(as.character(withRank_toplot_dt2$variable) %in% mutSamples, "withMut", "noMut")
  withRank_toplot_dt2$cond_sh <- factor(withRank_toplot_dt2$cond_sh, levels = c("noMut", "withMut"))
  stopifnot(!is.na(withRank_toplot_dt2$cond_sh))
  
  withRank_toplot_dt2$cond_border <- paste0(withRank_toplot_dt2$cond , "_", withRank_toplot_dt2$cond_sh)
  withRank_toplot_dt2$cond_border <- factor(withRank_toplot_dt2$cond_border, levels = c(
    paste0(cond1, "_noMut"),  paste0(cond1, "_withMut"), paste0(cond2, "_noMut"), paste0(cond2, "_withMut")
  ))
  stopifnot(!is.na(withRank_toplot_dt2$cond_border))
  
  
  
  
  
  
  
  
  
  p_var_boxplot <- ggplot(withRank_toplot_dt2, aes(x = symbol_lab, y = value_log10, fill = cond)) +
    geom_point(aes(color=cond_border), position=position_jitterdodge(), stroke=0.8, shape=21, alpha=0.8) +
    geom_boxplot(notch = TRUE, outlier.shape=NA)+
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
  
  outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tf_symbol, "_rank", tf_rank , "_allSamples_exprValues_boxplot_vCol.", plotType))
  ggsave(plot = p_var_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
  
  
  
  
  
  
  
  data.frame(
    hicds = hicds, 
    exprds = exprds,
    TF_entrez = tf,
    TF_symbol = tf_symbol,
    cond1 = cond1,
    cond2= cond2,
    meanExpr_cond1 = mean(withRank_toplot_dt2$value[withRank_toplot_dt2$cond == cond1]),
    meanExpr_cond2 = mean(withRank_toplot_dt2$value[withRank_toplot_dt2$cond == cond2]),
    regGenes = subTit,
    stringsAsFactors = FALSE
  )  
  
}
out_dt$diffMeanExpr <- out_dt$meanExpr_cond1-out_dt$meanExpr_cond2
out_dt$diffMeanExpr_abs <- abs(out_dt$diffMeanExpr)
out_dt <- out_dt[order(out_dt$diffMeanExpr_abs, decreasing=TRUE),]

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", symbol, "_allTF_allSamples_outDT.txt"))
# outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_KEAP1_NFE2L2_allSamples_outDT.txt"))
# outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_NRF2_allSamples_outDT.txt"))
write.table(out_dt, file=outFile, sep="\t", quote=F, append=F, col.names=TRUE, row.names=FALSE)
cat(paste0("... written: ", outFile, "\n"))






##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))
