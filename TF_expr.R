
startTime <- Sys.time()
script_name <- "TF_expr.R"

cat(">START ", script_name, "\n")

# Rscript TF_expr.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad
# Rscript TF_expr.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR

require(reshape2)
require(ggplot2)
require(ggsci)
require(foreach)
require(doMC)
registerDoMC(40)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 2) {
  hicds <- args[1] 
  exprds <- args[2] 
} else {
  hicds <-  "ENCSR489OCU_NCI-H460_40kb"
  exprds <- "TCGAluad_norm_luad"
  exprds <- "TCGAlusc_norm_lusc"
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

stopifnot(!duplicated(gff_dt$entrezID))
stopifnot(!duplicated(gff_dt$symbol))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)

outFolder <- "TF_EXPR"
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

akr1c_dt <- reg_dt[reg_dt$targetSymbol %in% c("AKR1C1", "AKR1C2", "AKR1C3"),]
sort(table(akr1c_dt$regSymbol), decreasing = TRUE)
akr1c_tf_dt <- unique(akr1c_dt[,c("regSymbol", "regEntrezID")])
stopifnot(!duplicated(akr1c_tf_dt$regEntrezID))
stopifnot(!duplicated(akr1c_tf_dt$regSymbol))
akr1c_tf <- setNames(akr1c_tf_dt$regSymbol, akr1c_tf_dt$regEntrezID)

hoxb_dt <- reg_dt[reg_dt$targetSymbol %in% c( "HOXB2","HOXB3","HOXB4","HOXB5","HOXB6","HOXB7"),]
sort(table(hoxb_dt$regSymbol), decreasing = TRUE)
hoxb_tf_dt <- unique(hoxb_dt[,c("regSymbol", "regEntrezID")])
stopifnot(!duplicated(hoxb_tf_dt$regEntrezID))
stopifnot(!duplicated(hoxb_tf_dt$regSymbol))
hoxb_tf <- setNames(hoxb_tf_dt$regSymbol, hoxb_tf_dt$regEntrezID)

sftpa_dt <- reg_dt[reg_dt$targetSymbol %in% c( "SFTPA2","SFTPA1"),]
sort(table(sftpa_dt$regSymbol), decreasing = TRUE)
sftpa_tf_dt <- unique(sftpa_dt[,c("regSymbol", "regEntrezID")])
stopifnot(!duplicated(sftpa_tf_dt$regEntrezID))
stopifnot(!duplicated(sftpa_tf_dt$regSymbol))
sftpa_tf <- setNames(sftpa_tf_dt$regSymbol, sftpa_tf_dt$regEntrezID)

mmp_dt <- reg_dt[reg_dt$targetSymbol %in% c( "MMP1", "MMP12","MMP13"),]
sort(table(mmp_dt$regSymbol), decreasing = TRUE)
mmp_tf_dt <- unique(mmp_dt[,c("regSymbol", "regEntrezID")])
stopifnot(!duplicated(mmp_tf_dt$regEntrezID))
stopifnot(!duplicated(mmp_tf_dt$regSymbol))
mmp_tf <- setNames(mmp_tf_dt$regSymbol, mmp_tf_dt$regEntrezID)

if(grepl("norm_", exprds)) {
  all_tf <- c( sftpa_tf, mmp_tf)
  
} else {
  all_tf <- c(akr1c_tf, hoxb_tf)
}


all_tf_dt <- data.frame(tf_symbol = as.character(all_tf), tf_entrez = names(all_tf), stringsAsFactors = FALSE)
all_tf_dt <- unique(all_tf_dt)
stopifnot(!duplicated(all_tf_dt$tf_symbol))
stopifnot(!duplicated(all_tf_dt$tf_entrez))

all_tf_unique <- setNames(all_tf_dt$tf_entrez, all_tf_dt$tf_symbol)



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

# all_tf_unique = setNames(c("9817", "4780"), c("KEAP1", "NFE2L2"))
# all_tf_unique = setNames(c("2551"), c( "GABPA"))

tf = all_tf_unique[1]
out_dt <- foreach(tf = all_tf_unique, .combine='rbind')%dopar%{

  
  stopifnot(tf %in% rownames(fpkm_dt))
    
  
  if(grepl("norm_", exprds)) {
    
    subTit <- paste0(
                     "MMP* TF: ", as.character(tf %in% names(mmp_tf)),
                     "; SFTPA* TF: ", as.character(tf %in% names(sftpa_tf)))
    
    
    
  } else {
    subTit <- paste0("AKR1C* TF: ", as.character(tf %in% names(akr1c_tf)),
                     "; HOXB* TF: ", as.character(tf %in% names(hoxb_tf)))
    
  }
  
  # subTit <- paste0("AKR1C* TF: ", as.character(tf %in% names(akr1c_tf)),
  #                "; HOXB* TF: ", as.character(tf %in% names(hoxb_tf)),
  #                "; MMP* TF: ", as.character(tf %in% names(mmp_tf)),
  #                "; SFTPA* TF: ", as.character(tf %in% names(sftpa_tf)))
  # 
  
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
  
  outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tf_symbol, "_allSamples_exprValues_boxplot.", plotType))
  ggsave(plot = p_var_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
  cat(paste0("... written: ", outFile, "\n"))
  
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

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_allTF_allSamples_outDT.txt"))
# outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_KEAP1_NFE2L2_allSamples_outDT.txt"))
# outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_NRF2_allSamples_outDT.txt"))
write.table(out_dt, file=outFile, sep="\t", quote=F, append=F, col.names=TRUE, row.names=FALSE)
cat(paste0("... written: ", outFile, "\n"))






##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))
