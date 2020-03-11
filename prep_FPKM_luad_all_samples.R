# Rscript prep_FPKM_luad_all_samples.R

# scaled estimates folder
setDir <- ""
cancerType <- "LUAD"
seFolderMain <- file.path(setDir, "/mnt/ndata/databank/TCGA/TCGA_gdac/TCGA_gdac_full_snapshots/", cancerType)
seSubFolders <- list.files(seFolderMain, pattern="^2")
data_release <- "20160128"
seFolder <- file.path(seFolderMain, data_release)
seFile <- list.files(seFolder, full.names = TRUE, pattern=paste0(".+illuminahiseq.+rnaseqv2.+RSEM_genes__data.Level_3.+gz$"))
stopifnot(length(seFile) == 1)


cmd <- paste0("tar -xzvf ", seFile)
cat(paste0(cmd, "\n"))
system(cmd)
stopifnot(file.exists(seFile))

extractFolder <- gsub(".tar.gz$", "", basename(seFile))
seInFile <- file.path(extractFolder, gsub(paste0("^.+", cancerType, ".Merge_"), paste0(cancerType, "."), basename(seFile)))
seInFile <- file.path(extractFolder, gsub("(.+__data\\.)Level_3.+", "\\1data.txt", basename(seInFile)))

seDT <- read.delim(seInFile, stringsAsFactors = FALSE)
colnames(seDT) <- gsub("\\.", "-", colnames(seDT))
stopifnot(colnames(seDT)[1] == "Hybridization-REF")
stopifnot(seDT[1,1] == "gene_id")

colsToKeep <- seDT[1,] %in% c("gene_id", "scaled_estimate")
stopifnot(length(colsToKeep) == ncol(seDT))
seDT <- seDT[, colsToKeep]
stopifnot(ncol(seDT) == sum(colsToKeep))
stopifnot(seDT[which(seDT[,1] == "gene_id"), 2:ncol(seDT)] == "scaled_estimate")
seDT <- seDT[-which(seDT[,1] == "gene_id"),]
seDT[1:3,1:3]
geneIDs <- strsplit(as.character(seDT[,"Hybridization-REF"]), split="\\|")
geneIDs <- sapply(geneIDs, function(x) x[[2]])
stopifnot(!duplicated(geneIDs))
# rownames(seDT) <- seDT[,"Hybridization-REF"]
rownames(seDT) <- geneIDs
seDT <- seDT[,-which(colnames(seDT) == "Hybridization-REF")]
seDT[1:3,1:3]

tmpDT <- as.data.frame(data.matrix(seDT))
cat("seDT[4,5] = ", seDT[4,5], "\n")
cat("tmpDT[4,5] = ", tmpDT[4,5], "\n")
stopifnot(as.numeric(as.character(tmpDT[4,5])) == as.numeric(as.character(seDT[4,5])))
stopifnot(is.numeric(tmpDT[4,5]))
seDT[1:3,1:3]
tmpDT[1:3,1:3]
seDT <- tmpDT

colnames(seDT) <- substr(x = colnames(seDT), start = 1, stop = 15)

all_fpkm_dt <- seDT

outFile <- "PREP_FPKM_LUAD_ALL_SAMPLES/all_luad_fpkm_dt.Rdata"
dir.create(dirname(outFile))
save(all_fpkm_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

# stopifnot(rownames(exprDT) %in% rownames(seDT))
# seDT <- seDT[rownames(exprDT),c(cond1_ID, cond2_ID)]
# stopifnot(dim(seDT) == dim(exprDT))

stopifnot(file.exists(seInFile))
stopifnot(file.exists(seFile))
stopifnot(grepl(".tar.gz$", seFile))  # the file should always exist !
stopifnot(grepl(".txt$", seInFile))  # the file should always exist !
folderToRemove <- dirname(seInFile)
stopifnot(folderToRemove != seFile)

cmd <- paste0("rm -rf ", folderToRemove)
cat(paste0(cmd, "\n"))
system(cmd)

stopifnot(file.exists(seFile))
stopifnot(! file.exists(seInFile))
stopifnot(! file.exists(folderToRemove))


stop("--ok\n")
