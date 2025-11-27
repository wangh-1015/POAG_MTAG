# gchromVAR
# more detailed instruction could be found in website:"https://github.com/caleblareau/gchromVAR"
# install
devtools::install_github("caleblareau/gchromVAR")

library(chromVAR)
library(gchromVAR)
library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
set.seed(123)
register(MulticoreParam(2))

### 1. import data to create RangedSummarizedExperiment
mat <- getMatrixFromProject(ArchRProj = proj_retina, useMatrix = "PeakMatrix")
saveRDS(mat, file = "PeakMatrix_SCE.rds")

counts <- assay(mat)
peaks <- rowRanges(mat)
cell_meta <- colData(mat)
celltype <- getCellColData(proj_retina, select = "Clusters2")$Clusters2
names(celltype) <- proj_retina$cellNames
celltypes <- unique(celltype)

# group according to cell type
grouped_counts <- sapply(celltypes, function(ct) {
  cells <- names(celltype)[celltype == ct]
  if(length(cells) == 1){
    count_mat[, cells]
  } else {
    Matrix::rowSums(count_mat[, cells, drop = FALSE])
  }
})
grouped_counts <- as.matrix(grouped_counts)

col_meta <- DataFrame(
  names = colnames(counts),
  cellType = proj_retina$Clusters2, 
  sample = proj_retina$Sample
)

SE <- SummarizedExperiment(
  assays = list(counts = counts),
  rowData = peaks,
  colData = col_meta
)
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)

### 2. Importing GWAS summary statistics
finemap_file <- "/home/wangh/project1_glaucoma/analysis/gchromVAR/finemap_pip"

# use importBedScore to map peak region
ukbb <- importBedScore(rowRanges(SE), files = finemap_file, colidx = 5)

### 3. calculate gchromVAR score
ukbb_wDEV <- computeWeightedDeviations(SE, ukbb)
zdf <- reshape2::melt(t(assays(ukbb_wDEV)[["z"]]))
zdf[,2] <- gsub("_PP001", "", zdf[,2])
colnames(zdf) <- c("ct", "tr", "Zscore")
head(zdf)

# Bonferroni test 
zdf$gchromVAR_pvalue <- pnorm(zdf$Zscore, lower.tail = FALSE)
zdf[zdf$gchromVAR_pvalue < 0.05/13, ]

