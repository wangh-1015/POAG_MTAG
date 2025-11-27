# Fig. 4: TWAS, PWAS manhattan plot, Locuszoom plot
# note: the .pdf fig will be saved in "your_directory"  

library(CMplot); library(data.table); library(dplyr)
setwd("your_directory")

TWAS <- fread("code_package/data/TWAS.txt", sep = "\t", header = T)
TWAS_highlight <- fread("code_package/data/TWAS_highlight.txt", sep = "\t", header = T)

PWAS <- fread("code_package/data/PWAS.txt", sep = "\t", header = T)
PWAS_highlight <- fread("code_package/data/PWAS_highlight.txt", sep = "\t", header = T)

# Fig. 4A TWAS plot
CMplot(TWAS,
       col = c("#FFBB78","#17BECF"),
       type="p",
       LOG10 = FALSE,
       plot.type="m",
       cex=0.4,
       lab.cex=1,
       ylab = "TWAS Z score",
       ylim=c(-15,15),
       band=0.2,
       threshold = c(4.68,-4.68),
       threshold.col="#f1948a",
       threshold.lwd = 1, 
       threshold.lty=2,
       amplify=F,
       highlight=TWAS_highlight[,1],
       highlight.cex= 0.5,
       highlight.col="#a569bd",
       highlight.text=TWAS_highlight[,1],
       highlight.text.col="#a569bd",
       highlight.text.cex=0.8,
       file.output=TRUE,
       file.name="TWAS plot",
       file="pdf",
       width=16, 
       height=10,
       chr.labels.angle = -30,
       verbose=FALSE)

# Fig. 4B PWAS plot
CMplot(PWAS,
       col = c("#FFBB78","#17BECF"),
       type="p",
       LOG10 = FALSE,
       plot.type="m",
       cex=0.4,
       lab.cex=1,
       ylab = "PWAS Z score",
       ylim=c(-11.1,11.1),
       band=0.2,
       threshold = c(4.12,-4.12),
       threshold.col="#f1948a",
       threshold.lwd = 1, 
       threshold.lty=2,
       amplify=FALSE,
       highlight=PWAS_highlight[,1],
       highlight.cex= 0.5,
       highlight.col="#a569bd",
       highlight.text=PWAS_highlight[,1],
       highlight.text.col="#a569bd",
       highlight.text.cex=0.8,
       file.output=TRUE,
       file.name="PWAS plot",
       file="pdf",
       width=16, 
       height=10,
       chr.labels.angle = -30,
       verbose=FALSE)


### Fig. 4C locuszoom plot, just one sample for gene: PIK3R3 
# read MTAG summary
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ensembldb")
BiocManager::install("AnnotationHub")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("GenomicFeatures")
install.packages("locuszoomr")

library(data.table)
library(locuszoomr)
library(AnnotationHub)
library(EnsDb.Hsapiens.v86)

MTAG <- fread("code_package/data/locuszoom_PIK3R3.txt", sep = "\t", header = T)
gene <- "PIK3R3"

# locus defination
# you need to get token first in "Obtaining LD information" part from website:
# "https://cran.r-project.org/web/packages/locuszoomr/vignettes/locuszoomr.html"

locus <- locus(data = MTAG, gene = gene, flank = 1e5,
              ens_db = "EnsDb.Hsapiens.v86",
              chrom="CHR",pos="BP",p="P")
locus <- link_LD(locus, token = "fe80a6e3fc78") 

# plot
pdf("PIK3R3.pdf", width = 9, height = 6)
locus_plot(locus,
             legend_pos = "topright",
             labels = c("index"),
             highlight = gene)
             
dev.off()