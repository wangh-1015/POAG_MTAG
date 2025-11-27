# Fig. 2: GWAS manhattan plot
# "POAG_MTAG.txt" is just a sample data (SNP with P<5e-8)
# note: the .pdf fig will be saved in "your_directory"  
install.packages("CMplot")
install.packages("data.table")
install.packages("dplyr")

library(CMplot); library(data.table); library(dplyr)
setwd("your_directory")
MTAG <- fread("code_package/data/POAG_MTAG.txt", sep = "\t", header = T)

CMplot(MTAG,
       col = c("#82e0aa","#5dade2"),
       type="p",
       LOG10 = TRUE,
       plot.type=c("m"),
#       points.alpha=100,
       cex=0.25,
       band=0.2,
       ylim=c(0,80),
       threshold = 5e-8,
       threshold.col="#f1948a",
       threshold.lwd = 1, 
       threshold.lty=1,
       amplify=FALSE,
#       highlight=SNP_list[,1],
#       highlight.cex= 0.25,
#       highlight.col="#a569bd",
#       highlight.text=SNP_list[,3],
#       highlight.text.col="#a569bd",
#       highlight.text.cex=0.27,
       file.output=TRUE,
       file.name="POAG",
       file="pdf",
       width=16, 
       height=10,
       chr.labels.angle = -30,
       verbose=FALSE)
