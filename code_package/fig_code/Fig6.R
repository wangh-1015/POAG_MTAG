# Fig. 6
# note: the .pdf fig will be saved in "your_directory"  

# Fig. 6A/B
# note: only choose panel B for example

library(CMplot); library(data.table); library(dplyr)
setwd("your_directory")

CM_plot <- fread("code_package/data/CM_plot_final_2M.txt.txt", sep = "\t", header = T)

View(CM_plot)


CM_plot$cell_type <- factor(CM_plot$cell_type,
                            levels = c("Pigment epithelial", "astrocyte", "Müller glia", "retina ganglion", 
                                       "starburst amacrine","microglia","GABA amacrine","AII-amacrine",
                                       "gly-amacrine","rod bipolar","cone", "Rod","OFF-cone bipolar", 
                                       "horizontal","ON-cone bipolar")
)
CM_plot <- CM_plot[order(CM_plot$cell_type), ]

# cell_type filter for 2M data
CM_plot <- CM_plot %>%
  filter(cell_type %in% c("astrocyte","microglia","Müller glia","AII-amacrine",
                          "GABA amacrine","cone","gly-amacrine","starburst amacrine",
                          "horizontal","Pigment epithelial","retina ganglion",
                          "OFF-cone bipolar","ON-cone bipolar",
                          "Rod","rod bipolar"))
# 15 colors for 2M
cell_colors <- c(
  "#D62728", "#1F77B4", "#2CA02C", "#F7B6D2", "#FF4B4B",
  "#FF7F0E", "#9467BD", "#17BECF", "#8C564B", "#E377C2",
  "#AEC7E8", "#BCBD22", "#98DF8A", "#FFBB78", "#7F7F7F")

CMplot(CM_plot,
       col = cell_colors,
       type="p",
       pch= 6,
       LOG10 = FALSE,
       plot.type="m",
       points.alpha=90,
       cex=0.8,
       ylab = "-log10_adj.P",
       band=0.8,
       mar=c(6,6,3,6),
       threshold = 2.653647,
       threshold.col="#f1948a",
       threshold.lwd = 2, 
       threshold.lty=2,
       amplify=FALSE,
       file.output=TRUE,
       file.name="scDRS plot 1M",
       file="pdf",
       dpi=300,
       width=12, 
       height=7,
       chr.labels.angle = -15,
       verbose=TRUE)

# Fig6. C
df <- data.frame(count=c(12034,3740,6980,4657,4498,5781),
                 group= c("0", "1", "2-5", "6-10", "11-20", ">20"))
df$group <- factor(
  df$group,
  levels = c("0", "1", "2-5", "6-10", "11-20", ">20")
)

pdf("Fig6C.pdf")
ggplot(df, aes(x = group, y = count, fill = "#17BECF")) +
  geom_bar(stat = "identity", 
           width = 0.6,       # Increase bar width to decrease gaps
           alpha = 0.8        # Set transparency
  ) +
  theme_minimal() +  # Minimal theme for clean background
  labs(title = "n = 37690 promoter peaks", 
       x = "# of co-accessible peaks",
       y = "# of promoter peaks") +
  theme(legend.position = "none", # Remove legend
        plot.title = element_text(hjust = 0.5, vjust = -25, size = 26),  # Larger title, centered
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.4, size = 20),  # Larger x-axis labels with angle
        axis.text.y = element_text(size = 20),  # Larger y-axis labels
        axis.title.x = element_text(size = 26, vjust = 2, hjust = 0.5, margin = margin(r = 0)),
        axis.title.y = element_text(size = 26, vjust = 2, hjust = 0.5, margin = margin(r = 0)),
        panel.grid.major.x = element_blank(),  # Hide vertical grid lines
        panel.grid.minor.y = element_blank()   # Hide minor grid lines
  )
dev.off()

### Fig. 6D
df2 <- data.frame(count=c(113433,71023,17960,5446,1232,2098),
                 group= c("0", "1", "2", "3", "4", ">4"))

df2$group <- factor(
  df2$group,
  levels = c("0", "1", "2", "3", "4", ">4")
)

pdf("Fig6D.pdf")
ggplot(df2, aes(x = group, y = count, fill = "#17BECF")) +
  geom_bar(stat = "identity", 
           width = 0.6,       # Increase bar width to decrease gaps
           alpha = 0.8        # Set transparency
  ) +
  theme_minimal() +  # Minimal theme for clean background
  labs(title = "n = 97759 peaks", 
       x = "# of target genes",
       y = "# of peaks") +
  theme(legend.position = "none", # Remove legend
        plot.title = element_text(hjust = 0.5, vjust = -25, size = 26),  # Larger title, centered
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.4, size = 20),  # Larger x-axis labels with angle
        axis.text.y = element_text(size = 20),  # Larger y-axis labels
        axis.title.x = element_text(size = 26, vjust = 2, hjust = 0.5, margin = margin(r = 0)),
        axis.title.y = element_text(size = 26, vjust = 2, hjust = 0.5, margin = margin(r = 0)),
        panel.grid.major.x = element_blank(),  # Hide vertical grid lines
        panel.grid.minor.y = element_blank()   # Hide minor grid lines
  ) 
dev.off()

# Fig. 6E-F
# note: example data was not available for the large storage of scATAC data. The below was only 
# a code example
library(GenomicRanges)

# 1. target gene/promoter gene track information
custom_df <- data.frame(
  start = c(46878268),
  end = c(46939917)
  #value=c(0.5,0.605,0.637)
)

# make GRanges object for target gene/promoter
gr_loops <- GRanges(
  seqnames = rep("chr11",1),
  ranges = IRanges(start = custom_df$start, end = custom_df$end),
  #gene = rep("TFAP2B",2)
  #value = custom_df$value
)

# 2. HiChiP loop track information
custom_df <- data.frame(
  seqnames = rep("chr11",2),
  start = c(46878268,46878268),
  end = c( 47012500,47155000),
  name = c("loop1")
)

# make GRanges object for HiChiP loop
hic_loops <- GRanges(
  seqnames = rep("chr11",2),
  ranges = IRanges(start = custom_df$start, end = custom_df$end),
  gene = rep("loop",2)
  # value = custom_df$value
)


grlist<- GRangesList("target gene" = gr_loops,"H3K27ac HiChIP loop"=hic_loops)

# region=GRanges(
#  seqnames = rep("chr1"),
#  ranges = IRanges(start = 61300000, end = 61560000),)

p <- plotBrowserTrack(
  ArchRProj = proj_retina,
  groupBy = "Clusters2",
  plotSummary = c("bulkTrack", "geneTrack","loopTrack"),
  sizes = c(10, 1.5, 2),
  threads = getArchRThreads(),
  #features = getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)[c("astrocyte","Muller_glia")],
  geneSymbol = "LRP4",
  upstream = 120000,
  downstream = 150000,
  loops = grlist,
  baseSize = 13,
  facetbaseSize = 13,
  title = "LRP4",
)

plotPDF(
  plotList = if (!is.list(p)) list(p) else p,
  name = "track.pdf", 
  ArchRProj = proj_retina, 
  addDOC = FALSE, 
  width = 10, 
  height = 10
)
