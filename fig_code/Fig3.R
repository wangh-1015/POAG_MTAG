# Fig. 3: common and rare variants converge on same genes
# note1: the .pdf fig will be saved in "your_directory" 
# note2: There's a little difference between this figure and fig 
# in our manuscript. We have used Adobe illustrator to re-draw it.
library(data.table); library(dplyr); library(ggplot2)
setwd("your_directory")

df <- fread("code_package/data/common_rare.txt", sep = "\t", header = T)

### Fig. 3A
pdf("rare_common.pdf")

ggplot(df, aes(x = MAF, y = BETA)) +
  
  geom_point(aes(color = study), size = 3.2) +
  
  scale_x_continuous(trans = 'log10',
                     breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1),
                     labels = scales::trans_format('log10', scales::math_format(10^.x))) +
  scale_y_continuous(breaks = seq(-1.5, 5.5, by = 1)) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1.3) +
  
  scale_color_manual(name = "study",
                     values = c("MTAG 5e-8" = "grey",
                                "WGS rare variants 5e-4" = "black", 
                                "WGS rare variants 1e-7" = "red",
                                "WGS rare variants 1e-5" = "orange",
                                "common genes" = "purple")) +
  
  guides(color = guide_legend(override.aes = list(size = 6,shape = 16))) +
  
  theme_minimal(base_size = 16) +
  
  labs(y = 'Minor allele effect size(beta)', color = 'Group') +

  # P<1e-7
  geom_point(data = subset(df, gene == "MYOC"), 
             aes(x = MAF, y = BETA),
             size = 7, shape = 21, color = "black", fill = "red") +

  geom_text(data = subset(df, gene == "MYOC"),
            aes(label = gene), 
            vjust = -1, hjust = 0.5, size = 5, color = "red") +

  # P<1e-5
  geom_point(data = subset(df, gene %in% c("ANGPTL7","ANKFY1","ADAM9")), 
             aes(x = MAF, y = BETA),
             size = 5, shape = 21, color = "black", fill = "orange") + 

  geom_text(data = subset(df, gene %in% c("ANGPTL7","ANKFY1","ADAM9")),
            aes(label = gene), 
            vjust = -1, hjust = 0.5, size = 5, color = "orange") +
  
  # MTAG 5e-8 
  geom_point(data = subset(df, gene %in% c("TMCO1","AFAP1")), 
             aes(x = MAF, y = BETA),
             size = 4, shape = 21, color = "black", fill = "grey") +
  
  geom_text(data = subset(df, gene %in% c("TMCO1","AFAP1")),
            aes(label = gene), 
            vjust = -1, hjust = 0.5, size = 5, color = "black") +

  # common gene examples
  geom_point(data = subset(df, gene %in% c("NXN","TMEM151A","AKAP6") & MAF < 0.01), 
             aes(x = MAF, y = BETA),
             size = 4, shape = 21, color = "black", fill = "purple") +
  
  geom_text(data = subset(df, gene %in% c("NXN","TMEM151A","AKAP6") & MAF < 0.01),
            aes(label = gene), 
            vjust = -1, hjust = 0.5, size = 5, color = "purple") +
  
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    aspect.ratio = 1
  )

dev.off()

### Fig. 3B burden-gene finemap_methods converged genes

a <- data.frame(method=c("MAGMA","FLAMES","Nearest gene","PoPS","cS2G","TWAS","Retina eQTLcMR"),
                num_overmap=c(28,17,14,13,12,12,5))

desired_order <- c("MAGMA","FLAMES","Nearest gene","PoPS","cS2G","TWAS","Retina eQTLcMR")
a$method <- factor(a$method, levels = desired_order)

a_layer1 <- a
a_layer1$Layer <- "Burden-gene only"
a_layer1$y <- 46

a_layer2 <- a
a_layer2$Layer <- "Shared gene"
a_layer2$y <- a$num_overmap

a_combined <- rbind(a_layer1, a_layer2)

pdf("burden_finemap_overlap.pdf")

ggplot(a_combined, aes(x=method)) +
  # bar for burden
  geom_bar(data = a_layer1, aes(y = y, fill = Layer), stat = "identity") +
  
  # bar for finemap
  geom_bar(data = a_layer2, aes(y = y, fill = Layer), stat = "identity") +
  
  scale_fill_manual(values = c("Burden-gene only" = "#AEC7E8", "Shared gene" = "#FFBB78")) +
  
  labs(x=" ", y="# burden-genes in WGS", fill = NULL) +
  theme_minimal() +
  theme(
    legend.position.inside = c(0.8, 0.8),
    legend.text = element_text(size = 14), 
    axis.text.x = element_text(angle=45, size = 12,color = "black"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16, vjust = 2, hjust = 0.5, margin = margin(r = 0)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )

dev.off()

