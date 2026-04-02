# Fig. 5: POAG genes prioritization
# note1: the .pdf fig will be saved in "your_directory" 
# note1: There's a little difference between this figure and fig 
# in our manuscript. We have used Adobe illustrator to re-draw it.

BiocManager::install("ComplexHeatmap") 
library(circlize)
library(ComplexHeatmap)
library(data.table); library(dplyr)
setwd("your_directory")

top50_genes <- fread("code_package/data/top50_genes.txt", sep = "\t", header = T)
top50_genes <- as.data.frame(top50_genes)
top50_genes <- top50_genes[1:50,-1]
rownames(top50_genes) <- top50_genes[,1]
top50_genes <- top50_genes[,-1]

###
col_fun = colorRamp2(c(0,1), c("white", "white"))

#column_split
column_split <- rep(1:6, times = c(9, 2, 7, 3, 2, 1))

pdf("top50_genes.pdf", width = 14, height = 14)

Heatmap(top50_genes, 
        #basic
        name = "POAG genes priority",
        col = col_fun,
        na_col = "black",
        border_gp = gpar(col = "black", lwy = 1.5), rect_gp = gpar(col = "white", lwd = 6),
        column_title = "POAG genes priority", column_names_rot = 70,
        cluster_rows = FALSE, show_column_dend = FALSE,
        show_heatmap_legend = FALSE,
        cluster_columns = FALSE,
        
        # row column
        row_names_side = "left", row_names_gp = gpar(fontsize = 8),
        column_names_side = "top", column_names_gp = gpar(fontsize = 8),
        
        #split columns/rows
        column_split = column_split, column_gap = unit(7, "mm"),
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#99ea99", "#ea9999", "#9999ea", "#eaea99", "grey80","#e3ccf9")),
                                                            labels = c("Genetic analysis", "Integrated", "cS2G", "single-cell","Score","Novel"), 
                                                            labels_gp = gpar(col = "black", fontsize = 9))),
        
        #layer control
        layer_fun = function(j, i, x, y, width, height, fill) {
          n_row = nrow(top50_genes)
          n_col = ncol(top50_genes)
          
          for(k in seq_along(i)) {
            cell_value = top50_genes[i[k], j[k]]
            
            if(j[k] >21 && j[k] < 24 ) {
              if(cell_value >= 13) {
                color = "grey80"
                
                grid.rect(
                  x = x[k] + unit(0.35, "mm"),
                  y = y[k],
                  width = width[k] * 0.9,
                  height = height[k] * 0.75,
                  gp = gpar(col = color, fill = color)
                )
                
                grid.text(
                  label = as.character(cell_value), 
                  x = x[k] + unit(0.55, "mm"), 
                  y = y[k], 
                  gp = gpar(col = "black", fontsize = 8)
                )
              }
              else{
              color = "white"
              grid.rect(
                x = x[k] + unit(0.1, "mm"),   
                y = y[k],   
                width = width[k] * 0.9,       
                height = height[k] * 0.75,     
                gp = gpar(col = color, fill = color)
              )
              
              grid.text(
                label = as.character(cell_value), 
                x = x[k] + unit(0.55, "mm"), 
                y = y[k], 
                gp = gpar(col = "black", fontsize = 8)
              )
              }
            }
            else {
              if(!is.na(cell_value) && cell_value != 0) {
                col_group = column_split[j[k]]
                color = switch(col_group,
                               "1" = "#99ea99",
                               "2" = "#ea9999",
                               "3" = "#9999ea",
                               "4" = "#eaea99",
                               "5" = "white",
                               "6" = "#e3ccf9")
                
                grid.rect(
                  x = x[k] + unit(0.1, "mm"),
                  y = y[k],
                  width = width[k] * 0.9,
                  height = height[k] * 0.75,
                  gp = gpar(col = color, fill = color)
                )
              }
          }
        }
      }
        ,
        width = unit(19, "cm"), height = unit(24, "cm")
)

dev.off()
