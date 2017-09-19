###
# Based on the sourcecode of plot_pseudotime_heatmap from the monocle package
# Author: Max Kaufmann, MD, 30th August 2017

GroupedPseudotimeHeatmap <- function(monoc.object, path.to.genelists, markergene = "Bsn_rat", pdf.folder = "./pdf"){
  
  library(monocle)
  library(RColorBrewer)
  dir.create(pdf.folder, showWarnings = FALSE)
  
  # load all genesets into one list, remove unexpressed genes and remove marker gene
  geneset.list <- list()
  file.names <- list.files(path.to.genelists)
  for (k in 1:length(file.names)) {
    geneset.list[[k]] <- scan(paste0(path.to.genelists, "/", file.names[k]), what="character", sep=NULL)
    geneset.list[[k]] <- geneset.list[[k]][geneset.list[[k]] %in% rownames(monoc.object)] 
    geneset.list[[k]] <- geneset.list[[k]][geneset.list[[k]]!=markergene] 
  }
  names(geneset.list) <- gsub(".txt","", file.names)
  
  # Save gap values for heatmap
  gaps <- cumsum(lengths(geneset.list))
  gaps <- c(1, gaps)
  
  # Create pseudotime heatmap matrix (based on the plot_pseudotime_heatmap from monocle)
  CreatePseudotimeMatrix <- function(input.m) {
    library(monocle)
    scale_max <- 3
    scale_min <- -3
    pseudocount <- NA
    newdata <- data.frame(Pseudotime = seq(min(pData(input.m)$Pseudotime), 
                                           max(pData(input.m)$Pseudotime), length.out = 100))
    m <- genSmoothCurves(input.m,
                       cores = 1,
                       trend_formula = "~sm.ns(Pseudotime, df=3)", 
                       relative_expr = T, 
                       new_data = newdata)
    m = m[!apply(m, 1, sum) == 0, ]
    norm_method <- c("vstExprs", "log")
    if (norm_method == "vstExprs" && is.null(input.m@dispFitInfo[["blind"]]$disp_func) == FALSE) {
      m = vstExprs(input.m, expr_matrix = m)
  } else if (norm_method == "log") {
      m = log10(m + pseudocount)
  }
    m = m[!apply(m, 1, sd) == 0, ]
    m = Matrix::t(scale(Matrix::t(m), center = TRUE))
    m = m[is.na(row.names(m)) == FALSE, ]
    m[is.nan(m)] = 0
    m[m > scale_max] = scale_max
    m[m < scale_min] = scale_min
  return(m)
  }
  
  heatmap_matrix <- CreatePseudotimeMatrix(monoc.object[unlist(geneset.list),])
  marker_matrix <- CreatePseudotimeMatrix(monoc.object[c(markergene, markergene),])
  
  # Plot heatmaps
  library(pheatmap)
  library(colorRamps)
  library(cowplot)
  
  bks <- seq(-3, 3, by = 0.01)
  hmcols <- colorRampPalette(c("blue","white", "red"))(length(bks) - 1)
  markercols <- colorRampPalette(c("#ffffff", "#006d2c"))(length(bks) - 1)
  
  my.plots <- vector(mode='list')
  
  for(k in 1:length(geneset.list)) {
    heatmap_matrix_subs <- heatmap_matrix[gaps[k]:gaps[k+1],]
    
    # compute distance
    dist <- as.dist(1-cor(t(heatmap_matrix_subs))) #pearson
    
    # cluster distance
    hr <- hclust(d = dist, method = "ward.D")
    
    #reorder matrix
    heatmap_matrix_subs <- heatmap_matrix_subs[hr$order,]
    
  pheatmap(
    mat = heatmap_matrix_subs,
    #
    scale = "row",
    color = hmcols,
    #breaks = bks,
    border_color = NA,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 3,
    cellheight = 100 / nrow(heatmap_matrix_subs),
    cellwidth = 600 / ncol(heatmap_matrix_subs),
    legend = FALSE,
    annotation_legend = FALSE,
    annotation_names_row = FALSE,
    annotation_colors = ann_colors,
    main = names(geneset.list)[k]
  )
    my.plots[[k]] <- recordPlot()
  }
  
# Create PDF with collected plots
  
  pdf(paste0(pdf.folder, "/marker.pdf"), onefile=FALSE, width = 11.69, height = 8.27)
  pheatmap(
    mat = marker_matrix,
    color = markercols,
    border_color = NA,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_colnames = FALSE,
    show_rownames = FALSE,
    cellheight = 10,
    cellwidth = 600 / ncol(heatmap_matrix_subs),
    legend = FALSE,
    main = markergene
  )
  dev.off()
  
  for (k in 1:length(my.plots)) {
    pdf(paste0(pdf.folder, "/", names(geneset.list)[k],".pdf"), onefile=TRUE, width = 11.69, height = 8.27)
    replayPlot(my.plots[[k]])
    dev.off()
  }
}
