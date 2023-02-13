# Compute regulon specificity score (RSS) 

## Shannon Entropy
shannon_entropy <- function(x){
  x = x[x>0]
  return(-1 * sum(x*log(x)))
}

## JSD and RSS for one cell type and one regulon
RSS <- function(x,y){
  JSD = shannon_entropy((x+y)/2) - (shannon_entropy(x)+shannon_entropy(y))/2
  RSS = 1 - sqrt(JSD)
  return(RSS)
}

## Compute RSS scores for all cell types and regulons 
RSS_vector <- function(a=RAS_norm, b=annot_onehot){
  RSS_res = NULL
  for (i in colnames(a)){
    tmp = NULL
    for (j in colnames(b)){
      tmp = c(tmp, RSS(a[,i], b[,j]))
    }
    RSS_res = cbind(RSS_res,tmp)
  }
  RSS_res = t(RSS_res)
  rownames(RSS_res) = colnames(a)
  colnames(RSS_res) = colnames(b)
  return(RSS_res)
}

# Compute top N specific/active regulons
## @param df, a data frame of aggregated regulon activities, with row = regulons; column = cell types
## @param x, target cell type 
## @param N, top N regulons to return 
topN_regulons <- function(df, x, N=10,
                          col_names = c("Regulon", "Metaprogram", "Activity")){
  tmp = df[order(df[,x], decreasing = T),x, drop=F]
  res = cbind.data.frame(rownames(tmp), rep(x, ncol(tmp)), tmp[,1])
  colnames(res) = col_names 
  return(res[1:N,])
  ##return(rownames(df[order(df[,x], decreasing = T),x, drop=F])[1:N])
}

# Plot heatmap of top N regulons
## @param mat, matrix of aggregated regulon activitiess (row = regulons; column = cell types)
## @param df, df with top N regulons for each cell type (output from topN_regulons)
## @param plot_type, type of top N regulons, usually either specific or active
## @param regulator_names, names of each cell type
## @param fig_dir, directory for output figures 
## @param palette, color palette for heatmap 
## @param num_top, total number of top regulons to plot 
plotTopRegulonHeatmap <- function(mat, df, plot_type, 
                                  regulator_names = unique(cellInfo$CellType2),
                                  fig_dir = fig_folder,
                                  palette = "RdBu", num_top=10, 
                                  cexRow = 2, cexCol = 1.5){
  top_regulons = df
  
  ## Subset top %num_top% regulons 
  top_regulons_selected = lapply(regulator_names, 
                                 
                                 function(x) as.character(top_regulons[top_regulons$Metaprogram
                                                                       == x, "Regulon"]))
  
  top_regulons_selected = lapply(top_regulons_selected, function(x) x[!is.na(x)])
  names(top_regulons_selected) = regulator_names
  
  ## Make sidebar for the heatmap
  sidebar_row = factor(unlist(sapply(names(top_regulons_selected), function(x) rep(x, length(top_regulons_selected[[x]])))))
  names(sidebar_row) = NULL
  
  ## Heatmap
  ## Need to update palette
  if (palette == "RdBu"){
    colors = colorRampPalette(rev(brewer.pal(n = 9, name = palette)))(100)
  } else{
    colors = colorRampPalette(brewer.pal(n = 9, name = palette))(100)
  }
  png(paste0(fig_dir, paste0("top_regulons_", plot_type, "_heatmap.png")), 
      width=1600, height=1200)
  aheatmap(mat[unlist(top_regulons_selected), names(top_regulons_selected)],
           color = colors,
           Rowv = NA, Colv = NA, annRow = sidebar_row, 
           cexAnn = 3, fontsize = 16, cexRow = cexRow, cexCol = cexCol,
           width=16, height=12)
  dev.off()
}
