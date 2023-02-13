## Function to derive pseudobulk (aggregated gene expressions for each cell subpopulation)
## Based on classification by NMF
## @para seurat_obj seurat object that stores cm
## @para nmf_score df that stores scores of each metaprogram; seurat_obj and nmf_score should have identical cells
## @para out_dir output directory to store results
## @para aggr_method method to aggregate results, either mean or median 
## @score_diff discard cells with metaprogram score less than this cutoff when aggregating 
makePseudobulk <- function(seurat_obj, nmf_score, 
                           out_dir = seurat_analysis_folder,
                           aggr_method="mean", score_cutoff=1){
    metagene_program_names = sort(unique(nmf_score$signature_1))
    pseudobulk = NULL
    
    for (metagene in metagene_program_names){
        cm = seurat_obj@raw.data[,nmf_score$signature_1 == metagene & 
                                     nmf_score$score_1 >= score_cutoff]
        if (aggr_method == "mean"){
            pseudobulk = cbind(pseudobulk, rowMeans(cm))
        }
        if (aggr_method == "median"){
            pseudobulk = cbind(pseudobulk, rowMedians(cm))
        }
    }
    colnames(pseudobulk) = metagene_program_names
    rownames(pseudobulk) = rownames(seurat_obj@raw.data)
    saveRDS(pseudobulk, file=paste0(out_dir, "pseudobulk_", aggr_method, ".RDS"))
}

## Run UMAP using 3 components instead of the standard 2
## Used for 3D UMAP plotting
## Defaul nDims set to 20, may want to adjust based on JackStraw
## Seurat object should already be clustered before running
Run3DUMAP<- function(seuratObj_clustered, nDims=20){
    # Re-run UMAPs that you have accurate calculations for all UMAP(s)
    yourseuratobject <- RunUMAP(seuratObj_clustered,
                                dims = 1:nDims,
                                n.components = 3L)
    return(yourseuratobject)
    
}

## Plot the 3 component UMAP object created above
## plotting performed using plotly- creates 3D html output that can be viewed from multiple angles
## Can save a static image from the html output if desired
## Can specify grouping variable- seurat clusters, sample, etc
Plot3DUMAP<-function(seurat_obj_3DUMAP,
                     outputDir=paste0(working_dir, "/", seurat_fig_folder, "3D_umap.html"),
                     GroupBy="seurat_clusters"){
    # Load plot_ly
    library(plotly)
    
    # Extract tSNE information from Seurat Object
    umap_1 <- seurat_obj_3DUMAP[["umap"]]@cell.embeddings[,1]
    umap_2 <- seurat_obj_3DUMAP[["umap"]]@cell.embeddings[,2]
    umap_3 <- seurat_obj_3DUMAP[["umap"]]@cell.embeddings[,3]
    
    
    # Prepare a dataframe for cell plotting
    plot.data <- FetchData(object = seurat_obj_3DUMAP, 
                           vars = c("UMAP_1", "UMAP_2", "UMAP_3", GroupBy))
    colnames(plot.data)<-c("UMAP_1", "UMAP_2","UMAP_3", "GroupBy")
    
    # Make a column of row name identities (these will be your cell/barcode names)
    plot.data$label <- paste(rownames(plot.data))
    
    # Plot your data, in this example my Seurat object had 21 clusters (0-20)
    p<-plot_ly(data = plot.data, 
               x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
               color = ~GroupBy, 
               colors = c("lightseagreen","gray50","darkgreen","red4","red","turquoise4","black",
                          "yellow4","royalblue1","lightcyan3","peachpuff3","khaki3","gray20",
                          "orange2","royalblue4","yellow3","gray80","darkorchid1","lawngreen",
                          "plum2","darkmagenta")[1:length(unique(plot.data[,4]))],
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 5, width=2), 
               text=~label, 
               hoverinfo="text")
    
    htmlwidgets::saveWidget(p, file=outputDir)
}
