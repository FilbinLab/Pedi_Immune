## Helper function to identify markers, top markers
my_get_markers<- function(seurat, n_top=10, ident, excludeProgram=NULL){
  seurat@meta.data[["ident"]]<-seurat@meta.data[[ident]]
  seurat<- SetIdent(seurat, value=seurat$ident)
  if(!is.null(excludeProgram)){
    seurat<- subset(seurat, ident != excludeProgram)
  }
  all_markers<- FindAllMarkers(seurat, min.pct = 0.2, logfc.threshold = 0.5)
  markers_filtered<- all_markers[all_markers$p_val_adj<0.05,]
  top_markers<- markers_filtered %>% group_by(cluster) %>% top_n(n=n_top, wt=avg_log2FC) %>% as.data.frame()
  
  return(list(filtered_markers=markers_filtered,
              top_markers=top_markers))
}

## Helper function to run standard seurat projection from reference --> query
my_project_seurat<- function(reference_seurat, query_seurat, reduction="pcaproject",
                             reference_annotation, prediction_colNames=NULL){
  seurat_list<-list(full=reference_seurat, pd1=query_seurat)
  
  ## Transfer anchors and predict
  transferAnchors<-FindTransferAnchors(reference = seurat_list$full, query=seurat_list$pd1,
                                       reduction = reduction)
  
  predictions <- TransferData(anchorset = transferAnchors, 
                              refdata =reference_seurat@meta.data[[reference_annotation]],
                              weight.reduction = reduction)
  
  predictions_final<- predictions$predicted.id; names(predictions_final)<- rownames(predictions)
  
  ## Add to query seurat object
  if(is.null(prediction_colNames)){
    query_seurat<- AddMetaData(query_seurat, predictions, col.name=c("Program_projected", 
                                                                     paste0(gsub("prediction.score.", "",
                                                                                 colnames(predictions)[2:ncol(predictions)]),
                                                                            "_Programscore_projected")))
  }else{
    query_seurat<- AddMetaData(query_seurat, predictions, col.name=prediction_colNames)
  }
  return(query_seurat)
  
  
}


## Helper for barchart by program proportion
## grouping variable is x axis, coloring variable is fill, wrapping variable (optional) is for facet wrap (FALSE or name of column)
my_barchart_programProp<- function(seurat, grouping_variable="sample", coloring_variable, wrapping_variable=FALSE,
                                   colors){
  df<- as.data.frame(table(seurat@meta.data[[grouping_variable]], seurat@meta.data[[coloring_variable]]))
  colnames(df)<- c("Sample", "Program", "NCells")
  df<- df %>% group_by(Sample) %>% mutate(perCells=NCells/sum(NCells)) %>% as.data.frame()
  
  if(wrapping_variable != FALSE){
    df$wrap<- plyr::mapvalues(df$Sample,seurat@meta.data[[grouping_variable]], seurat@meta.data[[wrapping_variable]], warn_missing = FALSE)
  }
  p<- ggplot(df, aes(x=Sample, y=perCells, fill=Program))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=colors[names(colors) %in% df$Program])+
    theme_classic()
  if(wrapping_variable !=FALSE){
    p+facet_grid(cols=vars(wrap), scales = "free_x", space="free_x")
  }else{p}
}
## Order features for plotting
## input should be of the format from seurat's FindAllMarkers
my_orderFeatures_forDotplot<- function(markers){
  markers$cluster<- as.character(markers$cluster)
  markers$cluster<- factor(markers$cluster, 
                           levels=unique(markers$cluster)[order(unique(markers$cluster))])
  markers$gene<- factor(markers$gene, levels=markers$gene[order(markers$cluster)])
  markers<- markers[order(markers$gene),]
  return(markers)
}



## Helper function for dotplot of markers of interest
my_dotplot<- function(seurat, grouping_variable, features){
  seurat@meta.data[["tmp"]]<- seurat@meta.data[[grouping_variable]]
  DotPlot(seurat, group.by = "tmp", features=features)+
    coord_flip()+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    scale_color_gradient2(low="blue", mid="white", high="red")
}

## Helper for plotting histogram of max score
## Helps to determine reasonable "Score too low" threshold
my_histMaxScores<- function(seurat, grouping_variable,score_variable, colors){
  seurat@meta.data[["grouping_variable"]]<- seurat@meta.data[[grouping_variable]]
  seurat@meta.data[["score_variable"]]<- seurat@meta.data[[score_variable]]
  cowplot::plot_grid(plotlist=lapply(unique(seurat$grouping_variable), function(x){
    tmp<- subset(seurat, grouping_variable==x)
    ggplot(tmp@meta.data, aes(x=score_variable, fill=grouping_variable))+
      geom_histogram()+
      scale_fill_manual(values=colors[x])+
      theme_classic()+
      theme(legend.position = "none")+
      ggtitle(x)
  }))
}

## function for sankey plot
my_sankeyPlot<- function(bd, colors){
  ggplot(bd,
         aes(y = NumberOfCells, axis1 = OriginalAnnotation , axis2 = ProjectedAnnotation)) +
    geom_alluvium(aes(fill = OriginalAnnotation), width = 1/12) +
    scale_fill_manual(values=colors)+
    theme_classic()+
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    theme(axis.text = element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          legend.position = "none")+
    ylab("")
}

## Helper function to map query (pd1 treated, 10x) seurat object onto same space as reference (full cohort) seurat object
## Note that this assumes that the reference seurat object is harmony integrated- harmony reduction used for re-running umap (necessary for returning model) + mapping query
my_mapSeurat<- function(ref_seurat, query_seurat, metadata_column_to_project){
  anchors <- FindTransferAnchors(reference = ref_seurat, 
                                 query = query_seurat,
                                 dims = 1:30, 
                                 reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, 
                              refdata =ref_seurat@meta.data[[metadata_column_to_project]],
                              dims = 1:30)
  query_seurat <- AddMetaData(query_seurat, metadata = predictions)
  query_seurat$Final_Annot<- query_seurat$Program_projected
  ref_seurat <- RunUMAP(ref_seurat, dims = 1:20,n.neighbors = 30, reduction = "harmony", return.model = TRUE)
  refdata_list<- list(metadata_column_to_project); names(refdata_list)<- metadata_column_to_project
  query_seurat <- MapQuery(anchorset = anchors, 
                           reference = ref_seurat, 
                           query = query_seurat,
                           refdata = refdata_list, 
                           reference.reduction = "harmony", 
                           reduction.model = "umap")
  return(query_seurat)
}

## Helper function to merge 2 seurat objects + integrate (seurat or harmony) if needed + rerun clustering  + add metadata back
## Integration options: None, Harmony, Seurat
my_mergeRerunClustering<- function(seurat1, seurat2, Integration, integrateBy="sample"){ 
  
  ## For running seurat's integration
  if(Integration=="Seurat"){
    seurat_list<- list(seurat1, seurat2)
    ifnb.list <- lapply(X = seurat_list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = seurat_list)
    anchors <- FindIntegrationAnchors(object.list =seurat_list, anchor.features = features)
    min_n<- min(ncol(seurat1), ncol(seurat2))
    k.weight<- ifelse(min_n-1 > 100, 100, min_n-1)
    seurat_merge <- IntegrateData(anchorset = anchors, k.weight = k.weight)
    
    DefaultAssay(seurat_merge) <- "integrated"
    
    # Run the standard workflow for visualization and clustering
    seurat_merge <- ScaleData(seurat_merge, verbose = FALSE) %>% 
      RunPCA(npcs = 30, verbose = FALSE) %>%
      RunUMAP(reduction = "pca", dims = 1:30)%>%
      FindNeighbors(reduction = "pca", dims = 1:30)%>%
      FindClusters( resolution = 0.5)
    
  } else{ ## For harmony or no integration
    RunHarmony<- ifelse(Integration=="Harmony", TRUE, FALSE)
    seurat_merge<- merge(seurat1, seurat2)
    seurat_merge_meta<- seurat_merge@meta.data
    seurat_merge<- RunFullSeurat(seurat_merge@assays$RNA@counts, RunHarmony = RunHarmony, sample=seurat_merge_meta[[integrateBy]])
    seurat_merge<- AddMetaData(seurat_merge, seurat_merge_meta[,!colnames(seurat_merge_meta) %in% colnames(seurat_merge@meta.data)])
  }
  
  
  return(seurat_merge)
}

## Helper function to normalize, score, and add scores/max program to seurat
my_scoreAddMax<- function(seurat, markers){
  cm_list<- NormCenter(seurat@assays$RNA@counts)
  cm_mean<- rowMeans(log2(cm_list$raw_data + 1))
  scores<- as.data.frame(t(scoreNmfGenes(cm_list$center_data, cm_mean, markers)))
  scores$max_Program<- apply(scores, 1, function(x){names(x)[which.max(x)]})
  scores$max_Programscore<- apply(scores[,colnames(scores) != "max_Program"], 1, function(x){x[which.max(x)]})
  colnames(scores)<- paste0(colnames(scores), "_scored")
  seurat<- AddMetaData(seurat, scores)
  return(seurat)
}