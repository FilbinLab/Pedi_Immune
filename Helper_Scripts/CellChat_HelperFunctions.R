
## Used to create cellchat input from seurat object for a single sample
## Input: mySample = sampleID
##        all_merged_seurat_obj = single seurat object for that sample alone
##        SpecifyCellAnnot = FALSE = just use annotations present in seurat, c("Cell1") = use only those specific cell types
##        FC/PC/pvaluethresh/onlyPos = settings for cellchat object (see cellchat vignette for details)
##        
CreateCellChatObject<-function(mySample, all_merged_seurat_obj, SpecifyCellAnnot=FALSE,
                               FCthresh=0, PCthresh=0, Pvaluethresh=0.05, 
                               only.pos=TRUE) {
  print(mySample)
  ## Subset to sample
  seurat_obj<-all_merged_seurat_obj[[mySample]]
  
  ## Set cells of interest
  ## If not specified, just select all cell types with at least 10 cells
  ## If specified, select those cell types with at least 10 cells in this sample
  if(SpecifyCellAnnot==FALSE){
    df<- table(seurat_obj$CellAnnot)
    cellTypesUse=names(df[df>10])
  }else{
    df<- table(seurat_obj$CellAnnot)
    tmp=names(df[df>10])
    cellTypesUse=SpecifyCellAnnot[SpecifyCellAnnot %in% tmp]
    }
  
  ## Subset to just these cells
  seurat_obj$CellAnnot<- as.character(seurat_obj$CellAnnot)
  seurat_obj<- subset(seurat_obj, CellAnnot %in% cellTypesUse)
  
  
  ## Create cellchat object and set database to use
  cellchat<-createCellChat(object=seurat_obj, group.by = "CellAnnot")
  cellchat@DB <- CellChatDB.human
  
  ## Downsample data for faster processing
  cellchat<-subsetData(cellchat)
  future::plan("multiprocess", workers = 4)
  
  ## preprocessing: identify OE genes/interactions
  cellchat <- identifyOverExpressedGenes(cellchat,  
                                         thresh.fc = FCthresh, 
                                         thresh.pc = PCthresh,
                                         thresh.p= Pvaluethresh,
                                         only.pos = only.pos)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  return(cellchat)
}

## Wrapper for running cellchat
## Input: cellchat = cellchat object from CreateCellChatObject
##        raw.use/population.size/min.cells = settings from cellchat- see vignette for more details
InferCellCellInteractions<- function(cellchat, raw.use=TRUE, population.size=FALSE,
                                     min.cells=10){
  message(timestamp())
  print(unique(cellchat@meta$sample))
  ## Compute the communication probability and infer cellular communication network
  message("Computing cellular communication network...")
  future::plan("multiprocess", workers = 4)
  cellchat <- computeCommunProb(cellchat, raw.use = raw.use, population.size = population.size)
  
  ## Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  
  ## Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  
  ## Calculate the aggregated cell-cell communication network
  ## sources.use/targets.use to use specific cell types
  cellchat <- aggregateNet(cellchat)
  return(cellchat)
}




## Determine appropriate figure width/height for export
DetermineFigureWidthHeight<- function(plotList, widthMultiplier, heightMultiplier){
  if(length(plotList)>=3){
    Ncol<- 3} else{
      Ncol<- length(plotList)}
  figureWidth<- Ncol * widthMultiplier
  
  if(length(plotList)<=3){
    figureHeight<- heightMultiplier} else{
      NumRows<- ceiling((length(plotList) - Ncol)/3) + 1
      figureHeight<- NumRows * heightMultiplier 
    }
  return(list(figureWidth=figureWidth,
              figureHeight=figureHeight,
              Ncol=Ncol))
}

## For plotting sig LR across samples: encode pvalue and prob in figure
## Can wrap just by sample (facetWrap="Sample")
## Or can wrap by both sample + sender cell (facetWrap=c("Sample", "sender))
## Input: sigLRs = signficiant ligand/receptor interactions to plot from Extract_res_allComparisons
##        senderCell = specific sender cell type used
##        xAxis = xaxis column- will usually be sender, but may change to sample, etc.
##        labelXAxis = FALSE = no X axis labeled
##        facetWrap = vector of length 1 or 2 with variables to wrap by
##        pathwaySubset = only plot certain pathways
##        maxProb = maximum probability limit, default to maximum present
PlotPvalueProb<- function(sigLRs, senderCell="", pt.size=16, xAxis="sender",
                          labelXAxis=FALSE, facetWrap= c("Sample", "sender"),
                          pathwaySubset="",
                          maxProb=max(sigLRs$Prob)){
  ncolForFacetWrap=length(unique(sigLRs$FacetWrap1))
 if(pathwaySubset!=""){sigLRs<-sigLRs[sigLRs$pathway %in% pathwaySubset,]}
 sigLRs$XAxis<- sigLRs[,xAxis]

 ## set upper limit for probability
 sigLRs$Prob<- ifelse(sigLRs$Prob>maxProb, maxProb, sigLRs$Prob)

 ## Dot plot of all LRs
 p<-ggplot(sigLRs, aes(x=XAxis, y=LR, size=pvalue, color=Prob))+
        geom_point(pch = 16)+
        theme_linedraw() + 
        scale_size(trans = 'reverse', range=c(1,4))+
        theme(panel.grid.major = element_blank(), 
              axis.text.x = element_blank(), 
              axis.title.x = element_blank(), 
              axis.title.y = element_blank()) + 
        scale_x_discrete(position = "bottom")+
        scale_colour_gradientn(colours = myPalette(100),limits=c(0, maxProb))
 
 ## Facetwrap- no facetwrap or single facetwrap or 2 variables
 if(length(facetWrap)==0){
   p<-p
 }
 else if(length(facetWrap)==1){
        #p<- p + facet_wrap(vars(FacetWrap1), ncol=ncolForFacetWrap, scales="free_x")+
        p<- p + facet_grid(cols=vars(FacetWrap1), scales="free_x", space="free_x")+
        ggtitle(paste0("Sender cell: ", senderCell))
 } else if(length(facetWrap)!=2){print("Can't wrap more than 2 facets!")
   }else{
     p<- p+facet_nested_wrap(vars(FacetWrap2, FacetWrap1), 
                             ncol=ncolForFacetWrap,
                             scales="free_x")
   }
 if(labelXAxis){
   p<- p+theme(axis.text.x = element_text(size=10, angle=45, 
                                          hjust=1, vjust = 1))
 }
 return(p)
}



## Helper function for merging multiple seurat objects
## sample column must be "sample"
## Cell annot column must be "CellAnnot"
MergeSeurat_RunSeuratPipeline<- function(s, seurat_list, samplesUse){
  seurat_1<- seurat_list[[1]]
  
  if (s %in% unique(seurat_1$sample) & s %in% samplesUse){
    print(paste0("Merging sample ", s))
    
    ## subset to sample, then get counts + metadata
    cm_meta_list<-lapply(1:length(seurat_list), function(x){
      seurat<- seurat_list[[x]]
      seurat<-subset(seurat, sample==s)
      
      meta<- seurat@meta.data[,c("sample", "CellAnnot")]
      meta$Type<- rep(names(seurat_list)[x], nrow(meta))
      
      cm<- as.data.frame(seurat@assays$RNA@counts)
      cm$gene<-rownames(cm)
      
      return(list(cm=cm,
                  meta=meta))})
    names(cm_meta_list)<- names(seurat_list)
    
    for(i in names(cm_meta_list)){
      print(paste0("cells in ", i, ": ", nrow(cm_meta_list[[i]]$meta)))
    }
    
    ## merge counts and metadata
    cm_list<- lapply(1:length(cm_meta_list), function(x){
      tmp<- cm_meta_list[[x]]$cm
      rownames(tmp)<- tmp$gene
      tmp<- tmp[,colnames(tmp)!="gene"]
      colnames(tmp)<- paste0(names(cm_meta_list[x]),
                             ".", colnames(tmp))
      tmp$gene<-rownames(tmp)
      return(tmp)
    })
    names(cm_list)<- names(cm_meta_list)
    
    meta_list<- lapply(1:length(cm_meta_list), function(x){
      tmp<- cm_meta_list[[x]]$meta
      rownames(tmp)<- paste0(names(cm_meta_list[x]),
                             ".", rownames(tmp))
      tmp$Cell<- rownames(tmp)
      return(tmp)
    })
    names(meta_list)<-names(cm_meta_list)
    
    
    merged_cm<- cm_list %>% purrr::reduce(left_join, by = "gene")
    rownames(merged_cm)<- merged_cm$gene
    merged_cm<-merged_cm[,colnames(merged_cm) !="gene"]
    
    merged_meta<- do.call("rbind", meta_list)
    rownames(merged_meta)<-merged_meta$Cell
    print(paste0("total cells: ", nrow(merged_meta)))
    
    merged_cm<- merged_cm[,order(colnames(merged_cm))]
    merged_meta<- merged_meta[order(rownames(merged_meta)),]
    print(paste0("Check: "))
    print(paste0("metadata/cm match: ",
                 sum(colnames(merged_cm)%in% rownames(merged_meta))))
    print(paste0("Total cells: ", nrow(merged_meta)))
    
    ## For small n: determine n dims to use downstream
    nDims<- ifelse(ncol(merged_cm)<50, ncol(merged_cm)-1, 50)
    
    ## Create seurat object from merged counts and merged meta
    seurat_merged<-preprocessSeuratObject(cm=merged_cm, project="")
    seurat_merged<- AddMetaData(seurat_merged, merged_meta)
    seurat_merged<- RunPCA(seurat_merged, features = seurat_merged@assays$RNA@var.features,
                           npcs=nDims)
    seurat_merged<-RunUMAP(seurat_merged, reduction="pca", dims = 1:nDims , 
                           n.neighbors = ifelse(nDims<30, nDims + 1, 30), verbose = FALSE)
    seurat_merged<- FindNeighbors(seurat_merged, reduction="pca", dims = 1:nDims) 
    seurat_merged<-FindClusters(seurat_merged, resolution = 0.8)
  } else{
    print("sample not found")
    seurat_merged<-NULL}
  return(seurat_merged)}

## Get significant LRs across each comparison, each sample, with pvalue + probability
## input: cellChat_obj= output from running CellChat O2 script / PathwaysRemove: any pathways you want to remove from analysis 
##        sigThreshold= pvalue thresh, generally set to 1 for now and filtered later 
##        senderCells= vector of all sender cells / receiverCells: vector of all receiver cells
##        subtype= add on subtype info if desired
Extract_res_allComparisons<- function(cellChat_obj, PathwaysRemove, sigThreshold, senderCells, receiverCells, subtypes=NULL){
  ## Cycling through each cell type comparison
  allComparisons_LRs<- lapply(cellChat_obj, function(comparisonCellChat){
    all_cellChat<- comparisonCellChat
    settings<- all_cellChat$Settings
    all_cellChat<- all_cellChat[names(all_cellChat)!= "Settings"]
    message(paste0("Identifying shared LRs in ", settings["cellTypesUsed",1]))
    if(length(comparisonCellChat)==1){return(NULL)}
    
    allSamples_sigLRs<-list()
    ## Cycling through each sample
    for (i in 1:length(all_cellChat)){
      cellChat<-all_cellChat[[i]]
      mySample<-names(all_cellChat)[i]
      print(mySample)
      
      ## Use CellChat function to extract all LRs 
      sigLRs<- subsetCommunication(cellChat, thresh = sigThreshold)
      
      ## Add sample info
      sigLRs$Sample<- mySample
      allSamples_sigLRs[[mySample]]<- sigLRs
    }
    ## Merge all samples for that comparison into 1 df
    allSamples_sigLRs_df<-do.call("rbind", allSamples_sigLRs)
    
    ## Add subtype info if needed
    if(!is.null(subtypes)){
      allSamples_sigLRs_df$Subtype<- plyr::mapvalues(allSamples_sigLRs_df$Sample, subtypes$sample, subtypes$MySubtype, 
                                                     warn_missing = FALSE)
    }
    
    return(allSamples_sigLRs_df)
  })
  
  ## Merge all comparisons into 1 df
  allComparisons_sigLRs_df<- do.call("rbind", allComparisons_LRs)
  
  ## Subset to sender/receiver cells
  if(!is.null(subtypes)){
    colnames(allComparisons_sigLRs_df)<- c("sender", "receiver", "ligand", "receptor", "Prob",
                                           "pvalue", "LR", "interaction_name_2", "pathway",
                                           "annotation", "evidence", "Sample", "Subtype")
  }else{
    colnames(allComparisons_sigLRs_df)<- c("sender", "receiver", "ligand", "receptor", "Prob",
                                           "pvalue", "LR", "interaction_name_2", "pathway",
                                           "annotation", "evidence", "Sample")
  }
  
  allComparisons_sigLRs_df<- allComparisons_sigLRs_df[allComparisons_sigLRs_df$sender %in%
                                                        senderCells &
                                                        allComparisons_sigLRs_df$receiver %in%
                                                        receiverCells,]
  
  ## remove any unwanted pathways
  if(length(PathwaysRemove)!=0){
    removeLR <- paste(PathwaysRemove,collapse="|")
    allComparisons_sigLRs_df <- allComparisons_sigLRs_df[grep(removeLR,
                                                              allComparisons_sigLRs_df$pathway, 
                                                              invert = TRUE),]
    
  }
  
  ## Remove any interactions where sender and receiver cell are the same
  allComparisons_sigLRs_df<- allComparisons_sigLRs_df[allComparisons_sigLRs_df$sender !=
                                                        allComparisons_sigLRs_df$receiver,]
  
  ## Convert all factor columns (default for CellChat subsetCommunications function) to char
  allComparisons_sigLRs_df<- allComparisons_sigLRs_df %>% 
    mutate_if(is.factor, as.character) %>% 
    as.data.frame()
  return(allComparisons_sigLRs_df)
}

## Preprocess dataframe of sig LRs for plotting sankey/river plots between cell types
## Input: data frame with following columns: "sender_major" = major cell type (cd4, tumor, etc) on ligand side
##                                           "receiver_major" = major cell type on receptor side
##                                           "sender" = minor cell type (inflammatory, ac-like, etc) on ligand side
##                                           "receiver" = minor cell type on receptor side
##                                           "LR" = ligand/receptor pair name
## major cell type for sender, major cell type for receiver (should be in sender_major and receiver_major columns)
## Normalization method for number of ligand receptors: "noNorm", "bySample", "byCell", "byCellBySample"
process_LR_forSankeyPlot<- function(sigLRs_plot, majorCellType1, majorCellType2, norm_method="noNorm"){
  ## Subset to sender/receiver major cell type
  sigLRs_plot<-sigLRs_plot[sigLRs_plot$sender_major %in% 
                             c(majorCellType1, majorCellType2) & sigLRs_plot$receiver_major %in% c(majorCellType1, majorCellType2) ,]
  
  ## Remove \n from cell types for better plotting
  sigLRs_plot$sender<- gsub("\n", "_", sigLRs_plot$sender)
  sigLRs_plot$receiver<- gsub("\n", "_", sigLRs_plot$receiver)
  
  ## For each unique sender/receiver combination, determine the total number of sig LRs
  sigLRs_plot$interaction_id<- paste0(sigLRs_plot$sender, "__", sigLRs_plot$receiver)
  sigLRs_nLR<- sigLRs_plot %>% group_by(interaction_id) %>% summarise(nLRs=length(unique(LR))) %>% as.data.frame()
  
  
  ## Split sender/receiver back up
  sigLRs_nLR$sender<- sapply(sigLRs_nLR$interaction_id, function(x){unlist(strsplit(x, "__"))[[1]]})
  sigLRs_nLR$receiver<- sapply(sigLRs_nLR$interaction_id, function(x){unlist(strsplit(x, "__"))[[2]]})
  
  ## Normalize based on set method
  ## Divide by number of samples tested in pair
  if(norm_method=="bySample"){
    ## For each unique sender/receiver combination, determine the total number of samples tested, add to number of LRs
    sigLRs_nSamples<-  sigLRs_plot %>% group_by(interaction_id) %>%
      summarise(nSamples_inInteraction=length(unique(Sample))) %>% as.data.frame()
    sigLRs_nLR$nSamples<- as.numeric(plyr::mapvalues(sigLRs_nLR$interaction_id, sigLRs_nSamples$interaction_id, sigLRs_nSamples$nSamples_inInteraction))
    sigLRs_nLR$nLRs_norm<- sigLRs_nLR$nLRs/sigLRs_nLR$nSamples
  }
  
  ## Divide by number of cells tested in pair
  if(norm_method=="byCell"){
    ## Remove sample/cell type with < 10- not input into cellchat
    nCells_bySample_filt<- nCells_bySample[nCells_bySample$NCells>10,]
    
    nCells_total<- nCells_bySample %>% group_by(CellType) %>% summarise(nCells=sum(NCells)) %>% as.data.frame()
    
    ## Cycle through each interaction, determining the number of cells that were input into cellchat
    ## Note that this needs to be done on sample-wise basis
    cells_perPair<- apply(sigLRs_nLR, 1, function(x){
      nCells_pair<- nCells_bySample_filt[nCells_bySample_filt$CellType %in% c(x["sender"],x["receiver"]),]
      nCells_pair_filt<- nCells_pair %>% group_by(Sample) %>% summarise(PF=sum(NCells>10)) %>% as.data.frame()
      samples_used<- as.character(nCells_pair_filt[nCells_pair_filt$PF==2, "Sample"])
      nCells_pair<- nCells_pair[nCells_pair$Sample %in% samples_used,]
      nCells_total<- sum(nCells_pair$NCells)
      return(nCells_total)
    })
    sigLRs_nLR$nCells<- cells_perPair
    sigLRs_nLR$nLRs_norm<- sigLRs_nLR$nLRs/sigLRs_nLR$nCells
  }
  
  ## Divide by number of cells tested in pair FOR EACH SAMPLE, then take mean across all samples
  if(norm_method=="byCellBySample"){
    
    ## For each unique sender/receiver/sample combination, determine the total number of sig LRs
    sigLRs_plot$interaction_id<- paste0(sigLRs_plot$Sample, "__", sigLRs_plot$sender, "__", sigLRs_plot$receiver)
    sigLRs_nLR<- sigLRs_plot %>% group_by(interaction_id) %>% summarise(nLRs=length(unique(LR))) %>% as.data.frame()
    
    ## Split sender/receiver back up
    sigLRs_nLR$sample<- sapply(sigLRs_nLR$interaction_id, function(x){unlist(strsplit(x, "__"))[[1]]})
    sigLRs_nLR$sender<- sapply(sigLRs_nLR$interaction_id, function(x){unlist(strsplit(x, "__"))[[2]]})
    sigLRs_nLR$receiver<- sapply(sigLRs_nLR$interaction_id, function(x){unlist(strsplit(x, "__"))[[3]]})
    
    ## Get number of cells per interaction
    sigLRs_norm<- apply(sigLRs_nLR, 1, function(x){
      nCells_pair<- nCells_bySample[nCells_bySample$CellType %in% c(x["sender"], x["receiver"]) & nCells_bySample$Sample==x["sample"],]
      nCells_pair<- min(nCells_pair$NCells)
      nLR_norm<- as.numeric(x["nLRs"])/nCells_pair
      return(nLR_norm)
    })
    sigLRs_nLR$nLRs_norm<- sigLRs_norm
    
    ## Group by interaction and take mean normalized nLRs_norm value
    sigLRs_nLR$interaction_id<- paste0(sigLRs_nLR$sender, "__", sigLRs_nLR$receiver)
    sigLRs_nLR<- sigLRs_nLR %>% group_by(interaction_id) %>% summarise(nLRs_norm=mean(nLRs_norm)) %>% as.data.frame()
    sigLRs_nLR$sender<- sapply(sigLRs_nLR$interaction_id, function(x){unlist(strsplit(x, "__"))[[1]]})
    sigLRs_nLR$receiver<- sapply(sigLRs_nLR$interaction_id, function(x){unlist(strsplit(x, "__"))[[2]]})
    
  }
  
  ## No norm
  if(norm_method=="noNorm"){
    sigLRs_nLR$nLRs_norm<- sigLRs_nLR$nLRs
  }
  
  ## Get vector of both cell types used
  majorCells1<- unique(sigLRs_plot$sender[sigLRs_plot$sender_major %in% majorCellType1])
  majorCells2<- unique(sigLRs_plot$sender[sigLRs_plot$sender_major %in% majorCellType2])
  
  
  ## Add on major cell type
  sigLRs_nLR$major_sender<- ifelse(sigLRs_nLR$sender %in% majorCells1, majorCellType1, majorCellType2 )
  sigLRs_nLR$major_receiver<- ifelse(sigLRs_nLR$receiver %in%majorCells1, majorCellType1, majorCellType2 )
  sigLRs_nLR$major_sr<- paste0(sigLRs_nLR$major_sender, "\nto\n", sigLRs_nLR$major_receiver)
  
  return(sigLRs_nLR)
}

## Helper function for nice sankey/river plot
## Input: output from process_LRs_forSankeyPlot. Any ordering of programs/cell types should be done prior to my_sankeyPlot
## facet_grid: whether to group rows by major sender
my_sankeyPlot<- function(sigLRs_plot, facet_grid=TRUE){
  p<- ggplot(data=sigLRs_plot, aes(axis1=sender, axis2=receiver,  y=nLRs_norm))+
    geom_alluvium(aes(fill = nLRs_norm)) +
    geom_stratum() +
    theme_classic()+
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Sender", "Receiver"),
                     expand = c(0.15, 0.05), position="top") +
    #scale_fill_gradientn(colors=c("blue", "green", "yellow", "orange", "red"))+
    #scale_fill_viridis_c()+
    paletteer::scale_fill_paletteer_c("viridis::plasma")+
    ylab("")+
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=12))
  if(facet_grid){
    p+theme(strip.text = element_blank())+
      facet_grid(rows=vars(major_sender), scale="free_y", space="free_y")
  }else{p}

}
