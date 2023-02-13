# Load libraries

## Libraries for basic preprocessing 
library(reshape2)
library(dplyr)

## Single cell libraries
library(Seurat)
#library(pagoda2)
##library(conos)
library(harmony)

## Libraries for plotting 
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)

## Annotation and pathway analysis library 
##library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

## Others
library(Rtsne)

# Define functions

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


## Functions for preprocessing count matrix 







## QC
## @param cm count matrix (raw or log transformed)
## @param hg_list list of house-keeping genes 
## @param scale_factor factor to scale tpm, default = 10
## @param log_base base of log, default = 2
## @param gene_detect detection limit for gene, default = 1
calculateQcMetrics <- function(cm, hg_list, scale_factor=10, log_base=2, gene_detect=1){
  ## Info for sample, plate and well
  sample = sapply(colnames(cm), function(x) unlist(strsplit(x, split=".", fixed=T))[1])
  plate = sapply(colnames(cm), function(x) unlist(strsplit(x, split=".", fixed=T))[2])
  well = sapply(colnames(cm), function(x) unlist(strsplit(x, split=".", fixed=T))[3])
  ## Total number of genes 
  gene = colSums(cm >= gene_detect)
  cm_norm = log2(cm/scale_factor+1)
  ## Average expression of house keeping genes
  hg_list = hg_list[hg_list %in% rownames(cm)]
  hk = colMeans(cm_norm[hg_list,])
  qc_df = cbind.data.frame(sample, plate, well, gene, hk)
  return(qc_df)
}






## Compute plain average expression or control corrected signature score 
## @param X.center centered relative expression
## @param X.mean average of relative expression of each gene (log2 transformed)
## @param n number of genes with closest average expression for control genesets, default = 100 
## @param simple whether use average, default  = FALSE
scoreSignature <- function(X.center, X.mean, s, n=100, simple = FALSE, verbose=FALSE) {
  if(verbose) {
    message("cells: ", ncol(X.center))
    message("genes: ", nrow(X.center))
    message("genes in signature: ", length(s))
    message("Using simple average?", simple)
    message("processing...")
  }
  
  s <- intersect(rownames(X.center), s)
  message("genes in signature, and also in this dataset: ", length(s))
  ##message("These genes are: ", s)
  
  if (simple){
    s.score <- colMeans(X.center[s,])
  }else{
    s.score <- colMeans(do.call(rbind, lapply(s, function(g) {
      # g <- s[2]
      # message(g)
      if(verbose) message(".", appendLF = FALSE)
      g.n <- names(sort(abs(X.mean[g] - X.mean))[2:(n+1)])
      X.center[g, ] - colMeans(X.center[g.n, ])
    })))
  }
  
  if(verbose) message(" done")
  return(s.score)
}

## Convert gene symbols to ensembl ID                
gene_symbol_to_ensembl_id <- function(gene, dataset){
  ensembl = useMart("ensembl")
  mart <- useDataset(dataset, useMart("ensembl"))
  ids <- getBM(filters="external_gene_name", attributes=c("ensembl_gene_id", "external_gene_name"), values=gene, mart=mart)
  non_duplicates <- which(duplicated(ids$external_gene_name) == FALSE)
  ids <- ids[non_duplicates, ] 
}

## Run GO term over-representation analysis                
go_analysis <- function(sigOE_genes, allOE_genes, keyType="ENSEMBL", 
                        OrgDb="org.Hs.eg.db", ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05){
  ego <- enrichGO(gene = sigOE_genes, 
                  universe = allOE_genes, 
                  keyType = keyType, 
                  OrgDb = OrgDb, 
                  ont = ont, 
                  pAdjustMethod = pAdjustMethod, 
                  qvalueCutoff = qvalueCutoff, 
                  readable = TRUE)
  
  cluster_summary <- subset(data.frame(ego), select = -c(geneID))
  cluster_summary <- data.frame(ego)
  return(list(ego=ego, cluster_summary=cluster_summary))
}

## Create seurat obj, normalize, find variable genes, and scale data 
## @param cm, count matrix
## @param project project name
## @param min.cells minimal number of cells required to express a gene
## @param min.genes minimal number of genes required to be expressed in a cell
## @param is.expr detection limit for a gene
## @param scale.factor scale factor normalization, default=1E5
## @param do.scale whether to scale the data (divided by sd), default=F
## @param do.center whether to center the data (substract mean), default=T
preprocessSeuratObject <-function(cm, project, min.cells=0, min.genes=0, 
                                  scale.factor=1E5, do.scale=F, do.center=T){
  ## Create Seurat Obj
  seurat_obj = CreateSeuratObject(cm, min.cells=min.cells, min.features = min.genes,  project=project)
  
  ## Normalize data
  seurat_obj <- NormalizeData(
    object = seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = scale.factor
  )
  
  ## Detection of variable genes across the single cells
  seurat_obj <- FindVariableFeatures(
    object = seurat_obj,
    mean.function = ExpMean, 
    dispersion.function = LogVMR 
    #x.low.cutoff = 0.0125, ## Default is 0.0125 
    #x.high.cutoff = 8, ## Default is 8
    #y.cutoff = 1 ## Default is 1 
  )
  ##length(x = seurat_obj@var.genes)
  
  ## Scaling the data and removing unwanted sources of variation
  seurat_obj <- ScaleData(
    object = seurat_obj,
    do.scale = do.scale,
    do.center = do.center
    #vars.to.regress = c("nUMI")
  )
  
  return(seurat_obj)
}

## Normalize and center count matrix, for input into ScoreSignature
NormCenter<-function(cm, scale_factor=10, log_base=2){
  cm = as.matrix(cm)
  cm_norm = log(cm/scale_factor+1, base=log_base)
  cm_norm_center = cm_norm-rowMeans(cm_norm)
  result = list()
  result[["raw_data"]] = cm
  result[["norm_data"]] = cm_norm
  result[["center_data"]] = cm_norm_center
  return (result)
}

## Convert v3 seurat object to v4 seurat object
## from raw data --> UMAP
ConvertV3ToV4<- function(seurat_v3, project, min.cells=0, min.genes=0, 
                         scale.factor=1E5, do.scale=F, do.center=T,
                         harmonyTheta=2, resolution=0.8, dims=20){
  ## remove any loaded seurat package, reload v4 seurat
  detach("package:Seurat", unload=TRUE)
  library(Seurat, lib.loc = "C:/Users/jenna/OneDrive/Documents/R/win-library/4.0/Seurat/V4")
  
  ## extract counts and meta from v3 seurat object
  cm<- seurat_v3@assays$RNA@counts
  cm<-as.matrix(cm)
  meta<-seurat_v3@meta.data
  meta<- meta[,!(colnames(meta) %in% c("orig.ident", "nCount_RNA", "nFeature_RNA"))]
  
  ## Create new seurat object
  seurat_v4<- CreateSeuratObject(counts=cm, min.features = 0, min.cells = 0,
                                 meta.data = meta)
  
  ## Normalize data
  seurat_v4 <- NormalizeData(
    object = seurat_v4,
    normalization.method = "LogNormalize",
    scale.factor = scale.factor
  )
  
  ## Detection of variable genes across the single cells
  seurat_v4 <- FindVariableFeatures(
    object = seurat_v4,
    mean.function = ExpMean, 
    dispersion.function = LogVMR 
  )
  ##length(x = seurat_obj@var.genes)
  
  ## Scaling the data and removing unwanted sources of variation
  seurat_v4 <- ScaleData(
    object = seurat_v4,
    do.scale = do.scale,
    do.center = do.center
    #vars.to.regress = c("nUMI")
  )
  
  ## PCA
  seurat_v4 <- RunPCA(
    object = seurat_v4, 
    features = seurat_v4@assays$RNA@var.features, 
    npcs = 100,
    verbose = TRUE, 
    ndims.print  = 1:5, 
    nfeatures.print = 5
  )
  
  ## Harmony
  seurat_v4 = RunHarmony(seurat_v4, "sample", theta = harmonyTheta, 
  max.iter.harmony = 50, plot_convergence = TRUE)
  
  ## Run UMAP
  seurat_v4 <- RunUMAP(seurat_v4, reduction = "harmony", dims = 1:dims)
  
  
  seurat_v4 <- FindNeighbors(seurat_v4,
                              reduction = "harmony",
                              dims = 1:dims,
                              force.recalc = TRUE)   %>% 
    FindClusters(resolution = resolution)
  
  detach("package:Seurat", unload=TRUE)
  return(seurat_v4)
}

## Full seurat v3 analysis
RunFullSeurat<- function(cm, do.scale=FALSE, do.center=TRUE, RunHarmony,samples,
                         resolution = .8, dims = 20, pca_dims=100, n.neighbors=30, theta=2){
  seurat_obj = preprocessSeuratObject(cm, project="", min.cells=0, min.genes=0, 
                                      scale.factor=1E5, do.scale=do.scale, do.center=do.center)
  seurat_obj@meta.data$sample = samples
  
  seurat_obj <- RunPCA(
    object = seurat_obj, 
    features = seurat_obj@assays$RNA@var.features, 
    npcs = pca_dims,
    verbose = TRUE, 
    ndims.print  = 1:5, 
    nfeatures.print = 5
  )
  
  if(RunHarmony){
    seurat_obj = RunHarmony(seurat_obj, "sample", theta = theta, 
                                         max.iter.harmony = 50, plot_convergence = TRUE)
    seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:dims, n.neighbors = n.neighbors)
    seurat_obj <- FindNeighbors(seurat_obj,
                              reduction = "harmony",
                              dims = 1:dims,
                              force.recalc = TRUE)   %>% 
    FindClusters(resolution = resolution)} else{
      ## Run UMAP
      seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:dims, n.neighbors = n.neighbors)
      seurat_obj <- FindNeighbors(seurat_obj,
                                  reduction = "pca",
                                  dims = 1:dims,
                                  force.recalc = TRUE)   %>% 
        FindClusters(resolution = resolution)
    }
  
  return(seurat_obj)
}


## Function for running seurat pipeline
## For subsetting --> rerunning
## Same as Orr's pipeline shared
RunFullSeurat_Immune<- function(cm, samples){
  gcdata<-CreateSeuratObject(counts = cm, min.cells = 0, min.features = 0, project = "")
  gcdata$sampleid<- samples
  nSamples<-table(gcdata$sampleid)
  gcdata<- subset(gcdata, sampleid %in% names(nSamples[nSamples>1]))
  gcdata <- NormalizeData(gcdata, normalization.method = "LogNormalize", scale.factor = 10000)
  gcdata_list<- SplitObject(gcdata, split.by = "sampleid")
  var.genes <- SelectIntegrationFeatures(gcdata_list, 
                                         nfeatures = 2000, verbose = TRUE, fvf.nfeatures = 2000, 
                                         selection.method = "vst")
  VariableFeatures(gcdata) <- var.genes
  gcdata <- ScaleData(gcdata, split.by = "sampleid", features = VariableFeatures(gcdata), 
                      do.center = T, do.scale = F)
  gcdata <- RunPCA(gcdata, features = VariableFeatures(gcdata), npcs = 40)
  gcdata <- FindNeighbors(gcdata, reduction = "pca", dims = 1:20, k.param = 20)
  gcdata <- FindClusters(gcdata, resolution = 1, algorithm = 4, random.seed = 100)  # conda activate /Users/jlabelle/Library/r-miniconda/envs/r-reticulate
  gcdata <- RunUMAP(gcdata, dims = 1:20, reduction = "pca", n.neighbors = 15, min.dist = 0.5, spread = 1, metric = "euclidean", seed.use = 1)  
  return(gcdata)
}



## Plot GO results (from runGO)
## Input: list of go_results, 1 for each cell type. Output from runGO.
## Input: n_terms: number of GO terms to plot
plotGO<- function(go_result, n_terms=20){
  library(stringr)
  
  geneSets<-names(go_result)
  allDotPlots<-list()
  for (i in geneSets){
    print(i)
    if (nrow(go_result[[i]]$cluster_summary)>0){
      print("Sig terms present")
      allDotPlots[[i]]<-dotplot(go_result[[i]]$ego, 
                                showCategory=n_terms, 
                                font.size = 15, 
                                title = i, 
                                label_format=10) +
        scale_y_discrete(labels=function(x) str_wrap(x, width=60))+
        theme(plot.title = element_text(hjust = 0.5, face = "bold", 
                                        color="black", size = 28),
              axis.title = element_text(face = "bold", color="black"), 
              axis.text.x = element_text(angle = 45, hjust = 1, color="black"),
              axis.text.y = element_text( color="black",face = "bold"))
    }
  }
  
  allDotPlots<- allDotPlots[!(unlist(lapply(allDotPlots, function(x){is.null(x)})))]
  return(allDotPlots)
}

## function for pseudobulking by variable
pseudobulk_byVariable<- function(cm, meta, variableToGroupBy){
  ## remove groups with only 1 cell
  tmp<- table(meta[,variableToGroupBy])
  group_names<- names(tmp[tmp!=1])
  
  pseudobulk = NULL
  if(any(colnames(cm) != rownames(meta))){print("cm / meta cell names do not match"); break}
  
  for (group in group_names){
    cm_tmp = cm[,meta[,variableToGroupBy] == group]
    pseudobulk = cbind(pseudobulk, rowMeans(cm_tmp))
  }
  pseudobulk<- as.data.frame(pseudobulk)
  colnames(pseudobulk) = group_names
  rownames(pseudobulk) = rownames(cm)
  return(pseudobulk)
}


# Function for heatmap plot, for set of GOI, Using ggplot
## Input count matrix should be pseudobulked counts, or grouped/averaged in some other way
## Column names should be 2 grouping variables, separated by "_". The first half of the column names will be the x axis (subtype usually)
## nameForCounts is just used for labeling legend appropriately
## add_lines = TRUE adds dashed lines between genesets. Requires "top_markers" which is df of all markers used, with 'cluster' column for each group
## addBoxes = for adding boxes (usually for sig) to heatmap cells. Should be df with 3-4 columns: XAxis=xaxis value, Gene=gene value, Boxes=whether to add box to that cell, if using facet wrapping, FacetWrap=variable used for facet wrapping
myHeatmap<- function(pseudobulked_cm, max.value=1e6, min.value=-1e6, nameForCounts="TPM", orderSubtypes="None", 
                     GOI, facetWrap=TRUE, orderFactors="None", add_lines=FALSE, top_markers=top_markers, addBoxes=FALSE){
  goi_cm<- as.data.frame(t(pseudobulked_cm[GOI,]))
  goi_cm$UniqueID<- rownames(goi_cm)
  goi_cm<- melt(goi_cm); colnames(goi_cm)<- c("UniqueID", "gene", "counts")
  goi_cm$gene<- factor(goi_cm$gene, levels=rev(unique(goi_cm$gene)))
  
  ## Split UniqueID into X axis (first half) and faceting variable (second half)
  tmp<- strsplit(goi_cm$UniqueID, split="_")
  goi_cm$XAxis<- unlist(lapply(tmp, function(x){x[1]}))
  goi_cm$FacetWrap<- unlist(lapply(tmp, function(x){x[2]}))
  
  ## Enforce min/max values
  goi_cm$counts<- ifelse(goi_cm$counts>max.value, max.value, ifelse(goi_cm$counts<min.value, min.value, goi_cm$counts))
  
  ## Order subtypes, if desired
  if(orderSubtypes[1] != "None"){
    goi_cm$XAxis<- factor(goi_cm$XAxis, levels=orderSubtypes)
  }
  
  ## Order factors, if desired
  if(orderFactors[1] != "None"){
    goi_cm$FacetWrap<- factor(goi_cm$FacetWrap, levels=orderFactors)
  }
  
  ## plot
  p<- ggplot(goi_cm, aes(x=XAxis, y= gene, fill=counts))+
    geom_tile()+
    theme_classic()+
    scale_fill_gradient2(high="red", low="yellow", mid="white", name=nameForCounts)+
    scale_x_discrete(position = "top")+
    theme(axis.text.x = element_text(angle=45, hjust=0, color="black", face="bold", size=16),
          axis.text.y = element_text(color="black", face="bold"))+
    xlab("")+ylab("")
  if(facetWrap){
    p<- p + facet_grid(cols=vars(FacetWrap), scale="free_x", space="free_x", switch="x")+
      theme(strip.text.x = element_text(size=16, color="black", face="bold"))
  }
  
  ## Add lines if needed
  if(add_lines){
    tmp<- split(top_markers, f=top_markers$cluster)
    lastGene<- lapply(tmp, function(x){
      x$tmp<- 1:nrow(x)
      x<- x[x$gene %in% rownames(pb), ]
      return(x[max(x$tmp), "gene"])
    })
    linePositions<- as.data.frame(pb_center[top_markers$gene,])
    linePositions$tmp<- 1:nrow(linePositions)
    linePositions<- linePositions[unlist(lastGene),"tmp"]
    linePositions<- linePositions[1:length(linePositions)-1]+ 0.5
    
    p<- p+geom_hline(yintercept = linePositions, linetype="dashed",size=1)
  }else{p<-p}
  
  ## add boxes (usually denoting significance) if needed
  if(addBoxes != FALSE){
    goi_cm$UniqueID<- paste0(goi_cm$XAxis, "_", goi_cm$FacetWrap, "_", goi_cm$gene)
    addBoxes$UniqueID<- paste0(addBoxes$XAxis, "_", addBoxes$FacetWrap, "_", addBoxes$Gene)
    goi_cm$boxes<- plyr::mapvalues(goi_cm$UniqueID, addBoxes$UniqueID, addBoxes$Box, warn_missing = FALSE)
    p<-p+geom_tile(data = goi_cm[goi_cm$boxes==TRUE,], fill = NA, color = "black", size = .5)
  }else{p<-p}
  return(p)
}

## Helper for barchart by program proportion
## grouping variable is x axis, coloring variable is fill, wrapping variable (optional) is for facet wrap (FALSE or name of column)
## grouping_variable_order = optional sample ordering. Number_or_proportion = plot by "proportion" of cells or raw "number" of cells
my_barchart_programProp<- function(seurat, grouping_variable="sample", coloring_variable, wrapping_variable=FALSE,
                                   colors, grouping_variable_order=NULL, number_or_proportion="proportion"){
  df<- as.data.frame(table(seurat@meta.data[[grouping_variable]], seurat@meta.data[[coloring_variable]]))
  colnames(df)<- c("Sample", "Program", "NCells")
  df<- df %>% group_by(Sample) %>% mutate(perCells=NCells/sum(NCells)) %>% as.data.frame()
  
  if(!is.null(grouping_variable_order)){
    df$Sample<- factor(df$Sample, levels=grouping_variable_order)
  }
  
  if(wrapping_variable != FALSE){
    df$wrap<- plyr::mapvalues(df$Sample,seurat@meta.data[[grouping_variable]], seurat@meta.data[[wrapping_variable]], warn_missing = FALSE)
  }
  
  if(number_or_proportion=="proportion"){
    df$cells<- df$perCells
  }else if(number_or_proportion=="number"){
    df$cells<- df$NCells
  }
  p<- ggplot(df, aes(x=Sample, y=cells, fill=Program))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=colors[names(colors) %in% df$Program])+
    theme_classic()+
    theme(axis.text.x = element_text(angle=45, hjust=1, face="bold", color="black"),
          axis.text.y = element_text(face="bold", color="black"),
          axis.title =element_text(face="bold", color="black") )+
    ylab(paste0(number_or_proportion, " of cells"))
  if(wrapping_variable !=FALSE){
    p+facet_grid(cols=vars(wrap), scales = "free_x", space="free_x")
  }else{p}
}
