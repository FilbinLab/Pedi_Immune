#!/usr/bin/env Rscript

## Libraries
library(CellChat)
library(Seurat) 

## Helper function
cellchatHelper<-paste0("~/scripts/R/CellChat_HelperFunctions.R")
source(cellchatHelper)


## args: seurat object name, column in metadata to use, output directory, cell type of interest #1, cell type of interest #2, whether or not to downsample (TRUE/FALSE)
## cell types #1 and #2 should be present in under annotColumn within the seurat object's metadata

args=commandArgs(trailingOnly=TRUE)

seurat_obj<- args[1]
annotColumn<- args[2]
output_dir<- args[3]
celltype_1<- args[4]
celltype_2<- args[5]
downsample<- args[6]


message(paste0("Processing comparison: ", celltype_1, " vs ", celltype_2))

## Load in list of seurat objects
## 1 object for every sample, with tumor + Tcell data combined
seurat_list<- readRDS(seurat_obj)

## Set column to use as cell annotations
## Can use "CellAnnot" for merged annotations (like OPC-like*) or "OriginalAnnot" for completely un-merged annots
seurat_list<- lapply(seurat_list, function(x){
  tmp<-x
  tmp$CellAnnot<- tmp@meta.data[,annotColumn]
  return(tmp)
})



## CellChat settings
Pvaluethresh<- 1
PCthresh<-0.2
raw.use=TRUE
population.size = FALSE
min.cells = 10
only.pos=TRUE



## remove any samples that don't have at least 10 cells for each cell type of interest
nCellsOfInterest<- lapply(seurat_list, function(x){
  tbl<- table(x$CellAnnot)
  c(tbl[celltype_1], tbl[celltype_2])
})

samplePF<- lapply(nCellsOfInterest, function(x){
  x[celltype_1]>10 & x[celltype_2] >10})
seurat_list<- seurat_list[unlist(samplePF) & !(is.na(samplePF))]

all_samples<-names(seurat_list)
message(paste0("Samples used: ", paste(all_samples, collapse=", ")))


## If downsample- randomly sample celltype 1 and celltype 2 so they have the same # of cells as the smallest program in their respective broad subtype (min 11)
if(downsample){
  seurat_list<- lapply(seurat_list, function(seurat){
    ## Set broad cell types used for celltype 1 and celltype 2
    celltype_1_broad<-as.character(plyr::mapvalues(celltype_1, seurat$CellAnnot, seurat$highLevel_annot, warn_missing = FALSE))
    celltype_2_broad<-as.character(plyr::mapvalues(celltype_2, seurat$CellAnnot, seurat$highLevel_annot, warn_missing = FALSE))
      
    ## Get number of cells per program, with high level
    nCells<- as.data.frame(table(seurat$CellAnnot))
    colnames(nCells)<-c("Program", "NCells")
    nCells$broad<- as.character(plyr::mapvalues(nCells$Program, seurat$CellAnnot, seurat$highLevel_annot, warn_missing = FALSE))
    
    ## Remove any programs without >10 cells- not running in CellChat
    nCells<- nCells[nCells$NCells>10,]
    
    ## Get min number of cells in a program, by broad group
    minCells<- sapply(unique(nCells$broad), function(x){
      min(nCells[nCells$broad==x,"NCells"])
    })
    
    ## Randomly sample celltype1 down to min for its broad program
    ds_celltype1<- sample(x=colnames(seurat)[seurat$CellAnnot==celltype_1],
                          size=minCells[celltype_1_broad], 
                          replace=FALSE)
    
    ## Randomly sample celltype2 down to min for its broad program
    ds_celltype2<- sample(x=colnames(seurat)[seurat$CellAnnot==celltype_2],
                          size=minCells[celltype_2_broad], 
                          replace=FALSE)
    
    ## Subset seurat to just these cells
    seurat<- subset(seurat, cells=c(ds_celltype1, ds_celltype2))
    return(seurat)
  })
}else{seurat_list<- lapply(seurat_list, function(x){ ## If not downsampling- just use all cells in celltype1 and celltype2
  subset(x, CellAnnot %in% c(celltype_1, celltype_2))
})
}


## Create object, preprocess
cellChat<- lapply(all_samples, function(x){
  CreateCellChatObject(mySample = x, 
                       all_merged_seurat_obj = seurat_list,
                       SpecifyCellAnnot = FALSE, ## already subset to cells of interest
                       Pvaluethresh = Pvaluethresh,
                       PCthresh =PCthresh,
                       only.pos=only.pos) 
})
names(cellChat)<- all_samples

message(paste0("Computing interaction network for ", celltype_1," vs ", celltype_2))


## skip if no cellchat objects for that subtype
if(class(cellChat) != "list"){
  message("Not enough cells present in either group- cellchat not run!")
  quit
}else{
  ## Compute interaction network
  cellChat_interact<- lapply(cellChat, function(x){
    InferCellCellInteractions(x, 
                              raw.use=raw.use, 
                              population.size = population.size, 
                              min.cells = min.cells)
  })
  
  
  saveRDS(cellChat_interact, file=paste0(output_dir, "/cellChat_res_", celltype_1,".vs.", celltype_2, ".Rds"))    
  
  ## Extract results
  PathwaysRemove<-c("HLA", "COL", "LAMININ", "MHC-I")
  sigLRs_list<- lapply(names(cellChat_interact), function(sample){
    cellchat<- cellChat_interact[[sample]]
    if(length(cellchat@netP$pathways)>1){
      sigLRs<-  subsetCommunication(cellchat, thresh = 1)
      sigLRs<- sigLRs[!sigLRs$pathway_name %in% PathwaysRemove,]
      sigLRs$sample<- sample
      return(sigLRs)
    }else{print(paste0("No LRs for ", sample)); return(NULL)}
  })
  names(sigLRs_list)<-names(cellChat_interact)
  sigLRs_list<- sigLRs_list[unlist(lapply(sigLRs_list, function(x){!is.null(x)}))]
  
  ## collapse to DF
  sigLRs_df<- do.call("rbind", sigLRs_list)
  
  
  saveRDS(sigLRs_df, file=paste0(output_dir, "/cellChat_LRs_",celltype_1, ".vs.", celltype_2, ".Rds"))
}
