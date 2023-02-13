library(Seurat)
library(ggplot2)
library(ggpubr)
library(Rtreemix)




## Function: for a given sample, subset clonal/nonclonal seurat into that sample, normalize
Exp_SubsetAndNormalize<- function(seurat_list, mySample){
  ## get clonal/nonclonal seurat for each sample
  clonal_tmp<- subset(seurat_list$clonal, sampleid==mySample)
  nonclonal_tmp<- subset(seurat_list$nonclonal, sampleid==mySample)
  
  ## Normalize- logTPM / logTP100k / TPM / TP100K
  clonal_norm<- log(clonal_tmp@assays$RNA@counts/10 + 1)
  nonclonal_norm<- log(nonclonal_tmp@assays$RNA@counts/10 + 1)
  
  norm_list<- list(clonal=clonal_norm, nonclonal=nonclonal_norm)
  subset_seurat_list<- list(clonal=clonal_tmp, nonclonal=nonclonal_tmp)
  return(list(norm=norm_list, seurat=subset_seurat_list))
}

## Function: get clonal cells in clonotype with > size 1
## For Exp: Return a table of # of cells per clonotype 
## For CV: Return a list of clonal cells, separated by clonotype
Exp_GetClonotypeCells<- function(subset_seurat_list, min_clonotype_size=1, mySample){
  clonal_tmp<- subset_seurat_list$clonal
  
  ## Get clonal cells
  clonotype_table<- table(clonal_tmp$clonal_group)
  clonotype_table<- clonotype_table[clonotype_table>min_clonotype_size] ## skipping clonotypes of less than min size
  clonotype_table<- clonotype_table[names(clonotype_table)!="NA"] ## remove clonotypes of "NA"
  
  ## if no clonal cells- return NULL
  if(length(names(clonotype_table))==0){print(paste0("No clonotypes in ", mySample)); return(NULL)} 
  
  ## Subset clonal seurat to cells that are in clonotype >1
  clonal_tmp<- subset(clonal_tmp, clonal_group %in% names(clonotype_table))
  
  ## Create simple flat vector of all clonal cells
  flat_clonal_cells<- unique(colnames(clonal_tmp))
  
  ## List of all clonal cells, grouped by clonotype- used for CV
  clonal_cells_tmp<- split(clonal_tmp@meta.data, f=clonal_tmp$clonal_group) 
  clonal_cells_list<- lapply(clonal_cells_tmp, function(x){rownames(x)}) 
  
  clonal_cells<- list(clonotype_table=clonotype_table, 
                      clonal_cells_list=clonal_cells_list,
                      flat_clonal_cells=flat_clonal_cells)
  return(clonal_cells)
  
}

## Function: randomly select non-clonal cells based on a table of clonal cells, repeat 10X for each clonotype
## For exp: return simple list of unique nonclonal cells selected
## For CV: return list of nonclonal cells selected, separated by the clonotype/iteration they were chosen for
Exp_SelectNonClonal<- function(clonotype_table, subset_seurat_list){
  set.seed(42)
  nonclonal_tmp<- subset_seurat_list$seurat$nonclonal
  ## For 10 iterations, select nonclonal cells for each clonotype
  all_nonclonal_cells_set<-lapply(1:10, function(y){
    nonclonal_cells<- lapply(clonotype_table, function(x){sample(colnames(nonclonal_tmp),x)})
  })
  
  ## Create flat vector of nonclonal cells for expression
  all_nonclonal_cells<- unname(unlist(unname(all_nonclonal_cells_set))) 
  
  ## Create list of nonclonal cells selected, separated by clonotype- for CV
  nonclonal_cells_vec<- c(); for(i in all_nonclonal_cells_set){nonclonal_cells_vec<- c(nonclonal_cells_vec, i)}
  
  nonclonal_cells<- list(flat_nonclonal_cells=all_nonclonal_cells,
                         list_nonclonal_cells=nonclonal_cells_vec)
  return(nonclonal_cells)
}

## Wrapper function that cycles through each sample and then:
## 1) normalizes
## 2) selects clonal cells (clonotype>1) + random nonclonal selection
## 3) gets clonal/nonclonal df based on cells in #2
## For expression: merges clonal into single df for all cells, same for nonclonal
## For CV: vector of all clonal cells separated by clonotype, same for nonclonal
Exp_WrapperForMultipleSamples<- function(seurat_list, min_clonotype_size){
  ## For expression
  clonal_counts<- list()
  nonclonal_counts<- list()
  
  ## For CV
  nonclonal_cells_list<- c() 
  clonal_cells_list<- c() 
  
  
  for(mySample in unique(seurat_list$clonal$sampleid)){
    print(mySample)
    
    subset_list<- Exp_SubsetAndNormalize(seurat_list, mySample)
    
    ## number of cells present in each clonotype
    clonal_cells<- Exp_GetClonotypeCells(subset_list$seurat, min_clonotype_size = min_clonotype_size, mySample = mySample)
    
    ## If there are any clonal cells: 
    if(class(clonal_cells)=="list"){
      ## randomly selected nonclonal cells
      nonclonal_cells<- Exp_SelectNonClonal(clonal_cells$clonotype_table, subset_list)
      
      ## Get count matrices based on clonal/nonclonal cells
      nonclonal_cm<- subset_list$norm$nonclonal[,unique(nonclonal_cells$flat_nonclonal_cells)]
      clonal_cm<- subset_list$norm$clonal[,clonal_cells$flat_clonal_cells]
      
      ## For expression
      nonclonal_counts[[mySample]]<- nonclonal_cm
      clonal_counts[[mySample]]<- clonal_cm
      
      
      ## For CV
      nonclonal_cells_list<- c(nonclonal_cells_list, nonclonal_cells$list_nonclonal_cells)
      clonal_cells_list<-c(clonal_cells_list, clonal_cells$clonal_cells_list)
    }else{next}
  }
  ## Bind counds into single df
  clonal_df<- do.call("cbind", clonal_counts)
  nonclonal_df<- do.call("cbind", nonclonal_counts)
  
  expr_list<- list(counts=list(clonal=clonal_df, nonclonal=nonclonal_df),
                   cell_list=list(clonal=clonal_cells_list, nonclonal=nonclonal_cells_list))
  return(expr_list)
  
}


## Function for processing clonal/nonclonal count matrix, calculating zscore
Exp_ProcessCounts<- function(expr_list_counts){
  clonal_df<- expr_list_counts$clonal
  nonclonal_df<- expr_list_counts$nonclonal
  
  ## Calculate mean expression across all cells
  clonal_mean<- as.data.frame(rowMeans(clonal_df)); colnames(clonal_mean)<- "clonal"
  nonclonal_mean<- as.data.frame(rowMeans(nonclonal_df)); colnames(nonclonal_mean)<- "nonclonal"
  
  ## Merge clonal/nonclonal counts together
  CNC_mean<- merge(clonal_mean, nonclonal_mean, by=0)
  rownames(CNC_mean)<- CNC_mean$Row.names; CNC_mean<- CNC_mean[,-1]
  
  ## Caluclate z score
  CNC_mean$ClonalMinusNonClonal<-  CNC_mean$clonal -CNC_mean$nonclonal
  #CNC_mean$zscore<- scale(CNC_mean$ClonalMinusNonClonal, center=TRUE, scale=FALSE)
  CNC_mean$zscore<- (CNC_mean$clonal - CNC_mean$nonclonal) / sd(CNC_mean$nonclonal)
  
  return(CNC_mean)
}

## Function to merge cells from the same clonotype but different samples
CV_MergeClonotypes<- function(expr_list){
  clonal_cells_list<- expr_list$cell_list$clonal
  
  clonal_cells_list<- sapply(unique(names(clonal_cells_list)), function(x) {
    unname(unlist(clonal_cells_list[names(clonal_cells_list)==x]))}, 
    simplify=FALSE)
  
  clonal_cells_list<- clonal_cells_list[order(names(clonal_cells_list))]
  cv_list<- expr_list
  cv_list$cell_list$clonal<- clonal_cells_list
  return(cv_list)
}

## Function for calculating l1 distance between all cell combinations across all clonotypes
## input count matrix/cell list should both be for either clonal or nonclonal cells
## In this function "clonotype" refers to both actual clonal clonotypes as well as their nonclonal counterparts
## cv_list object input just so that the results can be added directly to it
CV_CalculateL1Dist<- function(cv_list, clonal_or_nonclonal, seurat_list, verbose=TRUE, tpm_or_norm="tpm"){
  ## select count matrix/cell list for clonal or nonclonal subset- not sure if norm counts or tpm or appropriate
  if(tpm_or_norm=="tpm"){cm<- as.data.frame(seurat_list[[clonal_or_nonclonal]]@assays$RNA@counts)}
  if(tpm_or_norm=="norm"){cm<- as.data.frame(cv_list$counts[[clonal_or_nonclonal]])}
  
  cells_list<- cv_list$cell_list[[clonal_or_nonclonal]]
  
  ## initialize dataframes for adding distances and cell pairs used
  all_l1Dist_allClonotypes<- data.frame(row.names = rownames(cm), 
                                        dummy=rep("NA", length(rownames(cm))))
  all_cellPairs<- data.frame()
  
  ## Cycle through each clonotype, create list of all unique cell combinations, and calculate distance between all genes for all combinations
  for (clonotype in names(cells_list)){
    if(verbose){print(paste0("Processing clonotype ", clonotype))}
    
    ## Get counts for all cells in clonotype
    clonotype_cells<- cells_list[[clonotype]]
    if(length(clonotype_cells)==1){if(verbose){print("only 1 cell in clonotype; cannot calculate CV")}; next}
    cm_tmp<- cm[,unlist(clonotype_cells)]
    
    ## Create df with all possible pairwise combinations of cells within clonotype
    all_pairs<-as.data.frame(combn(unlist(clonotype_cells),2))
    
    ## Cycle through each pair, calculating l1 distance for all genes
    l1_dist_allCombinations<- apply(all_pairs, 2, function(x){
      cm_tmp2<- cm_tmp[,x]
      l1_dist<- apply(cm_tmp2, 1, function(y){L1.dist(y[1], y[2])})
      return(l1_dist)
    })
    colnames(l1_dist_allCombinations)<- paste0(clonotype, "_", 1:ncol(l1_dist_allCombinations))
    
    ## Merge with all other clonotypes
    all_l1Dist_allClonotypes<- merge(all_l1Dist_allClonotypes, l1_dist_allCombinations, by=0)
    rownames(all_l1Dist_allClonotypes)<- all_l1Dist_allClonotypes$Row.names; all_l1Dist_allClonotypes<- all_l1Dist_allClonotypes[,-1]
    
    ## make df with all cell combinations used to double check
    pairs_df<- data.frame(clonotype_number=clonotype, cell_pairs= paste0(all_pairs[1,], "_vs_", all_pairs[2,]) )
    all_cellPairs<- rbind(all_cellPairs, pairs_df)
  }
  
  ## Remove dummy column from results
  all_l1Dist_allClonotypes<- all_l1Dist_allClonotypes[,colnames(all_l1Dist_allClonotypes) != "dummy"]
  
  cv_list$CVres_l1Dist[[clonal_or_nonclonal]]<- all_l1Dist_allClonotypes
  cv_list$CVres_cellPairs[[clonal_or_nonclonal]]<- all_cellPairs
  
  return(cv_list)
}

## Function to process L1 distances and merge with counts results
CV_ProcessCV<- function(cv_list_CVres_l1Dist, clonal_or_nonclonal){
  all_l1Dist_allClonotypes<- cv_list_CVres_l1Dist[[clonal_or_nonclonal]]
  
  ## Average distance for all genes
  mean_l1Dist_allClonotypes<- as.data.frame(rowMeans(all_l1Dist_allClonotypes))
  colnames(mean_l1Dist_allClonotypes)<- c(paste0("MeanL1Distance_", clonal_or_nonclonal))
  
  return(mean_l1Dist_allClonotypes)
}

## For nonclonal, need to add Set# (1:10) to cell name
CV_AddSetToNonclonalCells<- function(cv_list){
  nonclonal_cells_list<- cv_list$cell_list$nonclonal
  nonclonal_cells_list<- nonclonal_cells_list[order(names(nonclonal_cells_list))]
  names(nonclonal_cells_list)<- paste0(names(nonclonal_cells_list), "_Set", 1:10)
  cv_list$cell_list$nonclonal<- nonclonal_cells_list
  return(cv_list)
}

## Function to merge and process counts and L1 distances:
## 1) Merge all results together: counts for CNC, clonal L1, nonclonal L1
## 2) Then calculate CV for for clonal and nonclonal L1
## 3) Then calculate Zscore for CV
MergeAndProcess_Counts.L1Dist<- function(CV_clonal_res, CV_nonclonal_res, expr_CNC_mean){
  ## Merge clonal L1 Dist with CNC counts
  df_final<- merge(expr_CNC_mean, CV_clonal_res, by=0)
  rownames(df_final)<- df_final$Row.names; df_final<- df_final[,-1]
  
  ## Merge with nonclonal L1 Dist
  df_final<- merge(df_final,CV_nonclonal_res, by=0)
  rownames(df_final)<- df_final$Row.names; df_final<- df_final[,-1]
  colnames(df_final)<- c("Expr_Clonal", "Expr_NonClonal", "Expr_CminusNC", "Expr_Zscore",
                         "L1Dist_Clonal", "L1Dist_NonClonal")
  
  ## Calculate CV for clonal/nonclonal
  df_final$CV_clonal<- df_final$L1Dist_Clonal/df_final$Expr_Clonal
  df_final$CV_NonClonal<-df_final$L1Dist_NonClonal / df_final$Expr_NonClonal
  
  ## Calculate CV zscore
  df_final<- df_final[!is.na(df_final$CV_NonClonal),]
  df_final$CV_CminusNC<- df_final$CV_clonal - df_final$CV_NonClonal
  df_final$CV_Zscore<- (df_final$CV_clonal - df_final$CV_NonClonal) / sd(df_final$CV_NonClonal)
  
  return(df_final)
}




