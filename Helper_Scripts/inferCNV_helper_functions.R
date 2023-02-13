library(infercnv)

addNormalControl <- function(ext_ctrl, inferCNV_analysis_folder, 
                             cm_raw, orig_samples, 
                             normal_cm_name="K27M_20190124normal.expr.txt"){
  message("Loading normal control count matrix...")
  ext_ctrl_cm = read.table(paste0(ext_ctrl, normal_cm_name), sep = "\t", header=T, stringsAsFactors = F, row.names = 1)
  ext_ctrl_cm = subset(ext_ctrl_cm, select=-c(id, coding))
  colnames(ext_ctrl_cm) = cleanCellName(colnames(ext_ctrl_cm))
  
  ## Read and process id of the ctrl data
  message("Loading normal control annotations...")
  microglia_cell = readLines(paste0(ext_ctrl, "Microglia.txt"))
  oc_cell = readLines(paste0(ext_ctrl, "Oligodendrocytes.txt"))
  microglia_cell = cleanCellName(microglia_cell)
  oc_cell = cleanCellName(oc_cell)
  normal_samples = c(rep("mg", length(microglia_cell)), rep("od", length(oc_cell)))
  names(normal_samples) = c(microglia_cell, oc_cell)
  
  ## Reorder the ext_ctrl cm and subset genes
  message("Reordering cell names for normal control cm and annotations")
  ext_ctrl_cm = ext_ctrl_cm[,names(normal_samples)]
  ext_ctrl_cm = ext_ctrl_cm[rownames(cm_raw),]
  
  ## Concatenate samples and normal controls 
  message("Concatenate cm and annotations of normal controls to patient data...")
  cm = cbind(cm_raw, ext_ctrl_cm)
  samples = c(orig_samples, normal_samples)
  message("Saving data...")
  saveRDS(cm, paste0(inferCNV_analysis_folder, "cm_exp_ctrl.rds"))
  saveRDS(samples, paste0(inferCNV_analysis_folder, "samples_exp_ctrl.rds"))
}

addNormalControlNuc <- function(ext_ctrl, inferCNV_analysis_folder, cm_raw, 
                                orig_samples, premrna = F){
  ## Read ctrl cm 
  message("Loading normal control count matrix...")
  if (premrna){
    suffix = "_premrna.rds"
  }else{
    suffix = ".rds"
  }
  mg_cm = readRDS(paste0(ext_ctrl, "mg_nuc", suffix))
  ##tc_cm = readRDS(paste0(ext_ctrl, "tc_nuc", suffix))
  od_cm = readRDS(paste0(ext_ctrl, "od_nuc", suffix))
  ext_ctrl_cm = cbind(mg_cm, od_cm)
  
  ## Generate control sample list 
  message("Loading normal control annotations...")
  normal_samples = c(rep("mg", dim(mg_cm)[2]), rep("od", dim(od_cm)[2]))
  names(normal_samples) = c(colnames(mg_cm), colnames(od_cm))
  
  ## Reorder the ext_ctrl cm and subset genes
  message("Reordering cell names for normal control cm and annotations")
  ext_ctrl_cm = ext_ctrl_cm[,names(normal_samples)]
  ext_ctrl_cm = ext_ctrl_cm[rownames(cm_raw),]
  
  ## Concatenate samples and normal controls 
  message("Concatenate cm and annotations of normal controls to patient data...")
  cm = cbind(cm_raw, ext_ctrl_cm)
  samples = c(orig_samples, normal_samples)
  message("Saving data...")
  saveRDS(cm, paste0(inferCNV_analysis_folder, "cm_exp_ctrl.rds"))
  saveRDS(samples, paste0(inferCNV_analysis_folder, "samples_exp_ctrl.rds"))
}

## Compute score for oligodendrocytes or immune cells 
## @param cm_norm count matrix with log transformed but uncentered data 
## @param cm_center count matrix with log transformed and centered data
## @param macrophage_marker marker genes for microglias
## @param t_cell_marker marker genes for T cells
## @param oc_marker markger genes for oligodendrocytes 
## @param simple whether use simple average, default = FALSE
## @return a df with scores for these cell types 
## How to deal with marker genes not found in this dataset?
computeNormalCellScore <- function(cm_norm, cm_center, gene_mean, scoreName=c("microglia","t_cell","oc"),
                                   microglia_marker=c("CD14", "AIF1", "FCER1G", "FCGR3A", "TROBP", "CSF1R"),
                                   t_cell_marker=c("CD2", "CD3D", "CD3E", "CD3G"),
                                   oc_marker=c("MBP", "TF", "PLP1", "MAG", "MOG", "CLDN11"),
                                   simple=FALSE,
                                   verbose=TRUE){
  ## Check if raw and normalized data have the same cell names
  try (if(any(rownames(cm_norm) != rownames(cm_center))) stop("normalized and centered cm should have identical cell names and orders"))
  
  ## Compute mean gene expression
  ##message("Compute mean gene expressions...")
  ##gene_mean = rowMeans(cm_norm)
  
  ## Compute each cell score
  message("Compute single cell scores...")
  microglia_score = scoreSignature(cm_center, gene_mean, microglia_marker, n=100, simple=simple, verbose=verbose)
  t_cell_score = scoreSignature(cm_center, gene_mean, t_cell_marker, n=100, simple=simple, verbose=verbose)
  oc_score = scoreSignature(cm_center, gene_mean, oc_marker, n=100, simple=simple, verbose=verbose)
  
  ## cbind all scores into a df and result
  message("Concatenate scores into a df...")
  result = data.frame(cbind(microglia_score, t_cell_score, oc_score))
  colnames(result) = c("mg", "tc", "od")
  rownames(result) = colnames(cm_norm)
  return(result)
}

## Add sample annotation 
## @param sc_scores single cell scores for microglia, t-cells and oligodendrocytes
## @param sample_names the vector of all patient/sample names
## @param threshold single cell threshold to call a cell any cell type 
## @returns a df with rownames = cell_name and one column = malignancy status 
generateMalignancyAnnotation <- function(sc_scores, sample_names, useCellType = F, threshold = 4){
  ## Check cell names and orders match
  try (if(any(rownames(sc_scores) != names(sample_names))) stop("sc scores and cell names should have identical cell names and orders"))
  ## Assign maligancy status based on scores and threshold 
  if (useCellType){
    sc_scores$malig_status = ifelse((apply(sc_scores, 1, max) < threshold),
                                    paste("malignant", sample_names, sep="_"),
                                    apply(sc_scores, 1, function(x) names(x)[x==max(x)])
    )
  }
  else{
    sc_scores$malig_status = ifelse((apply(sc_scores, 1, max) >= threshold), 
                                    "non-malignant", 
                                    paste("malignant", sample_names, sep="_"))
  }
  
  result = sc_scores[,c("malig_status")]
  names(result) = rownames(sc_scores)
  return(data.frame(result))
}

computeCnvScore <- function(cnv_values, squared=T){
  if (squared){
    return (colMeans(cnv_values^2))
  }
  else{
    return (colMeans(abs(cnv_values)))
  }
}

## Compute sample average CNV values
## @param cnv_values cnv_values from outputs of inferCNV
## @param samples a vector of sample names
## @param malig_status malignancy status of each cell (reuse inputs for inferCNV) 
## @param min_tumor_cell minium number of tumor cells required to calculate sample mean CNV value 
calculateAverageCnvValue <- function(cnv_values, samples, malig_status, min_tumor_cell=3){
  try (if(any(colnames(cnv_values) != names(samples) | any(colnames(cnv_values) != rownames(malig_status))))
    stop("CNV values, sample names, and malignancy status should have identical cell names and orders"))
  
  result = list()
  unique_sample_names = names(table(samples))
  message("Subsetting cnv_values and samples with only expressionally defined tumor cells...")
  ##cnv_values_tumor = cnv_values[,which(malig_status$result != "non-malignant")]
  putative_malig = sapply(malig_status, function(x) grepl("^malignant", x))
  cnv_values_tumor = cnv_values[,putative_malig]
  ##samples = samples[malig_status$result != "non-malignant"]
  samples = samples[putative_malig]
  
  try(if(dim(cnv_values_tumor)[2] != length(samples)) stop("Subsetted cnv values and sample names should have equal length"))
  
  message("Computing average CNV values for each sample...")  
  for (sample_name in unique_sample_names){
    tmp_cnv_values = cnv_values_tumor[,which(samples == sample_name)]
    message("The number of putative tumor cells in sample ", sample_name, " is ", dim(tmp_cnv_values)[2])
    if (dim(tmp_cnv_values)[2] >= min_tumor_cell){
      result[[sample_name]] = rowMeans(tmp_cnv_values)
    }
    else{
      result[[sample_name]] = NA
    }
  }
  return(result)
}

## Compute correlation between CNV values of each cell and average CNV values of all putative tumor cells from each sample 
## @param cnv_values cnv_values from outputs of inferCNV
## @param mean_cnv_values average CNV values of all putative tumor cells from a sample
## @param samples a vector of sample names
calculateCnvCor <- function(cnv_values, mean_cnv_values, samples){
  try (if(any(colnames(cnv_values) != names(samples)))
    stop("CNV values and sample names should have identical cell names and orders"))
  
  cor_coef = NULL
  for (i in seq(dim(cnv_values)[2])){
    sample = as.character(samples[i])
    #print(sample)
    #print(head(mean_cnv_values$sample))
    if (sum(is.na(mean_cnv_values[[sample]]))>0){
      cor_coef = c(cor_coef, 0)
    }else{ 
      ##cor_coef = c(cor_coef, cor(cnv_values[,i], mean_cnv_values[[sample]], method = "spearman"))
      cor_obj = cor.test(cnv_values[,i], mean_cnv_values[[sample]], 
                         alternative = "two.sided", method = "spearman", exact = FALSE)
      if (cor_obj$p.value > 0.05){
        cor_coef = c(cor_coef, 0)
      }
      else{
        cor_coef = c(cor_coef, cor_obj$estimate)
      }
    }
  }
  
  names(cor_coef) = names(samples)
  return(cor_coef)
}

## Distinguish malignant and non-malignant cells
## @param cnv_score a vector of cnv scores (output from calculateAverageCnvValue)
## @param cor_coef a vector of correlation coefficient (output from calculateCnvCor)
callTumorCell <- function(cnv_values, samples, malig_status, cnv_score_cutoff=0.03, cnv_coef_cutoff=0.2){
  ## Calculate CNV score
  message("Calculating CNV scores...")
  cnv_score = colMeans(cnv_values^2)
  ## Calculate CNV correlation 
  message("Calculating CNV correlations...")
  mean_cnv_value = calculateAverageCnvValue(cnv_values, samples, malig_status)
  cor_coef = calculateCnvCor(cnv_values, mean_cnv_value, samples)
  ## Call malignant cells
  message("Call malignant cells...")
  result = ifelse((cnv_score < cnv_score_cutoff & cor_coef < cnv_coef_cutoff), "non-malignant", "malignant")
  names(result) = names(samples)
  return(result)
}

plot_subcluster_10X <- function(folder, sample, cm, samples, normal_type, cutree_res,
                                gene_order_file, cutoff = 0.1, window_length = 201, HMM=FALSE){
  ## Create dir to store results
  if (!dir.exists(folder)){
    dir.create(folder)
  }
  
  ## store cell names in each cluster to a list
  group_num = unique(cutree_res[[sample]])
  group_name = lapply(group_num, function(x){
    tmp = names(cutree_res[[sample]][cutree_res[[sample]] == x])
    tmp = tmp[tmp %in% names(samples)[samples == sample]]})
  names(group_name) = group_num
  group_name = group_name[sapply(group_name, length)>1]
  
  
  cm_sub = cm[,c(unlist(group_name), names(samples)[samples %in% normal_type])]
  label = c(unlist(sapply(names(group_name), function(x) rep(paste0("malignant_", x), length(group_name[[x]])))), 
            unlist(sapply(normal_type, function(x) rep(x, sum(samples==x)))))
  names(label) = colnames(cm_sub)
  
  write.table(as.matrix(cm_sub), file = paste0(folder, "counts.matrix"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(label, file = paste0(folder, "cellAnnotations.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(folder, "counts.matrix"),
                                      annotations_file=paste0(folder, "cellAnnotations.txt"),
                                      delim="\t",
                                      gene_order_file=gene_order_file,
                                      ref_group_names=normal_type)
  
  # perform infercnv operations to reveal cnv signal
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=cutoff,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir=paste0(folder, "out_dir"),  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster
                               hclust_method="ward.D2",
                               window_length = window_length,
                               plot_steps = F,
                               denoise=T,
                               HMM=HMM
  )
  
}

call_cnv_hc <- function(cnv_cutree_table, cnv_cutree_each_sample, samples, obs_ref, inferCNV_analysis_folder ,norm_cutoff=3){
  ## Identify clusters of normal cells in each cluster 
  ## Defined as group of samples that clusters with >= norm_cutoff normal control cells
  normal_clusters = lapply(cnv_cutree_table, function(x) which((x["mg",] + x["od",])>=norm_cutoff))
  try (if(any(names(cnv_cutree_each_sample) != names(normal_clusters))) stop("cutree results list and normal clusters list should have identical names for elements"))
  
  ## Assign CNV to cells based on whether cells in normal clusters or not 
  call_cnv_res = NULL
  for (i in seq(length(normal_clusters))){
    ## Identify current sample and subset normal clusters and cutree results for that sample 
    sample = names(normal_clusters)[i]
    normal_clust = normal_clusters[[i]]
    hc_clust_res = cnv_cutree_each_sample[[i]]
    ## Note: this way of extracting sample names exclude normal controls from the final result 
    sample_names = sapply(names(hc_clust_res), function(x) unlist(strsplit(x, split = ".", fixed = T))[1])
    ## For each sample, only keep cells in that sample, not normal controls 
    hc_clust_res = hc_clust_res[sample_names == sample]
    tmp = ifelse(hc_clust_res %in% normal_clust, F, T)
    names(tmp) = names(hc_clust_res)[sample_names == sample]
    call_cnv_res = c(call_cnv_res, tmp)
  }
  
  ## Make a vector of mg and od; they are both normal cells 
  mg_names = names(samples)[samples == "mg"]
  mg_res = rep(F, length(mg_names))
  names(mg_res) = mg_names
  od_names = names(samples)[samples == "od"]
  od_res = rep(F, length(od_names))
  names(od_res) = od_names
  
  ## Concatenate samples and normal controls into a single vector 
  call_cnv_res_complete = c(call_cnv_res, mg_res, od_res)
  call_cnv_res_complete = call_cnv_res_complete[names(samples)]
  call_cnv_res = call_cnv_res_complete[names(samples)[obs_ref == "obs"]]
  
  ## Save CNV results (w/o normal controls)
  saveRDS(call_cnv_res, paste0(inferCNV_analysis_folder, "call_cnv.rds"))
  saveRDS(call_cnv_res_complete, paste0(inferCNV_analysis_folder, "call_cnv_w_ctrl.rds"))
}

# plot_subcluster <- function(folder, sample, group_1, group_2, group_1_num, group_2_num,
#                             cm, samples, gene_order_file, cnv_hc_each_sample, window_length = 201){
#     if (!dir.exists(folder)){
#         dir.create(folder)
#     }
#     
#     group_1 = cnv_hc_each_sample[[sample]]$labels[cnv_cutree_each_sample[[sample]] == group_1][1:group_1_num]
#     group_2 = cnv_hc_each_sample[[sample]]$labels[cnv_cutree_each_sample[[sample]] == group_2][1:group_2_num]
#     
#     cm = cm[,c(group_1, group_2, names(samples)[samples=="od" | samples == "mg"])]
#     label = c(rep("malignant_1", length(group_1)), rep("malignant_2", length(group_2)), rep("mg", 96), rep("od", 94))
#     names(label) = colnames(cm)
#     
#     write.table(cm, file = paste0(folder, "counts.matrix"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
#     write.table(label, file = paste0(folder, "cellAnnotations.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
#     
#     infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(folder, "counts.matrix"),
#                                         annotations_file=paste0(folder, "cellAnnotations.txt"),
#                                         delim="\t",
#                                         gene_order_file=gene_order_file,
#                                         ref_group_names=c("mg", "od"))
#     
#     # perform infercnv operations to reveal cnv signal
#     infercnv_obj = infercnv::run(infercnv_obj,
#                                  cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
#                                  out_dir=paste0(folder, "out_dir"),  # dir is auto-created for storing outputs
#                                  cluster_by_groups=T,   # cluster
#                                  hclust_method="ward.D2",
#                                  window_length = window_length,
#                                  plot_steps = F,
#                                  include.spike=T
#     )
# }

# plot_subcluster_nuc <- function(folder, sample, group_1, group_2, group_1_num, group_2_num,
#                                 cm, samples, gene_order_file, cnv_hc_each_sample, window_length = 201){
#     if (!dir.exists(folder)){
#         dir.create(folder)
#     }
#     
#     group_1 = cnv_hc_each_sample[[sample]]$labels[cnv_cutree_each_sample[[sample]] == group_1][1:group_1_num]
#     group_2 = cnv_hc_each_sample[[sample]]$labels[cnv_cutree_each_sample[[sample]] == group_2][1:group_2_num]
#     
#     cm = cm[,c(group_1, group_2, names(samples)[samples=="od" | samples == "tc" | samples == "mg"])]
#     label = c(rep("malignant_1", length(group_1)), rep("malignant_2", length(group_2)), 
#               rep("mg", 243), rep("tc", 22), rep("od", 38))
#     names(label) = colnames(cm)
#     
#     write.table(cm, file = paste0(folder, "counts.matrix"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
#     write.table(label, file = paste0(folder, "cellAnnotations.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
#     
#     infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(folder, "counts.matrix"),
#                                         annotations_file=paste0(folder, "cellAnnotations.txt"),
#                                         delim="\t",
#                                         gene_order_file=gene_order_file,
#                                         ref_group_names=c("mg", "tc", "od"))
#     
#     # perform infercnv operations to reveal cnv signal
#     infercnv_obj = infercnv::run(infercnv_obj,
#                                  cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
#                                  out_dir=paste0(folder, "out_dir"),  # dir is auto-created for storing outputs
#                                  cluster_by_groups=T,   # cluster
#                                  hclust_method="ward.D2",
#                                  window_length = window_length,
#                                  plot_steps = F,
#                                  include.spike=T
#     )
#     
#     infercnv::plot_cnv(infercnv_obj,
#                        out_dir=paste0(folder, "out_dir"),
#                        cluster_by_groups=T,
#                        color_safe_pal=FALSE,
#                        x.center=1,
#                        x.range=c(0.5,1.5),
#                        title="inferCNV",
#                        obs_title="Observations (Cells)",
#                        ref_title="References (Cells)",
#                        output_filename="heatmap",
#                        hclust_method = "ward.D2",
#                        write_expr_matrix = F
#     )
# }