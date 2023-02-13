calculatePairwiseCor <- function(ref_cm_norm, pb_cm_norm, input_genes){
  ## Subset both cm to have identical gene names 
  shared_genes = intersect(rownames(ref_cm_norm), rownames(pb_cm_norm))
  shared_genes = intersect(shared_genes, input_genes)
  pb_cm_norm = pb_cm_norm[shared_genes, ]
  ref_cm_norm = ref_cm_norm[shared_genes, ]
  
  ## Correlation
  pairwise_cor = cor(ref_cm_norm, pb_cm_norm)
  return(pairwise_cor)
}