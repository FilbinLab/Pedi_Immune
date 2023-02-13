
## Two pediatric glioma patients treated with neoadjuvant anti-PD1 and profiled by scRNAseq analyzed here and compared to their untreated counterparts

Both patients had 2 recurrences, both of which were profiled, for a total of 4 samples.

### For SS2 (01a_SS2/Immune):
*Patient 1 profiled by SS2*

1. 01a_Preprocess_QC_ToSeurat: standard qc, seurat pipeline 
2. 01b_PostTracer_Processing: process results from tracer, basic clonality plots. add info to seurat.
3. 02_ProjectBroad.DetailedAnnots: project myeloid/tcell, bc/dc, cd4/cd8 from full immune to 2 pd1 treated patients
4. 03_MergeGBM.SS2: merge all ss2 gbm- treated and untreated- into single seurat object
5. 04_ProjectPrograms: project program annotations for myeloid/cd4/cd8 from full cohort to merged ss2 gbm (treated and untreated).
6. 05_ProgramDE: run DE program-wise for pd1 vs untreated


Note that also did program projections for pd1 alone. As of 12/7/22, decided to use projections to ALL gbm rather than pd1 only- just to increase confidence in changes in proportion.

### For 10x (/01b_tenX/):
*Patient 2 profiled by 10X*

same general processing as for SS2:

1. QC, broad/detailed projection separate for treated and untreated
2. Program projection together for treated and untreated
3. Program DE between treated/untreated



### Merge ss2/10x into single seurat just for plotting changes in proportion, GOI expression, DE results, etc (/02_Merge_SS2.tenX/)
1. Merge all GBM into single seurat object (treated, untreated; SS2, TenX) (01_MergeGBM.Rmd)
2. Plot changes in program proportions, other basic vis, by treatment (02_PlotProgramProportions.Rmd)
3. Overlap b/w DE genes by program, separately for tenx/ss2 (03a_IdentifySharedDEGs.Rmd)
4. Expression of genes of interest (03b_GOIExpression.Rmd)
5. Changes in clonality, clonal group, etc by program/treatment (03c_ClonalAnalysis.Rmd)

### Compare to published adult glioma anti-pd1 treated data(/03_comparisonToAdult/)
*Data taken from GSE154795*
1. process adult data, annotate as needed (01_PreprocessAdult.Rmd)
2. project/score ped programs for adult (02a/b)
3. Combine ped/adult pd1 treated into single seurat object (03_MergePedAdult.Rmd)
4. Expression of genes of interest (04_GOIExpression.Rmd)

### Compare to published melanoma anti-PD1 treated data (04_comparisonToMelanoma/)
*Data taken from GSE120575*
1. Process melanoma data, annotate as Tcell, myeloid, CD4, CD8 (01_PreprocesMelanoma.Rmd)
2. Score for CD8 signatures to annotate (01_PreprocessMelanoma.Rmd)
3. Assign Tregs based on marker gene expression in CD4 (01_PreprocessMelanoma.Rmd)
4. Plot changes in programs, gene of interest expression (02_PlotChanges.Rmd)

