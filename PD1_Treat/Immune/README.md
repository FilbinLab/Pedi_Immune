## Processing and annotation of immune single cells from PD1 treated patient

*Preprocess_QC_ToSeurat.Rmd*: QC filtration and initial graph-based clustering, seurat object generation

*PostTracer_Preprocessing.Rmd*: Tracer used to identify clonally expanded T cells (run on Broad server, not locally run). Process results here for integration into downstream analysis

*Preprocess_BroadAnnotation.Rmd*: Annotation of myeloid/Tcells by marker gene expression, split into separate seurat objects and re-clustered

*Preprocess_DetailedAnnotation.Rmd*: Identify CD4/CD8 within T cells, DC and Bcells within Myeloid cells

*CompareToFullCohort.Rmd*: Score for programs identified in full cohort, compare broad proportions/expression patterns.
