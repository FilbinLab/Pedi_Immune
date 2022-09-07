
## Preprocessing immune cell data

*Note:*

*Tracer + clonal analysis was run for most immune samples by Orr Ashenberg (Regev lab). Not repeated by me. I use the output file for # clonal cells, clonotype, etc, from this analysis.  Output is in: data/data_fromOrr/2020_12_10_pedbrain_immune_no_batchcorrect.Rda and 2020_12_10_clonalTcells_cluster.Rda*

*For 5 samples (BT1873/BT1857/MUV92/MUV94/MUV95), tracer analysis was run by Filbin lab. For these samples, I ran all post-processing to determine clonal cells, clonotype, etc*

*Except for 03_PreprocessAdultDatasets directory, all files deal with Filbin pediatric samples ONLY. In downstream analyses, ped/adult are merged. Here, pediatric/adult are kept separate for now.*



# Pediatric immune  preprocessing

1. Tracer preprocessing for 5 new samples: 01_Tracer_Preprocess (Note: need to run 02_Counts_Preprocess/01_5Samples_ProcessToSeurat.Rmd just to have seurat object)\

2. Finish processing for these 5 samples: 02_Counts_Preprocess/01_5Samples_ProcessToSeurat.Rmd\

3. Merge these 5 samples with results from Orr for rest of cohort: Preprocessing/02_Counts_Preprocess/02_MergeSamples_QC_ToSeurat.Rmd (Note: have 4 cohort options here: allSamples = all filbin lab samples; pedOnly=samples originally thought to be ped; pedOnly_nomuv63=updated ped samples; pedOnly_nomuv63.withmu91=ped samples updated again. *pedOnly_nomuv63 used downstream*) In addition, have the option to use harmony integration or not. 

4. For ALL samples, cells are split by "broad" annotation (Myeloid vs Tcell) in Preprocessing/01_Counts_Preprocess/03_MergedSamples_BroadAnnotation.Rmd. *For "broad" annotation, non-integrated ped seurat object is used. This is the only time during entire final analysis that non-harmony-integrated seurat objects are used. Both (with and without) are shown here, but "without integration"-based annotations for myeloid/tcells are used downstream. For myeloid/tcell, integration doesn't seem to really matter as these populations are quite distinct. (Note: there is a small subset of myeloid cells that cluster with Tcells. Ran NMF on all Tcells at this point (originally just for doublet detection, based on methods in Cell 2021) and used this to identify. This NMF (In 02a_Tcells/NMF/NMF_DoubletID.RMd) needs to be run prior to finishing broad annotation file.)*

5. Originally, pediatric cells are then split by "detailed" annotation (CD4/CD8/Cycling for Tcells, DC/BCells/Myeloid for Myeloid) in 02_Counts_Preprocess/04_MergedSamples_DetailedAnnotation.Rmd. IN FINAL ANALYSIS, THESE ANNOTATIONS NO LONGER USED. Changed as of 3/22. Reasoning: CD4/CD8 annotations are not obvious and have a big effect on downstream NMF. Relatively small sample size of ped cohort also affects CD4/CD8 annotation. Instead, "detailed" annotation is as follows:
	1. Run NMF on pediatric T cells (ALL TCELLS) and adult (from Cell 2021 paper, NOT for filbin adult samples) Tcells (In *.Rmd)
	2. Identify "shared" and "specific" programs, re-generate genes from "shared" programs
	3. Using these program assignments, determine CD4/CD8 T cells separately for each NMF program, for both ped/adult simultaneously (in *.Rmd)

	Using this approach improves correlation between adult/ped cell programs, slightly improves accuracy of cd4/cd8 annotation, and increases our confidence in these annotations due to larger cell count provided by adult cells

	For myeloid cells, the original "detailed" annotations are used ( in 02_Counts_Preprocess/04_MergedSamples_DetailedAnnotation.Rmd.). This file edited to clarify that this approach used only for myeloid cells. (Note: no harmony integration is used for this, but these cell types are quite distinct and it doesn't seem to matter much)


Additional files: CompareSeuratPipelines.Rmd used for testing our "standard" pipeline vs Orr's vs integration

# Adult immune dataset preprocessing

Performed in 03_PreprocessAdultDatasets. Contains preprocessing for myeloid and cell adult datasets: 

	1. Myeloid: preprocessing for 3 publicly available myeloid datasets in Preprocess_AdultMyeloid. These datasets are: IDHmut (Science), GBM (Neftel), GBM (10X, Movahedi)

All 3 datasets are processed to seurat objects, any annotations as needed, any re-scaling as needed. Saved for comparison to pediatric myeloid cells downstream

Note: unlike for Tcells, downstream "detailed" annotation is done SEPARATELY for pediatric/adult. Any DC/Bcells (rare) are removed from adult myeloid datasets here, using same method as for pediatric

	2. Tcells: adult t cells used downstream are all from Cell 2021 and were provided as a seurat object by Orr. In Tcells_Ped.Adult_Merge, pediatric and adult cells are:

	1) merged into single seurat object 
	2) harmony integration 
	3) CD4/CD8 annotations- BUT these annotations are NOT used downstream. Compared to detailed annotation by NMF program, this method was slightly less accurate so was not used downstream.
