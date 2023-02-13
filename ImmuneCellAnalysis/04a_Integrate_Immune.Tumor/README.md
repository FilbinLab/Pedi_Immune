## Integration of immune and tumor scRNAseq

Two main approaches for comparing tumor/immune data:

1. ***CellChat/***: Identify receptor/ligand pairs using CellChat.
- ***01_PrepForCellChat.Rmd***: Create separate seurat objects for each sample
- ***02_CellChat_ForCluster.Rmd***: Run CellChat (run on O2 server for computational considerations- see scripts in /CellChat_HPC_scripts/)
- ***03_CellChat_Preprocess.Rmds***: Process cellchat object from cellchat for better visualizations- extract LR pairs
- ***04a/b_CellChat_Visualizations_*.Rmd***: Plot top ligand/receptor pairs between cell types (CD4, CD8, AC-like, myeloid, etc) OR programs (Cytotoxic CD8, Treg, AC-like, IFN_TAM, etc)

2. ***Correlate_Immune.Tumor_CellType.NMF***: Correlating proportion of programs/scores in tumor and immune

	This is done separately for Tcell/Myeloid and Tcell/Tumor + Myeloid/Tumor

	Essentially the same analyses run in each:
	1. Correlate proportion of programs- linear regression + pearson correlation
	2. Correlate scores of programs- linear regression + pearson correlation
