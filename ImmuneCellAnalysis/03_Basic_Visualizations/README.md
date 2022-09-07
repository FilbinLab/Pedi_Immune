## Integration of immune and tumor scRNAseq

Two main approaches for comparing tumor/immune data:

1. CellChat (currently working to rerun this on finalized annoys, will update later)

2. Correlating proportion of programs/scores in tumor and immune
	This is done separately for Tcell/Myeloid and Tcell/Tumor + Myeloid/Tumor
	Essentially the same analyses run in each:
		1. Correlate proportion of programs- linear regression + pearson correlation
		2. Correlate scores of programs- linear regression + pearson correlation
