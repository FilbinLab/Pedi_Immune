## All analyses for immune scRNAseq data

***01_Preprocessing/***: Initial QC performed by Orr Ashenberg (Regev lab, Broad Institute). Processing here: 
  - Post-**Tracer** processing for clonal T cell identifiction
  - **merging samples** (sequenced and added to cohort at different times)
  - **"Broad" annotation** (Tcell vs Myeloid)
  - **"Detailed" annotation** (CD4/CD8, DC/BCell/Myeloid)
  - preprocessing of **published adult** data
 
***02a_Tcells/***: analysis of T cells- pediatric and adult NMF, comparisons across age *(see README in directory for details)*

***02b_Myeloid/***: analysis of myeliod cells- pediatric and adult NMF, comparisons across age *(see README in directory for details)*

***03_Basic_Visualizations/***: Nice plots of broad/detailed annotation, sample, subtype, etc. Plotting genes of interest.

***04a_Integrate_Immune.Tumor/***: approaches to "integrating" across immune tumor: *(see READMEs in directory for details)*
  - **Correlate_Immune.Tumor_CellType.NMF**: correlating program proportions and scores in sample-wise manner
  - **CellChat**: Identification of receptor/ligand pairs across tumor cell types, myeloid programs, T cell programs
  -**SCENIC**: Transcription factor regulatory network analysis across myeloid/T cell programs

***04b_DeconvoluteBulk/***: Deconvolution of bulk RNAseq samples using single cell programs, probing for survival effect *(see README in directory for details)*

