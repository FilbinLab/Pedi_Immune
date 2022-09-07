## Deconvolute bulk RNAseq data using scRNAseq programs and probe for survival effect

#### Preprocessing (generating cibersort reference, processing bulk data, pseudobulking scRNAseq for validation) performed in ***01_Preprocessing.Rmd*** 

#### Two main deconvolution approache tested using pseudobulked scRNAseq data:
1. Cibersort- ***02a_Cibersort_validate.Rmd*** *(used downstream)*
2. Linear regression ***02b_lmDeconvolution_validate.Rmd*** *(not used downstream)*
 

#### Cibersort used downstream to deconvolute actual bulk data- ***03_DeconvoluteBulk_Preprocess.Rmd***

#### Two main survival approaches *(both used downstream)*
1. Kaplan meier: ***04a_KM_SurvivalAnalysis.Rmd***
2. Cox regression to control for confounders like subtype, histone mutation status, etc: ***04b_MultivariateCox_fromtemplate.Rmd***
