
##############
## PLOTTING ##
##############

## Plot heatmap of markers for each kmeans cluster
## To check how well clusters are stratified
## Input: cm = count matrix, genes x samples
##        marker_number = index of gene set to plot
##        cluster_name = name of geneset to plot in metadata
##        meta = metadata with scores, samples
##        marker_list = list of marker sets
plotHeatmapOfMarkers<-function(cm, marker_number, cluster_name, meta, marker_list){
  m<- names(marker_list)[[marker_number]]
  cm_subset<- cm[rownames(cm) %in% unname(unlist(marker_list[m])),]
  
  ## If gene list is just 1 gene: edit cm_subset for correct input
  if(length(unname(unlist(marker_list[m])))==1){
    cm_subset<- as.data.frame(cm_subset)
    cm_subset$ID<- rownames(cm_subset); cm_subset$Gene<- unname(unlist(marker_list[m]))
    colnames(cm_subset)<- c("counts", "ID", "Gene")
    df<- cm_subset
  } else{
    df<- (melt(t(cm_subset)))
    colnames(df)<-c("ID", "Gene", "counts")}
  
  ## df for heatmap input
  df2<-meta[,c("submitter_id", cluster_name)]; colnames(df2)<- c("ID", "Markers")
  df3<- merge(df, df2, by="ID")
  
  ## Plot heatmap
  ggplot(df3, aes(x=ID, y=Gene, fill=counts))+ 
    geom_tile()+
    scale_fill_viridis(discrete=FALSE)+
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    facet_grid(.~Markers,scales = "free_x")
}

#######################
## SURVIVAL ANALYSIS ## 
#######################

## Wrapper function for running survival analysis based on high/low marker split
## Input: input_df = data frame with: time=numeric column of time to death, censored=logical column with dead=1/alive=0, DepColumnName (geneset used to split)
##        DepColumnName = name of column with dependent variable to split survival analysis on- ex: High/Low expression of gene set
run_mySurvAnalysis<- function(input_df, DepColumnName){
  ## Run survival based on dependent variable
  input_df$Marker<-input_df[[DepColumnName]]
  fit <- survfit(Surv(time, censored) ~Marker,
                 data = input_df)
  
  ## Create dataframe of survival results
  pvalue<-surv_pvalue(fit, input_df)$pval.txt
  pvalue<-gsub("p = ", "", pvalue)
  myFit_median<-surv_median(fit)
  highMarker_median<-myFit_median[myFit_median$strata=="Marker=High", "median"]
  lowMarker_median<-myFit_median[myFit_median$strata=="Marker=Low", "median"]
  df<-data.frame(Marker=DepColumnName, pvalue=pvalue,
                 HighMarkerMedianSurvival=highMarker_median,
                 LowMarkerMedianSurvival=lowMarker_median)
  
  ## Return list of results
  res_list<- list(fit=fit, res_df=df )
  return(res_list)
}

## Visualize survival results with survminer
## Input: input_df = data frame with: time=numeric column of time to death, censored=logical column with dead=1/alive=0, DepColumnName (geneset used to split)
##        Surv_fit = results from running survfit(Surv(time, censored)~Dep, data=input_df). Output of run_mySurvAnalysis (res_list$fit)
##        DepColumnName = name of column with dependent variable to split survival analysis on- ex: High/Low expression of gene set
##        plotName = plot title. If NULL (default), just uses DepColumnName
plot_mySurvAnalysis<- function(surv_fit,input_df, DepColumnName, plotName=NULL, showPvalue=TRUE){
  ## Add on "Marker" to input df- so that it matches exactly from survival analysis
  input_df$Marker<-input_df[[DepColumnName]]
  
  plotTitle<- ifelse(is.null(plotName), DepColumnName, plotName )
  p<-ggsurvplot(surv_fit, pval=showPvalue, input_df)+ggtitle(plotTitle)
  return(p)
}


## Run multivariate (or univariate if controlFor=FALSE) cox regression
## Input: myVariable = variable name to plot in metadata
##        AllMeta = "finalMeta"- contains time, censored, and at least myVariable
##        controlFor = whether or not this will be univariate (FALSE) or multivariate (TRUE)
##        mySubset = which samples used in analysis- only used for saving analysis. useful when trying out multiple subsets
RunCox<- function(myVariable, AllMeta, controlFor=FALSE,mySubset="AllSamples"){
  if (controlFor[1] == FALSE){
    message("Running univariate...")
    df<- AllMeta[,c("time", "censored", myVariable)]
    runCox<- paste0("coxph(Surv(time, censored) ~", myVariable, ", data=df)")
    message(paste0("Cox formula:", runCox))
    cox<- eval(parse(text=runCox))
  }
    else if (length(controlFor)==1){
    message("Running multivariate: 1 variable to control for...")
    df<-AllMeta[,c("time", "censored", controlFor, myVariable)]
    runCox<- paste0("coxph(Surv(time, censored) ~", myVariable, "+ ", controlFor, ", data=df)")
    message(paste0("Cox formula:", runCox))
    cox<- eval(parse(text=runCox))
  }
    
   else{
    message("Running multivariate: >1 variables to control for")
    df<-AllMeta[,c("time", "censored", controlFor, myVariable)]
    controllingFor<-paste(controlFor, collapse="+")
    runCox<- paste0("coxph(Surv(time, censored) ~", myVariable, "+ ", controllingFor, ", data=df)")
    message(paste0("Cox formula:", runCox))
    cox<- eval(parse(text=runCox))
  }
  return(cox)
}

## Create df of univariate cox results for export- potential variables
## Input: input_results = cox results from RunCox
##        cohort/age = used for saving data, may not be necessary when only looking at a couple of cohorts
WriteOutCoxResults<- function(input_results, cohort="", age=""){
  cox_results<-as.data.frame(summary(input_results[[1]])$coefficients)
  cox_results$Variable<-names(input_results)[1]
  for (i in 2:length(input_results)){
    y<-input_results[[i]]
    coef<- as.data.frame(summary(y)$coefficients)
    coef$Variable<-names(input_results)[i]
    cox_results<- rbind(cox_results, coef)
  }
  return(cox_results)
}

## Wrapper function that applies RunCox function (and WriteOutCoxResults) to multiple genesets automatically
## Input: geneSet_groups = vector of genesets names to apply wrapper to
##        AllMeta = "FinalMeta"- contains time, censored, and all genesets in geneSet_groups
##        controlFor = whether or not this will be univariate (FALSE) or multivariate (TRUE)
##        mySubset = which samples used in analysis- only used for saving analysis. useful when trying out multiple subsets
##        exportName = analysis name to use when saving data
##        exportCsv = whether or not to export csv file of results
##        subsetToHighLow = when splitting 25/50/25, subset to just 25/25
##        resultsFolder = folder to export results to
##        appendToNumericVariable = if variable is numeric, append this to geneset name when subsetting column (Geneset1 --> Geneset1_Cluster)
RunCoxWrapper<-function(geneSet_groups, AllMeta, mySubset="AllSamples", 
                        controlFor=FALSE, 
                        exportName="",exportCsv=FALSE, subsetToHighLow=FALSE,
                        resultsFolder="", appendToNumericVariable="_Cluster"){
  AllCox<-list()
  for (i in geneSet_groups){
    AllMeta_Use<- AllMeta
    message(paste0("Running ", i))
    if(subsetToHighLow){
      ## for categorical
      if(class(AllMeta_Use[,i])%in% c("factor", "character")){
        useToSubset<- i
        HighLowClusters<- AllMeta_Use[,useToSubset]!="Mid"
        AllMeta_Use<- AllMeta_Use[HighLowClusters,]
        print(table(AllMeta_Use[,useToSubset]))
      } 
      ## For numeric
      if(class(AllMeta_Use[,i]) %in% c("numeric", "integer")){
        useToSubset<- paste0(i, appendToNumericVariable)
        HighLowClusters<- AllMeta_Use[,useToSubset]!="Mid"
        AllMeta_Use<- AllMeta_Use[HighLowClusters,]
        print(summary(AllMeta_Use[,useToSubset]))
      }

    }else{AllMeta_Use<- AllMeta}
    AllCox[[i]]<-RunCox(i, AllMeta_Use, controlFor = controlFor, mySubset=mySubset)
    
    message("")
  }
  
  cox_results<- WriteOutCoxResults(AllCox)
  if (exportCsv){
    write.csv(cox_results, paste0(resultsFolder, exportName, "_", mySubset, "_results.csv"))
  }
  
  names(AllCox)<- geneSet_groups
  return(list(cox_objs=AllCox, cox_df=cox_results))
}

########################
## PREPROCESSING TCGA ## 
########################

## Used to help convert GDC IDs (default for TCGA data) to TCGA IDs
## To match count matrices to survival data, need TGCA ID
## But the count matrices have GDC IDs
## Using the manifest as input, this function creates a "Payload" file that can then be used (see message) to 
## generate a txt file matching your GDC IDs
## This can then be used to convert your GDC IDs to TGCA IDs
GetPayloadForTGCAidConversion<- function(manifest_dir, cohortName){
  manifest_file<- paste0(manifest_dir, "gdc_manifest_", cohortName, ".txt")
  output_file<- paste0(manifest_dir, "Payload_", cohortName, ".txt")
  x=read.table(manifest_file,header = T)
  manifest_length= nrow(x)
  id= toString(sprintf('"%s"', x$id))
  Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '
  Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id","size":'
  Part3= paste(shQuote(manifest_length),"}",sep="")
  Sentence= paste(Part1,id,Part2,Part3, collapse=" ")
  write.table(Sentence,output_file,quote=F,col.names=F,row.names=F)
  message(paste0("curl --request POST --header \"Content-Type: application/json\" --data @Payload_",
                 cohortName,
                 ".txt \"https://api.gdc.cancer.gov/files\" > File_metadata_",
                 cohortName,
                 ".txt"))
}

## After creating payload file --> metadata file, use metadata file to convert GDC ID --> TCGA ID
ConvertGDCtoTCGA<- function(manifest_dir, cohortName, fileExtension, cm){
  
  ## Read in file created for converting IDs (after executing CL script in message)
  ConvertIDs<- read.delim(paste0(manifest_dir, "File_metadata_", cohortName, ".txt"))
  ConvertIDs$MyFileID<- gsub(fileExtension, "", ConvertIDs$file_name)
  ConvertIDs<- ConvertIDs[,colnames(ConvertIDs) %in% c("MyFileID",
                                                       "cases.0.submitter_id",
                                                       "cases.0.samples.0.sample_type")]
  colnames(ConvertIDs)<- c("SampleType", "TCGAID", "FileName")
  
  ## Remove recurrent tumors
  ConvertIDs<- ConvertIDs[ConvertIDs$SampleType=="Primary Tumor",]
  
  ## Remove recurrent tumors from count matrix, order rows/columns the same way, then rename columns
  ConvertIDs<- ConvertIDs[ConvertIDs$FileName %in% colnames(cm),]
  cm<- cm[,colnames(cm) %in% ConvertIDs$FileName]
  cm<- cm[,order(colnames(cm))]
  ConvertIDs<- ConvertIDs[order(ConvertIDs$FileName),]
  sum(colnames(cm) == ConvertIDs$FileName)
  colnames(cm)<- ConvertIDs$TCGAID
  cm<-cm[,order(colnames(cm))]
  return(cm)
  
}

##########
## MISC ##
##########

## Helper function for reading in all files for a cohort, merge together
## Input: path = path to data
##        pattern = pattern for files
ReadAndMerge<-function(path, pattern=""){
  files = list.files(path = path, pattern=pattern)
  print(head(files))
  cm = read.table(paste0(path, files[1]))
  colnames(cm)<- c("gene", gsub(pattern, "", files[1]))
  print(head(cm))
  
  for (f in 2:length(files)){
    tmp = read.table(paste0(path, files[f]))
    colnames(tmp)<- c("gene", gsub(pattern, "", files[f]))
    cm = merge(cm, tmp, by="gene")
  }
  
  rownames(cm)<- cm$gene; cm<-cm[,-1]
  return(cm)
}

## Convert ensembl IDs to gene symbols
## Note: no longer use this function- use at own risk
EnsemblToGeneSymbol<- function(cm_input){
  cm<-cm_input
  rownames(cm)<- gsub("\\..*", "", rownames(cm))
  geneSymbol<-na.omit(AnnotationDbi::select(org.Hs.eg.db, keys=rownames(cm), 
                                            columns=c("SYMBOL"), keytype="ENSEMBL"))
  
  geneSymbol<- geneSymbol[!(duplicated(geneSymbol$ENSEMBL)),]
  geneSymbol<- geneSymbol[!(duplicated(geneSymbol$SYMBOL)),]
  rownames(geneSymbol)<-geneSymbol$ENSEMBL
  
  ## Merge back with count matrix
  cm<- merge(cm, geneSymbol, by=0)
  rownames(cm)<- cm$SYMBOL
  cm<-cm[,!(colnames(cm) %in% c("SYMBOL", "ENSEMBL", "Row.names"))]
  return(cm)
}


## Convert raw counts to TPM using gene lengths
## Note: no longer use this function- use at own risk
CountsToTPM<- function(cm, id="ensGene"){
  ## Get lengths of genes, then remove any genes without this information (mostly non-coding)
  rownames(cm)<- gsub("\\..*", "", rownames(cm))
  geneLengths<-getlength(rownames(cm), "hg19", id)
  cm<-cm[!is.na(geneLengths),]
  
  ## Convert to TPM
  tmp <- cm/na.omit(geneLengths)
  cm<-(t(t(tmp)*1e6/colSums(tmp)))
  
  return(cm)
}


##########################################################
## HELPER FUNCTIONS HIGHLY SPECIFIC TO SYNAPTIC PROJECT ##
##########################################################

## saved here just in case, but unlikely to be useful for other projects

## Read in metadata, remove any samples that don't match between meta/tpm/fpkm. Switch meta censoring.
MatchMetaCounts<-function(fpkm, tpm, name){
  ## Get list of samples common to 3 dfs
  survival_meta<-read.delim(paste0("Metadata/RawSurvival_Clinical_Metadata/survival_", name, ".tsv"))
  samples<-survival_meta[survival_meta$submitter_id %in% colnames(fpkm), "submitter_id"]
  samples<-samples[samples %in% colnames(tpm)]
  
  
  ## Subset each df to common samples, order same way
  survival_meta<-survival_meta[survival_meta$submitter_id %in% samples,]
  survival_meta<- survival_meta[order(survival_meta$submitter_id),]
  
  fpkm<-fpkm[,colnames(fpkm) %in% samples]
  fpkm<-fpkm[,order(colnames(fpkm))]
  
  tpm<-tpm[,colnames(tpm) %in% samples]
  tpm<-tpm[,order(colnames(tpm))]
  
  ## Double check ordering
  print(sum(survival_meta$submitter_id == colnames(fpkm)))
  print(sum(colnames(fpkm)==colnames(tpm)))
  
  ## Switch metadata censoring data
  survival_meta$censored<- !(as.logical(survival_meta$censored))
  
  return(list(meta=survival_meta, fpkm=fpkm, tpm=tpm))
}

## Function to norm/center/mean counts
## If cm is already centered, skip center, make tpm slot NULL
centerCounts_meanCounts<- function(cm, alreadyCentered=FALSE){
  cm_list<- list()
  if(alreadyCentered){
    cm_list$tpm<- NA
    cm_list$cm_center<- cm
    cm_list$cm_mean<- rowMeans(cm_list$cm_center)
  } else{
    normCenter<- NormCenter(cm)
    cm_list$tpm<- cm
    cm_list$cm_center<- normCenter$center_data
    cm_list$cm_mean<- rowMeans(log2(cm_list$tpm + 1))
  }
  return(cm_list)
}

Check_SampleIDsMatch<- function(cm, survival, clinical){
  print(paste0("Samples in cm: ", ncol(cm)))
  print(paste0("Samples in survival: ", nrow(survival)))
  
  print(paste0("Number samples match between cm/survival: ",
               sum(colnames(cm)==survival$submitter_id)))
  
}

Subset_AndReCenter_ReMean<- function(cohort, samplesUse, tpmAvailable){
  cohort$survival<- cohort$survival[cohort$survival$submitter_id %in% samplesUse,]
  cohort$clinical<- cohort$clinical[cohort$clinical$submitter_id %in% samplesUse,]
  print(paste0("Samples in clinical and survival: ", nrow(cohort$survival), ", ", nrow(cohort$clinical)))
  
  
  ## If TPM is availale: re center
  if(tpmAvailable){
    cohort$tpm<- cohort$tpm[,colnames(cohort$tpm) %in% samplesUse]
    print(paste0("Samples in TPM: ", ncol(cohort$tpm)))
    cm_list<- centerCounts_meanCounts(cohort$tpm, alreadyCentered = FALSE)
    cohort$cm_center<- cm_list$cm_center
    cohort$cm_mean<- cm_list$cm_mean
  } else{
    cohort$cm_center<- cohort$cm_center[, colnames(cohort$cm_center) %in% samplesUse]
    print(paste0("Samples in centered counts: ", ncol(cohort$cm_center)))
    cm_list<- centerCounts_meanCounts(cohort$cm_center, alreadyCentered = TRUE)
    cohort$cm_mean<- cm_list$cm_mean
  }
  return(cohort)
}

## From a dataframe of scores for each geneset, split samples into high/low expressors
## For each module score, split samples into "high" or "low" expression
## Input: scores = a dataframe with genesets as columns, samples as rows 
##        splitBy = denotes how to split samples into high/low geneset scorers
SplitHighLow<- function(scores, splitBy="Median"){
  AllClusters<-data.frame(submitter_id=Scores$submitter_id)
  for (m in names(marker_list_InCohort)){
    print(m)
    df<- cbind(Scores[,"submitter_id"], Scores[,m]); df<-as.data.frame(df)
    ClusterName<- paste0(m, "_Cluster")
    df$V2<-as.numeric(df$V2)
    
    ## Alternative methods for splitting into high/low
    ## Split into top25%/mid50%/bottom25%
    if(SplitBasedOn=="SplitInto3Quartiles"){
      df[,"ClusterName"]<- df$V2>quantile(df$V2)["75%"]
      high<- df[df$V2>quantile(df$V2)["75%"],]; high[,"ClusterName"]<-"High"
      low<- df[df$V2<quantile(df$V2)["25%"],]; low[,"ClusterName"]<-"Low"
      mid<-df[!(df$V1 %in% c(high$V1, low$V1)),]; mid[,"ClusterName"]<-"Mid"
      df<-rbind(high,low); df<-rbind(df, mid)
    } 
    ## Split based on median
    if(SplitBasedOn=="Median"){
      df[,"ClusterName"]<- df$V2>median(df$V2)
      df[,"ClusterName"]<- gsub("FALSE", "Low", 
                                gsub("TRUE", "High", df[,"ClusterName"]))
    } 
    
    ## Convert to factor, with "Low" as comparison
    #df$ClusterName<- factor(df$ClusterName, levels=c("High", "Low"))
    colnames(df)<- c("submitter_id", m, ClusterName)
    AllClusters<-merge(AllClusters,df, by="submitter_id")
  }
}

## Run and plot survival analysis
## Input: meta = metadata with time to death, patient, vital status, and gene set group (high/low)
##        GeneSetName = name of geneset in meta
##        LowHighOnly = when splitting 25/50/25, use only 25/25
##        plotName = survival plot name- if NULL, will just use gene set name
RunSurvivalAnalysis<-function(meta, GeneSetName, LowHighOnly=FALSE, plotName=NULL){
  ## Prep data for survival analysis- need time to death, patient, vital status, and group (high/low)
  survivalInput<- meta
  survivalInput$censored<- as.logical(survivalInput$censored)
  survivalInput$time<- as.numeric(survivalInput$time)
  survivalInput$Marker<-survivalInput[,GeneSetName]
  if (LowHighOnly){
    survivalInput<-survivalInput[!(survivalInput$Marker=="Mid"),]
  }
  
  ## Run survival analysis
  fit <- survfit(Surv(time, censored) ~Marker,
                 data = survivalInput)
  
  ## Create dataframe of survival results
  pvalue<-surv_pvalue(fit, survivalInput)$pval.txt
  pvalue<-gsub("p = ", "", pvalue)
  myFit_median<-surv_median(fit)
  highMarker_median<-myFit_median[myFit_median$strata=="Marker=High", "median"]
  lowMarker_median<-myFit_median[myFit_median$strata=="Marker=Low", "median"]
  
  ## If samples have been split by quartile AND LowHighOnly=FALSE:
  ## Extract median for mid as well as high/low
  ## If not- just use median survival for high/low
  if("Marker=Mid" %in% myFit_median$strata & LowHighOnly==FALSE){
    midMarker_median<- myFit_median[myFit_median$strata=="Marker=Mid", "median"]
    df<-data.frame(Marker=GeneSetName, pvalue=pvalue,
                   HighMarkerMedianSurvival=highMarker_median,
                   MidMarkerMedianSurvival=midMarker_median,
                   LowMarkerMedianSurvival=lowMarker_median)
  } else {
    df<-data.frame(Marker=GeneSetName, pvalue=pvalue,
                   HighMarkerMedianSurvival=highMarker_median,
                   LowMarkerMedianSurvival=lowMarker_median)
  }
  
  ## Visualize with survminer
  plotTitle<- ifelse(is.null(plotName), GeneSetName, plotName )
  p<-ggsurvplot(fit, pval=TRUE, survivalInput)+ggtitle(plotTitle)
  return(list(plot=p, fit_df=df))
}