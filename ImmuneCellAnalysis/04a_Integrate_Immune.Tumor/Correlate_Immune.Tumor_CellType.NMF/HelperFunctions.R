

## Run lm on dependent variable of interest comparing to list of confounders
## input: df should be df with dependent values stored in dep_value_column
## potential confounders each stored in separate columns, 1 for each in confounder_columns
lm_checkForConfounders<- function(df, dep_value_column,  confounder_columns){
  all_test<- lapply(confounder_columns, function(x){
    message(paste0("Running lm for ", x))
    
    ## Variables that are originally numeric may be changed when removing NAs. save class here so you can reassign it
    original_class<- class(df[[x]])
    
    ## Remove any NAs from df
    df[[x]]<- gsub("NA", NA, df[[x]]) ## replace any character NAs with actual NAs
    na_values<- is.na(df[[x]])
    if(sum(na_values)>0){
      message(paste0(sum(na_values), " NA values removed"))
      df_tmp<- df[!na_values,] ## remove any NA values
    }else{
      df_tmp<- df
    }
    
    ## Reassign variable to numeric if needed
    if(original_class %in% c("numeric", "integer")){
      df_tmp[[x]]<- as.numeric(df_tmp[[x]])
    }
    
    ## Run lm for dependent vs potential confounder
    formula_use<- paste0(dep_value_column, " ~ ", x)
    print(paste0("Model used: ", formula_use))
    model<- summary(lm(formula_use, df_tmp))
    
    ## Extract p values. For categorical confounders, extract lowest pvalue
    pvalues<- as.data.frame(coef(model))
    pvalues<- pvalues[-1, "Pr(>|t|)"] ## first row is intercept- ignore
    return(list(model=model, min_pvalue=min(pvalues)))
    
  })
  names(all_test)<- confounder_columns
  return(all_test)
}


## Run linear regression for multiple programs (separately for each). Can also input programs to control for if desired
## Input: df with 2 columns for dependent: actual values (dep_value_column) + column to split by (usually programs; dep_splitter_column)
## same for independent. ControlFor: vector of 0-n variables to control for
linearRegression_multiplePrograms<- function(df, dep_value_column,ind_value_column,
                                             dep_splitter_column,ind_splitter_column,
                                             controlFor){
  ## Identify all indp/dep programs to use
  ind_programs<- as.character(unique(df[[ind_splitter_column]]))
  dep_programs<- as.character(unique(df[[dep_splitter_column]]))
  
  message(paste0("Running lm for: ", paste(ind_programs, collapse = ", "), " vs ",
                 paste(dep_programs, collapse = ", ")))
  message(paste0("Controlling for ", paste(controlFor, collapse = ", ")))
  
  ## Generate list of all dep/indp comparisons possible
  all_pairs<- expand.grid(ind_programs, dep_programs)
  colnames(all_pairs)<- c("ind_program", "dep_program")
  all_pairs$ind_program<- as.character(all_pairs$ind_program); all_pairs$dep_program<- as.character(all_pairs$dep_program)
  all_pairs<- lapply(1:nrow(all_pairs), function(x){all_pairs[x,]})
  
  ## Create model input, adding in variables to control for if needed
  formula_use_char<- paste0(dep_value_column, " ~ ", ind_value_column)
  if(length(controlFor)!= 0){
    formula_use_char<- paste0(formula_use_char, " + ", paste(controlFor, collapse = " + "))
  }
  print(paste0("Formula used: ", formula_use_char))
  formula_use<- as.formula(formula_use_char)
  
  ## Cycle through indp programs, running lm on each dep program
  all_lm<- lapply(all_pairs, function(x){
    df_tmp<- df[df[[ind_splitter_column]]==x$ind_program &
                  df[[dep_splitter_column]]==x$dep_program,]
    
    ## Run lm
    model<-lm(formula_use, df_tmp)
    model_summary<- summary(model)
    model_res<- data.frame(pvalue=model_summary$coefficients[ind_value_column, "Pr(>|t|)"],
                           r2=model_summary$r.squared, 
                           call=formula_use_char,
                           ind_variable=x$ind_program,
                           dep_variable=x$dep_program)
    return(list(model=model, model_res=model_res))
  })
  model_res<- do.call("rbind", lapply(all_lm, function(x){x$model_res}))
  models<- lapply(all_lm, function(x)x$model)
  return(list(model_res=model_res, model=models))
}
## Plot dotplot of % of ind variable vs % of dep variable
## Include lm res on plot if desired
## ind/dep_variable_plot_name just adds the desired name to plot labels
## program_variable_name denotes the column name within df that contains programs
## xAxis_lastPlotOnly: if TRUE, only add x axis label to 1 plot 
plot_Proportions_withLm.res<- function(df, programs, pvalue.x=34, pvalue.y=75, r.x=32, r.y=65, 
                                       colors_use, program_variable_name="Program", lm_res=FALSE,
                                       independent_variable, dependent_variable,
                                       ind_variable_plot_name="", dep_variable_plot_name="",
                                       xAxis_lastPlotOnly=TRUE){
  df$program_variable_name<-df[[program_variable_name]]
  all_plots<- lapply(programs, function(x){
    ## Subset to single program
    tmp<- df[df[[program_variable_name]]==x,]
    
    ## Add on plotting names
    tmp$color<- tmp[[program_variable_name]]
    tmp$dep_var<- tmp[[dependent_variable]]
    tmp$ind_var<- tmp[[independent_variable]]
    
    
    p<-ggplot(tmp, aes(x=ind_var, y=dep_var, color=color))+
      geom_point()+
      theme_classic()+
      geom_smooth(se=FALSE, method="lm")+
      scale_color_manual(values=colors_use)+
      theme(legend.position = "none",
            axis.text = element_text(color="black", face="bold"),
            axis.title = element_text(color="black", face="bold"))+
      xlab(paste0("Percentage ", ind_variable_plot_name, " cells"))+
      ylab(paste0("Percentage ", dep_variable_plot_name, " cells"))+
      facet_grid(rows=vars(program_variable_name))
    
    ## Add pvalue/r2 to plot if available
    if(class(lm_res) != "logical"){
      ## Extract lm results for this cell type
      model_res<- lm_res[lm_res$ind_variable==x,]
      
      ## Add to plot
      p<- p+
        ggtitle(paste0("p value: ",  round(model_res$pvalue, 4), ", R2: ", round(model_res$r2, 4)))
        #annotate("text",x=pvalue.x, y=pvalue.y,label=paste0("p value: ",  round(model_res$pvalue, 4)))+
        #annotate("text",x=r.x, y=r.y,label=paste0("R2: ", round(model_res$r2, 4) ))
    }
    
    ## Optionally: Only keep xlab for last program
    if(xAxis_lastPlotOnly){
      if(x!=programs[length(programs)]){
        p<-p+xlab("")
      }
    }
    p
  })
  return(all_plots)
}

## Correlate proportions of 2 cell types with one another
## input: df that contains, for both ind and dep variables:
## dep/indp_value_column: actual values to be correlated
## dep/ind_splitter_column: column with programs/cell types to be split by. All pairwise combinations will be run
correlateProportions<- function(df, dep_value_column,ind_value_column,
                                dep_splitter_column,ind_splitter_column){
  ## program lists for dependent/independent
  ind_programs<- as.character(unique(df[[ind_splitter_column]]))
  dep_programs<- as.character(unique(df[[dep_splitter_column]]))
  
  message(paste0("Running pearson correlation for: ", paste(ind_programs, collapse = ", "), " vs ",
                 paste(dep_programs, collapse = ", ")))
  
  ## Generate list of all dep/indp comparisons possible
  all_pairs<- expand.grid(ind_programs, dep_programs)
  colnames(all_pairs)<- c("ind_program", "dep_program")
  all_pairs$ind_program<- as.character(all_pairs$ind_program); all_pairs$dep_program<- as.character(all_pairs$dep_program)
  all_pairs<- lapply(1:nrow(all_pairs), function(x){all_pairs[x,]})
  
  ## Cycle through all these pairs, doing pearson correlation for their scores
  all_cor<- lapply(all_pairs,function(x){
    df_tmp<- df[df[[ind_splitter_column]]==x$ind_program &
                  df[[dep_splitter_column]]==x$dep_program,]
    cor<- cor(df_tmp[[ind_value_column]], df_tmp[[dep_value_column]])
    cor_df<- data.frame(PearsonCorrelation=cor, ind_variable=x$ind_program, dep_variable=x$dep_program)
    return(cor_df)
  })
  cor_df<- do.call("rbind", all_cor)
  return(cor_df)
}


## Function to calculate proportion of given program by sample within a seurat object
calculate_prop<- function(seurat_obj, sampleColumn="sampleid", programColumn){
  ## Add on columns to seurat object as needed
  seurat_obj$sampleid<- seurat_obj@meta.data[[sampleColumn]]
  seurat_obj$program<- seurat_obj@meta.data[[programColumn]]
  
  ## Remove samples with <10 total cells in programColumn
  samples_use<- table(seurat_obj$sampleid)
  samples_use<- samples_use[samples_use>=10]
  seurat_obj<- subset(seurat_obj,sampleid %in% names(samples_use))
  
  ## Calculate number of cells per sample per program
  prop_df<- as.data.frame(table(seurat_obj$sampleid, seurat_obj$program))
  colnames(prop_df)<- c("Sample", "Program", "Ncells")
  
  ## Calculate proportion of cells
  prop_df<- prop_df %>% group_by(Sample) %>% 
    mutate(perCells=Ncells/sum(Ncells) * 100) %>% as.data.frame()
  return(prop_df)
}