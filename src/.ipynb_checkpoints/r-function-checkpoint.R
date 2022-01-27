# library load
suppressPackageStartupMessages({
    library(survival)
    library(ranger)
    library(NbClust)
    library(survminer)
    library(tidyverse)
    library(TCGAbiolinks)
    library(SummarizedExperiment)
    library(DESeq2)
})

# KM - plot
km_survival <- function(df, file_name, path){
    fit <- survfit(Surv(time = OS.time, event = OS) ~ group, data = df)
    ggsurvplot(fit, data = df, pval = T, conf.int = T, title = file_name)
    ggsave(paste0(path, file_name, "_AE_km_plot.png"), )
}    

# log-rank test
log_rank_test <- function(df){   
  # column name extraction
  df_name <- colnames(df)
  df_name_filter <- df_name[str_detect(string = df_name, pattern = "Feature")]
  
  # log_rank test
  df_log_rank <- lapply(X = df_name_filter, function(col_name){
    cox <- coxph(
      formula = as.formula(paste0("Surv(time = OS.time, event = OS) ~", col_name)), 
      data = df)
    cox_result <- summary(cox)
    tibble(Features = col_name, log_p_value = cox_result$logtest["pvalue"]) %>% return()
  }) %>% bind_rows()
  
  # p-value 0.05 cuttoff
  df_log_rank %>% 
    filter(log_p_value < 0.05) %>% 
    return()
}

# random_forest test
nb_cluster_test <- function(df){ 
  
  nc <- NbClust(df,min.nc=2,max.nc=9,method="kmeans", index = c("silhouette","cindex"))
  nc$All.index %>% as_tibble() %>% select("Silhouette") %>% return()
}

# function
run_edgeR <- function(pr_name, sample_group_path){
    
  suppressMessages({
      sample_group <- read_delim(file = sample_group_path, delim = "\t", show_col_types = FALSE)
      query <- GDCquery(project = paste0("TCGA-", pr_name), 
                        data.category = "Gene expression",
                        data.type = "Gene expression quantification",
                        experimental.strategy = "RNA-Seq",
                        platform = "Illumina HiSeq",
                        file.type = "results",
                        barcode = sample_group %>% pull(1), 
                        legacy = TRUE)

      GDCdownload(query)
      RnaseqSE <- GDCprepare(query)
      Rnaseq_CorOutliers <- assay(RnaseqSE)

      # subgroup names
      sub0 <- sample_group %>% filter(group == 0) %>% pull(1)
      sub1 <- sample_group %>% filter(group == 1) %>% pull(1)

      sample_subgroup <- Rnaseq_CorOutliers %>% colnames() %>% lapply(X = ., FUN = function(value){

        value_trans <- str_extract_all(value, pattern = "TCGA-[:alnum:]+-[:alnum:]+-[:digit:]+") %>%  unlist()
        subgroup_df <- tibble(sample_barcode = NA, subgroup = NA, .rows = 0)


        if(value_trans %in% sub0){
          subgroup_df <- tibble(sample_barcode = value, subgroup = 0)
        } else {
          subgroup_df <- tibble(sample_barcode = value, subgroup = 1)
        }

        return(subgroup_df)

      }) %>% bind_rows()

      sub0_name <- sample_subgroup %>% filter(subgroup == 0) %>% pull(sample_barcode)
      sub1_name <- sample_subgroup %>% filter(subgroup == 1) %>% pull(sample_barcode)

      # normalization of genes
      dataNorm <- TCGAanalyze_Normalization(tabDF = Rnaseq_CorOutliers, geneInfo =  geneInfo)

      # quantile filter of genes
      dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                        method = "quantile", 
                                        qnt.cut =  0.25)

      #Diff.expr.analysis (DEA)
      dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,sub0_name],
                                  mat2 = dataFilt[,sub1_name],
                                  Cond1type = "Sub-0",
                                  Cond2type = "Sub-1",
    #                               fdr.cut = 0.01 ,
    #                               logFC.cut = 1,
                                  method = "glmLRT",
                                  paired = F)

      # DEGs table with expression values in normal and tumor samples
      dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Sub-0","Sub-1",
                                                dataFilt[,sub0_name],dataFilt[,sub1_name])
#   write_delim(dataDEGsFiltLevel, file = paste0(dea_result_path, "_", pr_name, "_EDGER_", file_name, ".txt"), delim = "\t")   
  })  
  
  return(dataDEGsFiltLevel)
}
run_edgeR_pancan <- function(sample_group_path){
  sample_group <- read_delim(file = sample_group_path, delim = "\t", show_col_types = FALSE)    
  study_list <- tcga_available()$Study_Abbreviation[-34] %>% 
    lapply(X = ., FUN = function(value){
      paste0("TCGA-", value) %>% return()
    }) %>% as.character()
  
  suppressMessages({
    tcga_mrna <- mclapply(X = study_list, FUN = function(type){
    
    error_occur <- FALSE
    tryCatch(
      expr = {
        query <- GDCquery(project = type, 
                          data.category = "Gene expression",
                          data.type = "Gene expression quantification",
                          experimental.strategy = "RNA-Seq",
                          platform = "Illumina HiSeq",
                          file.type = "results",
                          barcode = sample_group %>% pull(1), 
                          legacy = TRUE)
        # Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
        GDCdownload(query, directory = dir)
        
        # rsem.genes.results as values
        RnaseqSE <- GDCprepare(query, directory = dir)
        
      },
      error = function(e) { 
        error_occur <<- TRUE
      }
    )
    
    if(error_occur){
      return(NULL)
    } else {
      RnaseqSE %>%
        assay() %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "gene") %>% return()
    }
    
  }, mc.cores = 20)
  
  # inner join
  reduce_join <- function(left, right){
    if(is.null(right)){
      return(left)
    } else{
      inner_join(left, right, by = "gene") %>% return()
    }
  }
  tcga_mrna_join <- purrr::reduce(tcga_mrna, reduce_join) %>% 
    column_to_rownames(var = "gene") %>% as.matrix()
  
  # normalization of genes # ~ 10 min
  dataNorm <- TCGAanalyze_Normalization(tabDF = tcga_mrna_join, geneInfo =  geneInfo)
  
  # subgroup names
  sub0 <- sample_group %>% filter(km_tsne1 == 0) %>% pull(1)
  sub1 <- sample_group %>% filter(km_tsne1 == 1) %>% pull(1)
  
  sample_subgroup <- tcga_mrna_join %>% colnames() %>% lapply(X = ., FUN = function(value){
    
    value_trans <- str_extract_all(value, pattern = "TCGA-[:alnum:]+-[:alnum:]+-[:digit:]+") %>%  unlist()
    subgroup_df <- tibble(sample_barcode = NA, subgroup = NA, .rows = 0)
    
    
    if(value_trans %in% sub0){
      subgroup_df <- tibble(sample_barcode = value, subgroup = 0)
    } else {
      subgroup_df <- tibble(sample_barcode = value, subgroup = 1)
    }
    
    return(subgroup_df)
    
  }) %>% bind_rows()
  
  sub0_name <- sample_subgroup %>% filter(subgroup == 0) %>% pull(sample_barcode)
  sub1_name <- sample_subgroup %>% filter(subgroup == 1) %>% pull(sample_barcode)
  
  # quantile filter of genes
  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                    method = "quantile", 
                                    qnt.cut =  0.25)
  
  #Diff.expr.analysis (DEA) # ~ 20min
  dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,sub0_name],
                              mat2 = dataFilt[,sub1_name],
                              Cond1type = "Sub-0",
                              Cond2type = "Sub-1",
                              fdr.cut = 0.01 ,
                              logFC.cut = 1,
                              method = "glmLRT",
                              paired = F)
  
  # DEGs table with expression values in normal and tumor samples
  dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Sub-0","Sub-1",
                                            dataFilt[,sub0_name],dataFilt[,sub1_name])    
  })  
  

  
  return(dataDEGsFiltLevel)
}
run_deseq <- function(pr_name, sample_group_path){
    
  suppressMessages({
      sample_group <- read_delim(file = sample_group_path, delim = "\t", show_col_types = FALSE)
      query <- GDCquery(project = paste0("TCGA-", pr_name), 
                        data.category = "Gene expression",
                        data.type = "Gene expression quantification",
                        experimental.strategy = "RNA-Seq",
                        platform = "Illumina HiSeq",
                        file.type = "results",
                        barcode = sample_group %>% pull(1), 
                        legacy = TRUE)

      GDCdownload(query)
      RnaseqSE <- GDCprepare(query)
      Rnaseq_CorOutliers <- assay(RnaseqSE) # to matrix

      # normalization of genes, # quantile filter of genes
      dataNorm <- TCGAanalyze_Normalization(tabDF = Rnaseq_CorOutliers, geneInfo =  geneInfo)
      dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                        method = "quantile", 
                                        qnt.cut =  0.25)

      # metadata
      metadata <- tibble(sample = RnaseqSE$sample_submitter_id, barcode = RnaseqSE$barcode) %>% 
        mutate(sample = str_extract_all(sample, pattern = "TCGA-[:alnum:]+-[:alnum:]+-[:digit:]+") %>% unlist()) %>% 
        left_join(x = ., y = sample_group, by = "sample") %>% 
        mutate(group = as.factor(ifelse(group == 0, "Sub0", "Sub1")))

      tcga_se <- DESeqDataSetFromMatrix(countData = dataFilt, colData = metadata, design = ~ group)
      tcga_deseq <- DESeq(tcga_se)

      tcga_deseq_result <- results(tcga_deseq, tidy = T)

    #   write_delim(tcga_deseq_result, file = paste0(dea_result_path, "_", pr_name, "_DESEQ2_", file_name, ".txt"), delim = "\t")
        
  })  

  return(tcga_deseq_result)
}