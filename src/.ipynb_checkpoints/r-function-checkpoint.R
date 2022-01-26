# install.package
requiredPackages <- c('survival','ranger', 'NbClust', 'survminer', 'tidyverse', 'BiocManager')
for(p in requiredPackages){
    if(!require(p,character.only = TRUE)){install.packages(p)}}

for(p in c("TCGAbiolinks", "SummarizedExperiment")){
    if(!require(p,character.only = TRUE)){BiocManager::install(p)}}

# library load
suppressPackageStartupMessages({
    library(survival)
    library(ranger)
    library(NbClust)
    library(survminer)
    library(tidyverse)
    library(TCGAbiolinks)
    library(SummarizedExperiment)
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

run_deg <- function(project, sample){
  # Query platform Illumina HiSeq with a list of barcode 
  query <- GDCquery(project = paste0("TCGA-", pr_name), 
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    experimental.strategy = "RNA-Seq",
                    platform = "Illumina HiSeq",
                    file.type = "results",
                    barcode = sample %>% pull(1), 
                    legacy = TRUE)
  
  # Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
  GDCdownload(query)
  
  # Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
  # rsem.genes.results as values
  LIHCRnaseqSE <- GDCprepare(query)
  
  # For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
  LIHCRnaseq_CorOutliers <- TCGAanalyze_Preprocessing(LIHCRnaseqSE)
  
  # subgroup names
  sub0 <- sample %>% filter(group == 0) %>% pull(1)
  sub1 <- sample %>% filter(group == 1) %>% pull(1)
  
  sample_subgroup <- LIHCRnaseq_CorOutliers %>% colnames() %>% lapply(X = ., FUN = function(value){
    
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
  dataNorm <- TCGAanalyze_Normalization(tabDF = LIHCRnaseq_CorOutliers, geneInfo =  geneInfo)
  
  # quantile filter of genes
  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                    method = "quantile", 
                                    qnt.cut =  0.25)
  
  #Diff.expr.analysis (DEA)
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
  
  return(dataDEGsFiltLevel)
}