suppressPackageStartupMessages({
    library(tidyverse)
    library(survival)
    library(ranger)
    library(NbClust)
    library(survival)
    library(survminer)
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