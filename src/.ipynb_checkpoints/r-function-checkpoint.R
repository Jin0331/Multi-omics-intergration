library(tidyverse)
library(survival)

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