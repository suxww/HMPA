rm(list=ls())
library(ggplot2)
library(ggsci)
library(gridExtra)
library(tidyr)
library(ggpubr)
library(survival)
library(survminer)
library(read_excel)
survival_data <- read_excel("~/Desktop/DB_files/tools_clinical/clinical_ccrcc.xlsx")
gene_columns <- c(colnames(survival_data[,10:797]))
all_p_values <- c()

# transfer number
survival_data_exp <- sapply(survival_data[, 10:797], as.numeric)
survival_data <- cbind(survival_data_exp,survival_data[,1:9])
all_time <- list()
all_surv <- list()

for (gene_column in gene_columns) {

  # select column
  s_fig_data <- survival_data[, c("time", "event", gene_column)]
  
  # split expression
  s_fig_data$hmpa_expression_group <- ifelse(s_fig_data[, gene_column] > median(s_fig_data[, gene_column]), 'high', 'low')
  
  # Kaplan-Meier
  km_fit <- survfit(Surv(time, event) ~ hmpa_expression_group, data = s_fig_data)
  
  test_result <- survdiff(Surv(time, event) ~ hmpa_expression_group, data = s_fig_data)
  # p value
  p_value <- test_result$pvalue
  all_p_values <- c(all_p_values, p_value)
  
  all_time[[gene_column]] <- km_fit[["time"]]
  all_surv[[gene_column]] <- km_fit[["surv"]]
}

result_data <- data.frame(
  gene_column = character(),
  p_value = numeric(),
  time_data = numeric(),
  surv_data = numeric(),
  group = character()
)

for (i in 1:length(gene_columns)) {
  result_data <- rbind(result_data, data.frame(
    gene_column = gene_columns[i],
    p_value = all_p_values[i],
    time_data = unlist(all_time[[gene_columns[i]]]),
    surv_data = unlist(all_surv[[gene_columns[i]]])
  ))
}

write.csv(result_data,"result_data.csv",row.names = F)

# group
grouped_result <- result_data %>%
  group_by(gene_column, group) %>%
  summarize(
    time_surv_data = paste(time_data, surv_data, sep = ", ", collapse = "; ")
  ) %>%
  ungroup()


final_result <- data.frame(
  gene_column = grouped_result$gene_column,
  p_value = rep(unique(result_data$p_value), each = n_distinct(grouped_result$group)),
  group = grouped_result$group,
  surv = grouped_result$time_surv_data,
  stringsAsFactors = FALSE
)

write.csv(final_result,"result_data.csv",row.names = F)
