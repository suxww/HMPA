if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")# 
#BiocManager::install("TCGAbiolinks", force = TRUE) 
rm(list = ls())

library(TCGAbiolinks)
library(tidyverse)

setwd("~/cptac_reanalysis/UCEC/results/03r_exp")
library("rjson")
json <- jsonlite::fromJSON("metadata.cart.2024-02-19.json")
View(json)
id <- json$associated_entities[[1]][,1]
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
case_id <- sapply(json$associated_entities,function(x){x[,3]})
file_sample <- data.frame(case_id,file_name=json$file_name,sample_id)  

#match ID
clinical <- read_delim("clinical.cart.2024-02-19/clinical.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
#clinical <- read_delim("~/Desktop/cptac_reanalysis/LUAD/source/clinical.cart.2023-12-17/clinical.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
file_sample <- data.frame(case_id,file_name=json$file_name,sample_id)  
file_sample<-merge(file_sample,clinical,by = "case_id")

file_sample <- data.frame(sample_id,file_name=json$file_name) 


write.csv(file_sample,'file_sample.csv',row.names = F)

#count_file <- list.files('gdc_download_20231217_122821.032513',pattern = '*.tsv',recursive = TRUE)
count_file <- list.files('gdc_download_20240218_225657.919585',pattern = '*gene_counts.tsv',recursive = TRUE)
count_file_name <- strsplit(count_file,split='/')
count_file_name <- sapply(count_file_name,function(x){x[2]})

matrix = data.frame(matrix(nrow=60660,ncol=0))
for (i in 1:length(count_file_name)){
  path = paste0('gdc_download_20240218_225657.919585/',count_file[i])
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  data <- data[7] 
  colnames(data) <- file_sample$case_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}

write.csv(matrix,'case_COUNT_matrix.csv',row.names = TRUE)

