
library(tidyverse)
library(DMwR2)


data1 = read.csv('Cohort1_Search.xlsx - 208-20221118.pg_matrix.csv', header = T)
data2 = read.csv('Cohort2_result.xlsx - Cohort2_result.csv', header = T)
data3 = read.csv('Cohort3_Search.xlsx - report.pg_matrix.csv', header = T)

# 数据清洗函数定义
data.cleaning = function(df){

  anno_cols = c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description")
  intensity.names = setdiff(names(df), anno_cols)
  

  df[intensity.names][df[intensity.names] == 0] <- NA
  
  LOG.names = sub("^", "LOG2.", intensity.names) 
  df[LOG.names] = log2(df[intensity.names])
  
  filter_valids = function(df, conditions, at_least_one = FALSE) {
    log2.names = grep("^LOG2", names(df), value = TRUE)  
    
    cond.filter = sapply(1, function(i) {
      df2 = df[log2.names]   
      df2 = as.matrix(df2)   
      sums = rowSums(is.na(df2)) 
      sums <= length(log2.names)*0.3   
    })
    
    if (at_least_one) {
      df$KEEP = apply(cond.filter, 1, any)
    } else {
      df$KEEP = apply(cond.filter, 1, all)
    }
    
    return(df)  
  }
  

  df.F = filter_valids(df, conditions = NULL, at_least_one = TRUE)
  
  impute_data = function(df, width = 0.3, downshift = 1) {
    LOG2.names = grep("^LOG2", names(df), value = TRUE)
    impute.names = sub("^LOG2", "impute", LOG2.names)
    
    df[impute.names] = lapply(LOG2.names, function(x) !is.finite(df[, x]))
    
    # 确保数据中没有 Inf/-Inf 干扰 KNN 插补
    for(col in LOG2.names){
      df[!is.finite(df[, col]), col] <- NA
    }
    
    set.seed(1)
    df[LOG2.names] = knnImputation(df[LOG2.names], k=3)
    return(df)
  }
  
  df.FNI = impute_data(df.F[df.F$KEEP,])
  return(df.FNI)
}


data1 = data.cleaning(data1)
data1 = data1[data1$KEEP, c('Protein.Ids', grep("^LOG2", names(data1), value = TRUE))]

data2 = data.cleaning(data2)
data2 = data2[data2$KEEP, c('Protein.Ids', grep("^LOG2", names(data2), value = TRUE))]

data3 = data.cleaning(data3)
data3 = data3[data3$KEEP, c('Protein.Ids', grep("^LOG2", names(data3), value = TRUE))]


data1_t = t(data1)
data2_t = t(data2)
data3_t = t(data3)

write.table(data1_t, 'Cohort1_processed.csv', col.names = FALSE, sep = ',')
write.table(data2_t, 'Cohort2_processed.csv', col.names = FALSE, sep = ',')
write.table(data3_t, 'Cohort3_processed.csv', col.names = FALSE, sep = ',')


data1_final = read.csv('Cohort1_processed.csv', header = T)
colnames(data1_final)[1] = 'label'
data1_final$label = gsub('^LOG2.','',data1_final$label)

data1_final$label = gsub('[0-9_\\.].*','',data1_final$label) 

data2_final = read.csv('Cohort2_processed.csv', header = T)
colnames(data2_final)[1] = 'label'
data2_final$label = gsub('^LOG2.','',data2_final$label)
data2_final$label = gsub('[0-9_\\.].*','',data2_final$label)

data3_final = read.csv('Cohort3_processed.csv', header = T)
colnames(data3_final)[1] = 'label'
data3_final$label = gsub('^LOG2.','',data3_final$label)
data3_final$label = gsub('[0-9_\\.].*','',data3_final$label)

data1_final$label = factor(data1_final$label)
data2_final$label = factor(data2_final$label)
data3_final$label = factor(data3_final$label)


common_cols = intersect(names(data1_final), intersect(names(data2_final), names(data3_final)))

data1_final = data1_final[, common_cols]
data2_final = data2_final[, common_cols]
data3_final = data3_final[, common_cols]

