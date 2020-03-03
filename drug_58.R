rm(list = ls())
setwd("/cloud/project")
load('/cloud/project/environment/newdata2.RData')
source("/cloud/project/code/function_drug_58.R")

library(ggplot2)
library(reshape2)
library(ggpubr)
library(Hotelling)
library(knitr)
library(dendextend)
library(circlize)
library(VennDiagram)
library(HDtest)

# range of y-axis
plot.range = c(-0.16, 0.32) 
# significance levels 
# for Q344X
sig.lv1 = 0.05
# for Rho
sig.lv2 = 0.05
# folder 1, 2 and old, including normalized data
folder1 = list()
folder2 = list()
folder.old = list()

# 3. Light-Off 30 secs

# 3.1 Time slot
t_interval = 30
range_time = c(-t_interval:(t_interval-1))

# 3.2 Normalization
workingData = subset(current1.lightoff, current1.lightoff$time >=  - t_interval & current1.lightoff$time < t_interval)
folder1$off30data = normal.diy(workingData = workingData, baseline = baseline1, 
                               current.lightoff = current1.lightoff, current.lighton = current1.lighton)
workingData = subset(current2.lightoff, current2.lightoff$time >=  - t_interval & current2.lightoff$time < t_interval)
folder2$off30data = normal.diy(workingData = workingData, baseline = baseline2, 
                               current.lightoff = current2.lightoff, current.lighton = current2.lighton)
workingData = subset(old.lightoff, old.lightoff$time >=  - t_interval & old.lightoff$time < t_interval)
folder.old$off30data = normal.diy(workingData = workingData, baseline = baseline.old, 
                                  current.lightoff = old.lightoff, current.lighton = old.lighton)
# append the old dataset to two folders
folder1$off30dataAll = rbind(folder1$off30data, folder.old$off30data)
folder2$off30dataAll = rbind(folder2$off30data, folder.old$off30data)

rm(workingData)

# 3.3 Hotelling T-squared test
workingData = folder1$off30data
data.old = folder.old$off30data
OnOff = 'Light-Off'
# first pass folder
folder1$off30test = test.diy(workingData = folder1$off30data, data.old = folder.old$off30data, OnOff = 'Light-Off')
# second pass folder
folder2$off30test = test.diy(workingData = folder2$off30data, data.old = folder.old$off30data, OnOff = 'Light-Off')

# 3.4 Find consistent drugs

#### The First Pass Folder

# drugs significantly different from Q344X
folder1$off30sigQ344X = subset(folder1$off30test[,2], folder1$off30test[,2] < sig.lv1 )
# drugs significantly not different from RHO
folder1$off30sigRHO = subset(folder1$off30test[,4], folder1$off30test[,4] >= sig.lv2 )
# drugs significantly different from Q344X AND not different from RHO
folder1$off30sig = subset(folder1$off30test[,c(2,4)], folder1$off30test[,2] < sig.lv1 & folder1$off30test[,4] >= sig.lv2 )
# Normalized Mean Curves per drug per second
folder1$off30meanDrugTime = tapply(as.numeric(unlist(folder1$off30data$mean_int_normalized)), list(folder1$off30data$genotype, folder1$off30data$time), mean)    


#### The Second Pass Folder
# drugs significantly different from Q344X
folder2$off30sigQ344X = subset(folder2$off30test[,2], folder2$off30test[,2] < sig.lv1 )
# drugs significantly not different from RHO
folder2$off30sigRHO = subset(folder2$off30test[,4], folder2$off30test[,4] >= sig.lv2 )
# drugs significantly different from Q344X AND not different from RHO
folder2$off30sig = subset(folder2$off30test[,c(2,4)], folder2$off30test[,2] < sig.lv1 & folder2$off30test[,4] >= sig.lv2 )
# Normalized Mean Curves per drug per second
folder2$off30meanDrugTime = tapply(as.numeric(unlist(folder2$off30data$mean_int_normalized)), list(folder2$off30data$genotype, folder2$off30data$time), mean)    


#### Intersection of Two Folders
# drugs in common
drugs.common = intersect(folder1$off30meanDrugTime %>% rownames, folder2$off30meanDrugTime %>% rownames)
# differences between two folders per common drug
drugs.diff = folder1$off30meanDrugTime[drugs.common,] - folder2$off30meanDrugTime[drugs.common,]
# Euclidean distances between two folders per common drug
drugs.dist.2reps = sort(apply(drugs.diff^2, 1, sum) %>% sqrt, decreasing = F)

data1 = folder1$off30data
data2 = folder2$off30data
drugs.common = drugs.common
test.consistency.result = test.consistency(data1 = folder1$off30data, data2 = folder2$off30data, drugs.common = drugs.common)

drugs.consistent = list()
drugs.consistent$off30sigQ344X = intersect(names(folder1$off30sigQ344X), names(folder2$off30sigQ344X))
drugs.consistent$off30sigRHO = intersect(names(folder1$off30sigRHO), names(folder2$off30sigRHO))
drugs.consistent$off30sig = intersect(rownames(folder1$off30sig), rownames(folder2$off30sig))
drugs.consistent$off30ttest = sort(rownames(test.consistency.result)[test.consistency.result$fdr >= 0.05], decreasing = F)

write.csv(folder1$off30test, file = "Drugs_hotelling_rep1.csv")
write.csv(folder2$off30test, file = "Drugs_hotelling_rep2.csv")
write.csv(test.consistency.result, file = "Drugs_hdtest.csv")
