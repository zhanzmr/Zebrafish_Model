rm(list = ls())
setwd("/cloud/project")

library(ggplot2)
library(reshape2)
library(ggpubr)
library(Hotelling)
library(knitr)
library(HDtest)
library(kableExtra)

#setwd("~/zebrafish/analysis1218")
load("/cloud/project/environment/old60.RData")
source("/cloud/project/code/function_new.R")

# average value from -30 to 59 seconds
# folder 1, 2 and old, including normalized data
folder <- list()

#add offset
offset = 0.13
t_interval = 30
range_time = c(-t_interval:(2*t_interval-1))
plot.range = c(-0.03, 0.43)

workingData = subset(current1.lightoff, current1.lightoff$time >=  - t_interval & current1.lightoff$time < 2*t_interval)
folder$data = normal.diy(workingData = workingData, baseline = baseline1, current.lightoff = current1.lightoff)

# plot
chosen1 = c('Q344X','Rho') 
plot = list()
plot$fig = plot.diy(workingData = folder$data, plot.range = plot.range, chosen = chosen1, OnOff = 'Light-Off', rep = 18)

ggarrange(plot$fig$mean, plot$fig$mean_light_normalized,
          plot$fig$mean_baseline_normalized, plot$fig$mean_int_normalized, 
          ncol = 2, nrow = 2, labels = c("a)", "b)","c)","d)"))

Rho <- folder$data[folder$data$genotype == "Rho", ]
Q344X <- folder$data[folder$data$genotype == "Q344X", ]


### Rho
meanDrugTime = tapply(as.numeric(unlist(Rho$mean_int_normalized)), 
                      list(Rho$genotype, Rho$time), mean)
SE = tapply(as.numeric(unlist(Rho$mean_int_normalized)), 
            list(Rho$genotype, Rho$time), sd)

SEM = SE/sqrt(48*max(Rho$rep))
result_Rho <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result_Rho) <- c("mean", "SE", "SEM")

write.csv(result_Rho, file = "result_Rho.csv")
### Q344X
meanDrugTime = tapply(as.numeric(unlist(Q344X$mean_int_normalized)), 
                      list(Q344X$genotype, Q344X$time), mean)
SE = tapply(as.numeric(unlist(Q344X$mean_int_normalized)), 
            list(Q344X$genotype, Q344X$time), sd)

SEM = SE/sqrt(48*max(Q344X$rep))
result_Q344X <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result_Q344X) <- c("mean", "SE", "SEM")
write.csv(result_Q344X, file = "result_Q344X.csv")

#### 1 second after light onset ##

# 1 second analysis
folder1 = list()
offset = 0.13
t_interval = 1
range_time = c(-t_interval:(t_interval-1))

workingData = subset(current1.lightoff, current1.lightoff$time >=  - t_interval & current1.lightoff$time < t_interval)
folder1$off30data = normal.diy(workingData = workingData, baseline = baseline1, current.lightoff = current1.lightoff)


plot.range = c(-0.03, 0.43)
chosen1 = c('Q344X','Rho') 
plot = list()
plot$fig = plot.diy(workingData = folder1$off30data, plot.range = plot.range, chosen = chosen1, OnOff = 'Light-Off', rep = 18)

ggarrange(plot$fig$mean, plot$fig$mean_light_normalized,
          plot$fig$mean_baseline_normalized, plot$fig$mean_int_normalized, 
          ncol = 2, nrow = 2, labels = c("a)", "b)","c)","d)"))



result = data.frame()
for (i in 1:18){
  rep <- folder1$off30data[folder1$off30data$rep == i,]
  rho = rep[rep$genotype == "Rho" & rep$time == 0,]
  q344 = rep[rep$genotype == "Q344X" & rep$time == 0,]
  ave_rho = mean(rho$mean_int_normalized)
  ave_q344 = mean(q344$mean_int_normalized)
  result[1,i] = ave_rho
  result[2,i] = ave_q344
  print(i)
}

colnames(result) = c("rep1","rep2","rep3","rep4","rep5","rep6","rep7","rep8","rep9",
                     "rep10","rep11","rep12","rep13","rep14","rep15","rep16","rep17","rep18")
rownames(result) = c("Rho", "Q344X")
write.csv(result,"Time_0_original_Rho_Q344X.csv")

















