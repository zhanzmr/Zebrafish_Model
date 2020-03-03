library(ggplot2)
library(reshape2)
library(ggpubr)
library(Hotelling)
library(knitr)
library(HDtest)
library(kableExtra)

rm(list = ls())
#setwd("~/zebrafish/analysis1218")
load("~/zebrafish/analysis1218/ADCY60.Rdata")
load("~/zebrafish/analysis1218/dmso_new.Rdata")
source("~/zebrafish/analysis1218/code/function_new.R")

folder_adcy <- list()
folder_dmso <- list()
#add offset
offset <- 0.06
t_interval <- 30
range_time <- c(-t_interval:(t_interval-1))
# range of y-axis
plot.range = c(-0.03, 0.43)
# significance levels 
# for Q344X
#sig.lv1 = 0.05

# ADCY
workingData = subset(current1.lightoff, current1.lightoff$time >=  - t_interval & current1.lightoff$time < t_interval)
folder_adcy$off30data = normal.diy(workingData = workingData, baseline = baseline1, 
                               current.lightoff = current1.lightoff, current.lighton = current1.lighton)
chosen1 = c('SQ1', 'SQ2') 
folder_adcy$off30fig = plot.diy(workingData = folder_adcy$off30data, plot.range = plot.range, chosen = chosen1, OnOff = 'Light-Off', rep = 1)

ggarrange(folder_adcy$off30fig$mean, folder_adcy$off30fig$mean_light_normalized,
          folder_adcy$off30fig$mean_baseline_normalized, folder_adcy$off30fig$mean_int_normalized, 
          ncol = 2, nrow = 2, labels = c("a)", "b)","c)","d)"))

# DMSO
offset <- 0.13
workingData = subset(lightoff_dmso, lightoff_dmso$time >=  - t_interval & lightoff_dmso$time < t_interval)
folder_dmso$off30data = normal.diy(workingData = workingData, baseline = baseline1, 
                                   current.lightoff = lightoff_dmso, current.lighton = lighton_dmso)

chosen1 = c('Rho', 'Q344X') 
folder_dmso$off30fig = plot.diy(workingData = folder_dmso$off30data, plot.range = plot.range, chosen = chosen1, OnOff = 'Light-Off', rep = 9)

ggarrange(folder_dmso$off30fig$mean, folder_dmso$off30fig$mean_light_normalized, folder_dmso$off30fig$mean_batch_normalized,
          folder_dmso$off30fig$mean_baseline_normalized, folder_dmso$off30fig$mean_int_normalized, 
          ncol = 2, nrow = 3, labels = c("a)", "b)","c)","d)","e)"))

######## Hotelling test with DMSO dataset

# Compare SQ1 and Q344XDMSO1
workingData = folder_adcy$off30data[folder_adcy$off30data$genotype == "SQ1",]
data.old = folder_adcy$off30data[folder_adcy$off30data$genotype == "Q344XDMSO1",]
result = test.diy2(workingData = workingData, data.old = data.old, OnOff = 'Light-Off')
result1 = as.data.frame(t(result))

# Compare SQ2 and Q344XDMSO2
workingData = folder_adcy$off30data[folder_adcy$off30data$genotype == "SQ2",]
data.old = folder_adcy$off30data[folder_adcy$off30data$genotype == "Q344XDMSO2",]
result = test.diy2(workingData = workingData, data.old = data.old, OnOff = 'Light-Off')
result2 = as.data.frame(t(result))

# Compare SQ1 and DMSO_Q344X
workingData = folder_adcy$off30data[folder_adcy$off30data$genotype == "SQ1",]
data.old = folder_dmso$off30data[folder_dmso$off30data$genotype == "Q344X",]
result = test.diy2(workingData = workingData, data.old = data.old, OnOff = 'Light-Off')
result3 = as.data.frame(t(result))

# Compare SQ2 and DMSO_Q344X
workingData = folder_adcy$off30data[folder_adcy$off30data$genotype == "SQ2",]
data.old = folder_dmso$off30data[folder_dmso$off30data$genotype == "Q344X",]
result = test.diy2(workingData = workingData, data.old = data.old, OnOff = 'Light-Off')
result4 = as.data.frame(t(result))

# output excel
hotel_result = rbind(result1, result2, result3, result4)
write.csv(hotel_result, file = "result_SQ_and_dmso.csv")

# average value from -30 to 59 seconds

# folder 1, 2 and old, including normalized data
folder_adcy <- list()
folder_dmso <- list()
#add offset
offset = 0.06
t_interval = 30
range_time = c(-t_interval:(t_interval-1))
plot.range = c(-0.03, 0.43)

workingData = subset(current1.lightoff, current1.lightoff$time >=  - t_interval & current1.lightoff$time < 2*t_interval)
folder_adcy$data = normal.diy(workingData = workingData, baseline = baseline1, 
                               current.lightoff = current1.lightoff, current.lighton = current1.lighton)

offset = 0.13
workingData = subset(lightoff_dmso, lightoff_dmso$time >=  - t_interval & lightoff_dmso$time < 2*t_interval)
folder_dmso$data = normal.diy(workingData = workingData, baseline = baseline1, 
                                   current.lightoff = lightoff_dmso, current.lighton = lighton_dmso)

SQ1 <- folder_adcy$data[folder_adcy$data$genotype == "SQ1", ]
SQ2 <- folder_adcy$data[folder_adcy$data$genotype == "SQ2", ]
Q344XDMSO1 <- folder_adcy$data[folder_adcy$data$genotype == "Q344XDMSO1", ]
Q344XDMSO2 <- folder_adcy$data[folder_adcy$data$genotype == "Q344XDMSO2", ]

Q344X_DMSO <- folder_dmso$data[folder_dmso$data$genotype == "Q344X", ]
Rho_DMSO <- folder_dmso$data[folder_dmso$data$genotype == "Rho", ]

# plot Q344X_SQ2, Q344X_DMSO # Q344X_DMSO[,-11]
SQ2_DMSO = list()
workingData = rbind(SQ1, Q344XDMSO1)

chosen1 = c('SQ1','Q344XDMSO1') 
SQ2_DMSO$fig = plot.diy(workingData = workingData, plot.range = plot.range, chosen = chosen1, OnOff = 'Light-Off', rep = 1)

ggarrange(SQ2_DMSO$fig$mean, SQ2_DMSO$fig$mean_light_normalized,
          SQ2_DMSO$fig$mean_baseline_normalized, SQ2_DMSO$fig$mean_int_normalized, 
          ncol = 2, nrow = 2, labels = c("a)", "b)","c)","d)"))

SQ2_DMSO$fig$mean_int_normalized

### SQ1
meanDrugTime = tapply(as.numeric(unlist(SQ1$mean_int_normalized)), 
                      list(SQ1$genotype, SQ1$time), mean)
SE = tapply(as.numeric(unlist(SQ1$mean_int_normalized)), 
            list(SQ1$genotype, SQ1$time), sd)

SEM = SE/sqrt(48*max(SQ1$rep))
result_SQ1 <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result_SQ1) <- c("mean", "SE", "SEM")

### SQ2
meanDrugTime = tapply(as.numeric(unlist(SQ2$mean_int_normalized)), 
                      list(SQ2$genotype, SQ2$time), mean)
SE = tapply(as.numeric(unlist(SQ2$mean_int_normalized)), 
            list(SQ2$genotype, SQ2$time), sd)

SEM = SE/sqrt(48*max(SQ2$rep))
result_SQ2 <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result_SQ2) <- c("mean", "SE", "SEM")

### Q344XDMSO1
meanDrugTime = tapply(as.numeric(unlist(Q344XDMSO1$mean_int_normalized)), 
                      list(Q344XDMSO1$genotype, Q344XDMSO1$time), mean)
SE = tapply(as.numeric(unlist(Q344XDMSO1$mean_int_normalized)), 
            list(Q344XDMSO1$genotype, Q344XDMSO1$time), sd)

SEM = SE/sqrt(48*max(Q344XDMSO1$rep))
result_Q344XDMSO1 <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result_Q344XDMSO1) <- c("mean", "SE", "SEM")

### Q344XDMSO2
meanDrugTime = tapply(as.numeric(unlist(Q344XDMSO2$mean_int_normalized)), 
                      list(Q344XDMSO2$genotype, Q344XDMSO2$time), mean)
SE = tapply(as.numeric(unlist(Q344XDMSO2$mean_int_normalized)), 
            list(Q344XDMSO2$genotype, Q344XDMSO2$time), sd)

SEM = SE/sqrt(48*max(Q344XDMSO2$rep))
result_Q344XDMSO2 <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result_Q344XDMSO2) <- c("mean", "SE", "SEM")

### combine
result = rbind(result_SQ1, result_SQ2, result_Q344XDMSO1, result_Q344XDMSO2)
write.csv(result, file = "-30_to_59_ADCY.csv")


### DMSO_rho
meanDrugTime = tapply(as.numeric(unlist(Rho_DMSO$mean_int_normalized)), 
                      list(Rho_DMSO$genotype, Rho_DMSO$time), mean)
SE = tapply(as.numeric(unlist(Rho_DMSO$mean_int_normalized)), 
            list(Rho_DMSO$genotype, Rho_DMSO$time), sd)

SEM = SE/sqrt(48*max(Rho_DMSO$rep))
result_Rho_DMSO <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result_Rho_DMSO) <- c("mean", "SE", "SEM")

### DMSO_Q344X
meanDrugTime = tapply(as.numeric(unlist(Q344X_DMSO$mean_int_normalized)), 
                      list(Q344X_DMSO$genotype, Q344X_DMSO$time), mean)
SE = tapply(as.numeric(unlist(Q344X_DMSO$mean_int_normalized)), 
            list(Q344X_DMSO$genotype, Q344X_DMSO$time), sd)

SEM = SE/sqrt(48*max(Q344X_DMSO$rep))
result_Q344X_DMSO <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result_Q344X_DMSO) <- c("mean", "SE", "SEM")

### combine
result = rbind(result_Rho_DMSO, result_Q344X_DMSO)
write.csv(result, file = "-30_to_59_DMSO.csv")

# # hotelling t-test
# lightoff <- folder_adcy$off30data
# SQ <- lightoff[lightoff$genotype == "SQ1" | lightoff$genotype == "SQ2", ]
# dmso <- lightoff[lightoff$genotype == "Q344XDMSO1" | lightoff$genotype == "Q344XDMSO2", ]
# 
# # Compare original Q344x and Q344x + NTR
# result = test.diy2(workingData = SQ, data.old = dmso, OnOff = 'Light-Off')
# result = as.data.frame(t(result))
# write.csv(result, file = "result_SQ_and_dmso.csv")
# 
# # Compare original Q344x and Q344x + NTR
# workingData = folder1$off30data[folder1$off30data$genotype == "Neg",]
# data.old = folder.old$off30data[folder.old$off30data$genotype == "Q344X",]
# result = test.diy2(workingData = workingData, data.old = data.old, OnOff = 'Light-Off')
# result = as.data.frame(t(result))

