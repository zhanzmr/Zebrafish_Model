library(ggplot2)
library(reshape2)
library(ggpubr)
library(Hotelling)
library(knitr)
library(HDtest)
library(kableExtra)

rm(list = ls())
#setwd("~/zebrafish/analysis1218")
load("~/zebrafish/analysis1218/chokh60.Rdata")
load("~/zebrafish/analysis1218/chokh_car_60.Rdata")
source("~/zebrafish/analysis1218/code/function_new.R")

folder_chokh <- list()
folder_car <- list()
#add offset
offset <- 0.06
t_interval <- 30
range_time <- c(-t_interval:(t_interval-1))
# range of y-axis
plot.range = c(-0.03, 0.43)
# significance levels 
# for Q344X
#sig.lv1 = 0.05

# chokh

Chokh_lightoff$genotype = gsub("Cho","Chokh",Chokh_lightoff$genotype)
workingData = subset(Chokh_lightoff, Chokh_lightoff$time >=  - t_interval & Chokh_lightoff$time < t_interval)
folder_chokh$off30data = normal.diy(workingData = workingData, baseline = baseline1, 
                                   current.lightoff = Chokh_lightoff, current.lighton = Chokh_lighton)
# chosen1 = c('Chokh','1B6') 
# folder_chokh$off30fig = plot.diy(workingData = folder_chokh$off30data, plot.range = plot.range, chosen = chosen1, OnOff = 'Light-Off', rep = 1)
# 
# ggarrange(folder_chokh$off30fig$mean, folder_chokh$off30fig$mean_light_normalized,
#           folder_chokh$off30fig$mean_baseline_normalized, folder_chokh$off30fig$mean_int_normalized, 
#           ncol = 2, nrow = 2, labels = c("a)", "b)","c)","d)"))

# chokh_car
workingData = subset(Chokh_Car_lightoff, Chokh_Car_lightoff$time >=  - t_interval & Chokh_Car_lightoff$time < t_interval)
folder_car$off30data = normal.diy(workingData = workingData, baseline = baseline_car, 
                                    current.lightoff = Chokh_Car_lightoff, current.lighton = Chokh_Car_lighton)
# chosen1 = c('Chokh_Car','1D3','1D4') 
# folder_car$off30fig = plot.diy(workingData = folder_car$off30data, plot.range = plot.range, chosen = chosen1, OnOff = 'Light-Off', rep = 1)
# 
# ggarrange(folder_car$off30fig$mean, folder_car$off30fig$mean_light_normalized,
#           folder_car$off30fig$mean_baseline_normalized, folder_car$off30fig$mean_int_normalized, 
#           ncol = 2, nrow = 2, labels = c("a)", "b)","c)","d)"))

######## Hotelling test Chokh with Chokh_car

# Compare Chokh and Chokh_car
workingData = folder_chokh$off30data[folder_chokh$off30data$genotype == "Chokh",]
data.old = folder_car$off30data[folder_car$off30data$genotype == "Chokh_Car",]
result = test.diy2(workingData = workingData, data.old = data.old, OnOff = 'Light-Off')
result1 = as.data.frame(t(result))

# output excel
write.csv(result1, file = "result_Chokh_and_Chokh_car.csv")


###################################### average value from -30 to 59 seconds

# folder 1, 2 and old, including normalized data
folder_chokh <- list()
folder_car <- list()
#add offset
offset = 0.06
t_interval = 30
range_time = c(-t_interval:(t_interval-1))
plot.range = c(-0.03, 0.43)

Chokh_lightoff$genotype = gsub("Cho","Chokh",Chokh_lightoff$genotype)
workingData = subset(Chokh_lightoff, Chokh_lightoff$time >=  - t_interval & Chokh_lightoff$time < 2*t_interval)
folder_chokh$data = normal.diy(workingData = workingData, baseline = baseline1, 
                                   current.lightoff = Chokh_lightoff, current.lighton = Chokh_lighton)

workingData = subset(Chokh_Car_lightoff, Chokh_Car_lightoff$time >=  - t_interval & Chokh_Car_lightoff$time < 2*t_interval)
folder_car$data = normal.diy(workingData = workingData, baseline = baseline_car, 
                                    current.lightoff = Chokh_Car_lightoff, current.lighton = Chokh_Car_lighton)

chokh <- folder_chokh$data[folder_chokh$data$genotype == "Chokh", ]
chokh_car <- folder_car$data[folder_car$data$genotype == "Chokh_Car", ]

#data = folder_chokh$data
# plot chokh chokh_car
chokh_chokhcar = list()
workingData = rbind(chokh, chokh_car)

chosen1 = c('Chokh','Chokh_Car') 
chokh_chokhcar$fig = plot.diy(workingData = workingData, plot.range = plot.range, chosen = chosen1, OnOff = 'Light-Off', rep = 1)

ggarrange(chokh_chokhcar$fig$mean, chokh_chokhcar$fig$mean_light_normalized,
          chokh_chokhcar$fig$mean_baseline_normalized, chokh_chokhcar$fig$mean_int_normalized, 
          ncol = 2, nrow = 2, labels = c("a)", "b)","c)","d)"))

chokh_chokhcar$fig$mean_int_normalized

### chokh
meanDrugTime = tapply(as.numeric(unlist(chokh$mean_int_normalized)), 
                      list(chokh$genotype, chokh$time), mean)
SE = tapply(as.numeric(unlist(chokh$mean_int_normalized)), 
            list(chokh$genotype, chokh$time), sd)

SEM = SE/sqrt(48*max(Chokh_lightoff$rep))
result_chokh <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result_chokh) <- c("mean", "SE", "SEM")

### chokh_car
meanDrugTime = tapply(as.numeric(unlist(chokh_car$mean_int_normalized)), 
                      list(chokh_car$genotype, chokh_car$time), mean)
SE = tapply(as.numeric(unlist(chokh_car$mean_int_normalized)), 
            list(chokh_car$genotype, chokh_car$time), sd)

SEM = SE/sqrt(48*max(Chokh_Car_lightoff$rep))
result_chokh_car <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result_chokh_car) <- c("mean", "SE", "SEM")

### combine
result = rbind(result_chokh, result_chokh_car)
write.csv(result, file = "-30_to_59_chokh.csv")
