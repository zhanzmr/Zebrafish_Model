load("C:/Users/mengr/Dropbox/VMR_Project/data0625/currentNew60.RData")
source("~/zebrafish/analysis1218/code/function_new.R")


folder_rep1 <- list()
folder_rep2 <- list()
#add offset
offset <- 0.14
t_interval <- 30
range_time <- c(-t_interval:(t_interval-1))
# range of y-axis
plot.range = c(-0.03, 0.43)

# rep1
workingData = subset(current1.lightoff, current1.lightoff$time >=  - t_interval & current1.lightoff$time < t_interval)
folder_rep1$off30data = normal.diy(workingData = workingData, baseline = baseline1, 
                                   current.lightoff = current1.lightoff, current.lighton = current1.lighton)
chosen1 = c('Drug_58', 'Drug_14') 
folder_rep1$off30fig = plot.diy(workingData = folder_rep1$off30data, plot.range = plot.range, chosen = chosen1, OnOff = 'Light-Off', rep = 21)

ggarrange(folder_rep1$off30fig$mean, folder_rep1$off30fig$mean_light_normalized,
          folder_rep1$off30fig$mean_baseline_normalized, folder_rep1$off30fig$mean_int_normalized, 
          ncol = 2, nrow = 2, labels = c("a)", "b)","c)","d)"))

# rep2
workingData = subset(current2.lightoff, current2.lightoff$time >=  - t_interval & current2.lightoff$time < t_interval)
folder_rep2$off30data = normal.diy(workingData = workingData, baseline = baseline2, 
                                   current.lightoff = current2.lightoff, current.lighton = current2.lighton)

chosen1 = c('Drug_58', 'Drug_14') 
folder_rep2$off30fig = plot.diy(workingData = folder_rep2$off30data, plot.range = plot.range, chosen = chosen1, OnOff = 'Light-Off', rep = 19)

ggarrange(folder_rep2$off30fig$mean, folder_rep2$off30fig$mean_light_normalized, folder_rep2$off30fig$mean_batch_normalized,
          folder_rep2$off30fig$mean_baseline_normalized, folder_rep2$off30fig$mean_int_normalized, 
          ncol = 2, nrow = 3, labels = c("a)", "b)","c)","d)","e)"))

folder_rep1$off30fig$mean_int_normalized
folder_rep2$off30fig$mean_int_normalized
ggarrange(folder_rep1$off30fig$mean_int_normalized, folder_rep2$off30fig$mean_int_normalized, 
          ncol = 2, nrow = 3, labels = c("a)", "b)","c)","d)","e)"))

## Hoteling t-test
workingData = folder_rep1$off30data[folder_rep1$off30data$genotype == "Drug_58",]
data.old = folder_rep1$off30data[folder_rep1$off30data$genotype == "Drug_14",]
result = test.diy2(workingData = workingData, data.old = data.old, OnOff = 'Light-Off')
result1 = as.data.frame(t(result))

