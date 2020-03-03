load("C:/Users/mengr/Dropbox/VMR_Project/data0806/NTR_METZ3_60.RData")
source("C:/Users/mengr/Dropbox/VMR_Project/data0806/function.R")

# folder 1, 2 and old, including normalized data
folder1 = list()
folder2 = list()
folder.old = list()
compare_fold = list()
#add offset
offset = 0.13
t_interval = 30
range_time = c(-t_interval:(t_interval-1))

workingData = subset(current1.lightoff, current1.lightoff$time >=  - t_interval & current1.lightoff$time < 2*t_interval)
folder1$off30data = normal.diy(workingData = workingData, baseline = baseline1, 
                               current.lightoff = current1.lightoff, current.lighton = current1.lighton)

Pos_NTR <- folder1$off30data[folder1$off30data$genotype == "Pos", ]
Neg_NTR <- folder1$off30data[folder1$off30data$genotype == "Neg", ]


meanDrugTime = tapply(as.numeric(unlist(Pos_NTR$mean_int_normalized)), 
                      list(Pos_NTR$genotype, Pos_NTR$time), mean)
# SEM
SE = tapply(as.numeric(unlist(Pos_NTR$mean_int_normalized)), 
            list(Pos_NTR$genotype, Pos_NTR$time), sd)

SEM = SE/sqrt(48*max(current1.lightoff$rep))
result.Pos_NTR <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result.Pos_NTR) <- c("mean.pos", "SE.pos", "SEM.pos")
######################
meanDrugTime = tapply(as.numeric(unlist(Neg_NTR$mean_int_normalized)), 
                      list(Neg_NTR$genotype, Neg_NTR$time), mean)
# SEM
SE = tapply(as.numeric(unlist(Neg_NTR$mean_int_normalized)), 
            list(Neg_NTR$genotype, Neg_NTR$time), sd)

SEM = SE/sqrt(48*max(current1.lightoff$rep))
result.Neg_NTR <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result.Neg_NTR) <- c("mean.neg", "SE.neg", "SEM.neg")
NTR = rbind(result.Pos_NTR, result.Neg_NTR)
write.csv(NTR, file = "C:/Users/mengr/Dropbox/VMR_Project/analysis0903/NTR_-30_to_59_Seconds.csv")


