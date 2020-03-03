# WT+DMSO
load("C:/Users/mengr/Dropbox/VMR_Project/data0806/DMSO60.RData")
source("C:/Users/mengr/Dropbox/VMR_Project/analysis0903/function.R")

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

Rho_DMSO <- folder1$off30data[folder1$off30data$genotype == "Rho", ]
q344x_DMSO <- folder1$off30data[folder1$off30data$genotype == "Q344X", ]

# SEM Rho
meanDrugTime = tapply(as.numeric(unlist(Rho_DMSO$mean_int_normalized)), 
                      list(Rho_DMSO$genotype, Rho_DMSO$time), mean)
SE = tapply(as.numeric(unlist(Rho_DMSO$mean_int_normalized)), 
            list(Rho_DMSO$genotype, Rho_DMSO$time), sd)

SEM = SE/sqrt(48*9)
result.Rho_DMSO <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result.Rho_DMSO) <- c("mean.rho", "SE.rho", "SEM.rho")
# SEM  q344X
meanDrugTime = tapply(as.numeric(unlist(q344x_DMSO$mean_int_normalized)), 
                      list(q344x_DMSO$genotype, q344x_DMSO$time), mean)
SE = tapply(as.numeric(unlist(q344x_DMSO$mean_int_normalized)), 
            list(q344x_DMSO$genotype, q344x_DMSO$time), sd)

SEM = SE/sqrt(48*9)
result.q344x_DMSO <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result.q344x_DMSO) <- c("mean.q344x", "SE.q344x", "SEM.q344x")
DMSO = rbind(result.Rho_DMSO, result.q344x_DMSO)
write.csv(DMSO, file = "C:/Users/mengr/Dropbox/VMR_Project/analysis0903/DMSO_-30_to_59_Seconds.csv")

##################################################################################################
# 1 second analysis
folder1 = list()
folder2 = list()
folder.old = list()
compare_fold = list()
offset = 0.13
t_interval = 1
range_time = c(-t_interval:(t_interval-1))

workingData = subset(current1.lightoff, current1.lightoff$time >=  - t_interval & current1.lightoff$time < t_interval)
folder1$off30data = normal.diy(workingData = workingData, baseline = baseline1, 
                               current.lightoff = current1.lightoff, current.lighton = current1.lighton)

result = data.frame()
for (i in 1:9){
  rep <- folder1$off30data[folder1$off30data$rep == i,]
  rho = rep[rep$genotype == "Rho" & rep$time == 0,]
  q344 = rep[rep$genotype == "Q344X" & rep$time == 0,]
  ave_rho = mean(rho$mean_int_normalized)
  ave_q344 = mean(q344$mean_int_normalized)
  result[1,i] = ave_rho
  result[2,i] = ave_q344
  print(i)
}

colnames(result) = c("rep1","rep2","rep3","rep4","rep5","rep6","rep7","rep8","rep9")
rownames(result) = c("Rho", "Q344X")
write.csv(result,"C:/Users/mengr/Dropbox/VMR_Project/analysis0903/Time 0 for each of the 9 DMSO biological replicates.csv")


#################################################################################################
Rho_DMSO <- folder1$off30data[folder1$off30data$genotype == "Rho", ]
q344x_DMSO <- folder1$off30data[folder1$off30data$genotype == "Q344X", ]

### Rho
meanDrugTime = tapply(as.numeric(unlist(Rho_DMSO$mean_int_normalized)), 
                      list(Rho_DMSO$genotype, Rho_DMSO$time), mean)
SE = tapply(as.numeric(unlist(Rho_DMSO$mean_int_normalized)), 
            list(Rho_DMSO$genotype, Rho_DMSO$time), sd)

SEM = SE/sqrt(48*9)
result.Rho_DMSO <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result.Rho_DMSO) <- c("mean.rho", "SE.rho", "SEM.rho")

### Q344x
meanDrugTime = tapply(as.numeric(unlist(q344x_DMSO$mean_int_normalized)), 
                      list(q344x_DMSO$genotype, q344x_DMSO$time), mean)
SE = tapply(as.numeric(unlist(q344x_DMSO$mean_int_normalized)), 
            list(q344x_DMSO$genotype, q344x_DMSO$time), sd)

SEM = SE/sqrt(48*9)
result.q344x_DMSO <- as.data.frame(rbind(meanDrugTime,SE,SEM))
rownames(result.q344x_DMSO) <- c("mean.q344x", "SE.q344x", "SEM.q344x")
DMSO = rbind(result.Rho_DMSO, result.q344x_DMSO)

write.csv(DMSO, file = "C:/Users/mengr/Dropbox/VMR_Project/analysis0903/DMSO_1_Seconds.csv")
