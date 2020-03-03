library(dplyr)
library(readr)
rm(list = ls())
setwd("~/zebrafish/analysis1218")



func1 <- function(file_tmp) {
  read_delim(file_tmp, "\t", escape_double = FALSE, trim_ws = TRUE)
}

# Chokh 
raw1 <- sapply("~/zebrafish/Chokh/20171024 Green Plate Results.XLS", func1, simplify = FALSE) %>%  bind_rows(.id = "id")
raw1$genotype %>% unique 
raw1.new = raw1[1:(nrow(raw1)/2)*2, c('location', 'genotype', 'start', 'end', 
                                      'inadist', 'smldist', 'lardist')]
raw1.new$mean = apply(raw1.new[,c('inadist', 'smldist', 'lardist')], 1, sum)
raw1.new = raw1.new[, -4]
raw1.new$start = floor(raw1.new$start)

# 30 secs before and after lighton / lightoff
lighton.shot = 1800
lightoff.shot = 5400
t_interval = 60

Chokh_Car_lighton = subset(raw1.new, raw1.new$start >= lighton.shot - t_interval & raw1.new$start < lighton.shot + t_interval)
Chokh_Car_lightoff = subset(raw1.new, raw1.new$start >= lightoff.shot - t_interval & raw1.new$start < lightoff.shot + t_interval)

# add time variable to light on/off dataset
Chokh_Car_lighton$time = rep(rep(c(-t_interval:(t_interval-1)), each = 96), nrow(Chokh_Car_lighton)/96/2/t_interval)
Chokh_Car_lightoff$time = rep(rep(c(-t_interval:(t_interval-1)), each = 96), nrow(Chokh_Car_lightoff)/96/2/t_interval)

# add batch id to light on/off dataset
Chokh_Car_lighton$rep = rep(1:(nrow(Chokh_Car_lighton)/96/2/t_interval), each = 96*2*t_interval)
Chokh_Car_lightoff$rep = rep(1:(nrow(Chokh_Car_lightoff)/96/2/t_interval), each = 96*2*t_interval)

# delete location variable in light on/off dataset, 9 variables in total
Chokh_Car_lighton = Chokh_Car_lighton[,-1]
Chokh_Car_lightoff = Chokh_Car_lightoff[,-1]

# add light intensity to light on/off dataset
light_intensity <- read.csv("~/zebrafish/light_intensity1.csv",header = T)

Chokh_Car_lighton$lightintensity = rep(light_intensity$Intensity, nrow(Chokh_Car_lighton)/96)
Chokh_Car_lightoff$lightintensity = rep(light_intensity$Intensity, nrow(Chokh_Car_lightoff)/96)

# null light intensity of data points in dark periods
Chokh_Car_lighton$lightintensity[Chokh_Car_lighton$start < lighton.shot] = 0
Chokh_Car_lightoff$lightintensity[Chokh_Car_lightoff$start >= lightoff.shot] = 0

# create baselines of 4 variables
light1on30 = subset(Chokh_Car_lighton, Chokh_Car_lighton$start >= lighton.shot - t_interval & Chokh_Car_lighton$start < lighton.shot )
# baseline, a data frame that stores the baseline for each genotype (drug)
baseline_car <- as.data.frame(cbind(
  with(light1on30, tapply(inadist, genotype, mean)),
  with(light1on30, tapply(smldist, genotype, mean)),
  with(light1on30, tapply(lardist, genotype, mean)),
  with(light1on30, tapply(mean, genotype, mean))))
colnames(baseline_car) = c('inadist', 'smldist', 'lardist', 'mean')

save(Chokh_Car_lightoff, baseline_car,Chokh_Car_lighton,file="chokh_car_60.Rdata")
