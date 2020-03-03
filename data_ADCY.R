library(dplyr)
library(readr)
rm(list = ls())
setwd("/mnt/c/Users/mengr/Documents/zebrafish/analysis1218")
# 21 files in the 1st pass folder
files1st <- list.files(path = "/mnt/c/Users/mengr/Documents/zebrafish/ADCY/",
                      pattern = "*.XLS", full.names = T)

func1 <- function(file_tmp) {
  read_delim(file_tmp, "\t", escape_double = FALSE, trim_ws = TRUE)
}

# load data
raw1 <- sapply(files1st, func1, simplify = FALSE) %>%  bind_rows(.id = "id")
raw1$genotype %>% unique # only Q344x and Rho

raw1.new = raw1[1:(nrow(raw1)/2)*2, c('location', 'genotype', 'start', 'end', 
                                      'inadist', 'smldist', 'lardist')]
raw1.new$mean = apply(raw1.new[,c('inadist', 'smldist', 'lardist')], 1, sum)
raw1.new = raw1.new[, -4]

raw1.new$start = floor(raw1.new$start)

# 30 secs before and after lighton / lightoff
lighton.shot = 1800
lightoff.shot = 5400
t_interval = 60

current1.lighton = subset(raw1.new, raw1.new$start >= lighton.shot - t_interval & raw1.new$start < lighton.shot + t_interval)
current1.lightoff = subset(raw1.new, raw1.new$start >= lightoff.shot - t_interval & raw1.new$start < lightoff.shot + t_interval)
# add time variable to light on/off dataset
current1.lighton$time = rep(rep(c(-t_interval:(t_interval-1)), each = 96), nrow(current1.lighton)/96/2/t_interval)
current1.lightoff$time = rep(rep(c(-t_interval:(t_interval-1)), each = 96), nrow(current1.lightoff)/96/2/t_interval)

# add batch id to light on/off dataset
current1.lighton$rep = rep(1:(nrow(current1.lighton)/96/2/t_interval), each = 96*2*t_interval)
current1.lightoff$rep = rep(1:(nrow(current1.lightoff)/96/2/t_interval), each = 96*2*t_interval)

# delete location variable in light on/off dataset, 9 variables in total
current1.lighton = current1.lighton[,-1]
current1.lightoff = current1.lightoff[,-1]

# add light intensity to light on/off dataset
light_intensity <- read.csv("/mnt/c/Users/mengr/Documents/zebrafish/light_intensity1.csv",header = T)

current1.lighton$lightintensity = rep(light_intensity$Intensity, nrow(current1.lighton)/96)
current1.lightoff$lightintensity = rep(light_intensity$Intensity, nrow(current1.lightoff)/96)

# null light intensity of data points in dark periods
current1.lighton$lightintensity[current1.lighton$start < lighton.shot] = 0
current1.lightoff$lightintensity[current1.lightoff$start >= lightoff.shot] = 0

levels(factor(current1.lightoff$genotype))
# create baselines of 4 variables
light1on30 = subset(current1.lighton, current1.lighton$start >= lighton.shot - t_interval & current1.lighton$start < lighton.shot )
# baseline, a data frame that stores the baseline for each genotype (drug)
baseline1 <- as.data.frame(cbind(
  with(light1on30, tapply(inadist, genotype, mean)),
  with(light1on30, tapply(smldist, genotype, mean)),
  with(light1on30, tapply(lardist, genotype, mean)),
  with(light1on30, tapply(mean, genotype, mean))))
colnames(baseline1) = c('inadist', 'smldist', 'lardist', 'mean')

save(current1.lightoff, baseline1,current1.lighton,file="ADCY60.Rdata")