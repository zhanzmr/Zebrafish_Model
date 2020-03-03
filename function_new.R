# normalization method
normal.diy <- function(workingData, baseline, current.lightoff){
  ## 1. light normalization ----
  lm.light <- lm(mean ~ lightintensity, data = workingData)
  workingData$mean_light_normalized = lm.light$residuals + offset
  rm(lm.light)
  if (max(workingData$rep) > 1) {
      ## 2. Batch normalization ----
      lm.batch <- lm(mean ~ factor(rep), data = workingData)
      workingData$mean_batch_normalized = lm.batch$residuals + offset
      rm(lm.batch)
  }
  ## 3. Baseline normalization ----
  grandmean <- apply(baseline, 2, mean)
  lighton30before <- subset(current.lightoff, current.lightoff$time < 0)
  baseline30beforeOff = as.data.frame(with(lighton30before, tapply(mean, genotype, mean)))
  colnames(baseline30beforeOff) <- 'mean'
  # create mean normalized variable
  workingData$mean_baseline_normalized = workingData$mean
  # normalize Mean variable according to drug types
  for (ij in rownames(baseline)) {
    workingData$mean_baseline_normalized[workingData$genotype == ij] = workingData$mean[workingData$genotype == ij] - baseline30beforeOff[ij, 'mean'] + grandmean['mean'] + offset
  }
  rm(lighton30before)
  ## 4. Integrated normalization ----
  # Integrated normlalization for 4 variables
  workingData$mean_int_normalized <- workingData$mean_light_normalized
  if (max(workingData$rep) > 1) {
    lm.integrated = lm(mean_light_normalized ~ factor(rep), data = workingData)
    workingData$mean_int_normalized = lm.integrated$residuals
    rm(lm.integrated)
  }
  # normalize mean variable according to drug types
  for (ij in rownames(baseline30beforeOff)) {
    workingData$mean_int_normalized[workingData$genotype == ij] = workingData$mean_int_normalized[workingData$genotype == ij] - baseline30beforeOff[ij, 'mean'] + grandmean['mean'] + offset
  }
  
  rm(baseline30beforeOff)
  return(workingData)
}

# 1.2.1 function to draw a plot with ribbons
drawPlot = function(meanDF, sdDF, title, limits, rep) {
  melted.df = melt(meanDF)
  names(melted.df)[1:3] = c('genotype', 'time', 'average')
  melted.df$sd = as.vector( sdDF /sqrt(48*rep))
  
  
  # alternative denominator
  # sqrt( nrow(current.lighton)/t_interval/2 )
  
  # draw a plot
  gra_light = ggplot(melted.df, aes(x=time, y=average, color=genotype)) + 
    geom_line() + geom_ribbon(aes(x=time, ymin=average-sd, ymax=average+sd), 
                              alpha = 0.1, linetype=0) +
    scale_color_discrete(name="Genotype") +
    scale_x_continuous(name="Time(s)") +
    scale_y_continuous(name="Activity(cm)", limits=limits) +
    theme_set(theme_grey(base_size = 15)) +
    ## change input here  ==========------------------------------>
    labs(title = title)+
    theme(legend.position=c(.8,0.7)) + theme_bw()
  return(gra_light)
}

# 1.2.2 function to draw a plot with ribbons only with selected drugs
drawVar = function(dataset, variable, chosen, title, limits, rep) {
  ## Calculate mean and standard deviation for mean variable 
  meanDrugTime = tapply(as.numeric(unlist(dataset[,variable])), list(dataset$genotype, 
                                                                     dataset$time), mean)
  sdDrugTime = tapply(as.numeric(unlist(dataset[,variable])), list(dataset$genotype, 
                                                                   dataset$time), sd)
  # indices of selected drugs
  chosen.idx = sort(match(chosen, row.names(meanDrugTime)))
  
  # selected variable
  drawPlot(meanDrugTime[chosen.idx,], sdDrugTime[chosen.idx,], 
           title = title, limits=limits, rep)
}

# 1.2.3 function plot.diy, outputting figures
plot.diy = function(workingData, plot.range, chosen, OnOff = 'Light-Off', rep){
  fig = list() # a list storing figures
  # raw Mean variable
  fig$mean = drawVar(workingData, 'mean', limits=plot.range, 
                     chosen = chosen, title = paste(OnOff, '(Mean Original)'), rep)
  fig$mean_light_normalized = drawVar(workingData, 'mean_light_normalized', 
                                      limits=plot.range, chosen = chosen, 
                                      title = paste(OnOff, '(Mean Light Normalization)'), rep)
  if (rep > 1) {
    fig$mean_batch_normalized = drawVar(workingData, 'mean_batch_normalized', 
                                      limits=plot.range, chosen = chosen, 
                                      title = paste(OnOff, '(Mean Batch Normalization)'), rep)
  }
  fig$mean_baseline_normalized = drawVar(workingData, 'mean_baseline_normalized', 
                                         limits=plot.range, chosen = chosen, 
                                         title = paste(OnOff, '(Mean Baseline Normalization)'), rep)
  fig$mean_int_normalized = drawVar(workingData, 'mean_int_normalized', 
                                    limits=plot.range, chosen = chosen, 
                                    title = paste(OnOff, '(Mean Integrated Normalization)'), rep)
  return(fig)
}

### Hotelling t-test

func_wide = function(data){
  matrix(data[,2], ncol = t_interval*2)
}

test.diy2 = function(workingData, data.old, OnOff = 'Light-Off'){
  
  before.time = 1:t_interval
  after.time = 1:t_interval + t_interval
  
  # list of drugs for each folder, Q344X and Rho included at the end
  
  # p-value results for Mean variable
  pvalue.mat = rep(0,8)
  names(pvalue.mat) = c(paste('Before', OnOff),paste('BeforeT', OnOff), paste('BeforeDF', OnOff), paste('BeforeDF2', OnOff), paste('After', OnOff),paste('AfterT', OnOff), paste('AfterDF', OnOff), paste('AfterDF2', OnOff))
  # pvalue.mat = array(0, dim = c(1, 2), dimnames = list(c(paste('Before', OnOff), paste('After', OnOff))))
  
  ji = 'mean'
  jii = paste0(ji, '_int_normalized') # normalized
  
  # Normalized
  control1 = data.old[, c('time', jii)]
  # control1 is first Q344X and then Rho in the ijk loop
  dt.n1 = func_wide(as.data.frame(control1))
  drug.tmp = workingData[, c('time', jii)]
  dt.n2 = func_wide(as.data.frame(drug.tmp))
  
  # before lightoff
  mytest = hotelling.test(as.matrix(dt.n1[,before.time]), 
                          as.matrix(dt.n2[,before.time]), shrinkage = T)
  pvalue.mat[paste0('Before ', OnOff)] = mytest$pval
  pvalue.mat[paste('BeforeDF', OnOff)] = mytest$stats$df[1]
  pvalue.mat[paste('BeforeDF2', OnOff)] = mytest$stats$df[2]
  pvalue.mat[paste('BeforeT', OnOff)] = mytest$stats$statistic
  # after lightoff
  mytest = hotelling.test(as.matrix(dt.n1[,after.time]), 
                          as.matrix(dt.n2[,after.time]), shrinkage = T)
  pvalue.mat[paste0('After ', OnOff)] = mytest$pval
  pvalue.mat[paste('AfterDF', OnOff)] = mytest$stats$df[1]
  pvalue.mat[paste('AfterDF2', OnOff)] = mytest$stats$df[2]
  pvalue.mat[paste('AfterT', OnOff)] = mytest$stats$statistic
  
  # colnames(pvalue.mat) = paste(c('Before', 'After'), OnOff)
  # print out
  rm(mytest, drug.tmp, dt.n1, dt.n2)
  return(pvalue.mat)
}


plot_int = function(workingData, plot.range, chosen, OnOff = 'Light-Off', rep){
  fig$mean_int_normalized = drawVar(workingData, 'mean_int_normalized', 
                                    limits=plot.range, chosen = chosen, 
                                    title = paste(OnOff, '(Mean Integrated Normalization)'), rep)
}
  