# 1.1 function normal.diy, performing normalization
normal.diy = function(workingData, baseline, current.lightoff, current.lighton){
  
  ## 1. light normalization ----
  lm.light = lm(mean ~ lightintensity, data = workingData)
  workingData$mean_light_normalized = lm.light$residuals
  rm(lm.light)
  
  ## 2. Batch normalization ----
  lm.batch = lm(mean ~ factor(rep), data = workingData)
  workingData$mean_batch_normalized = lm.batch$residuals
  rm(lm.batch)
  
  ## 3. Baseline normalization ----
  grandmean = apply(baseline, 2, mean)
  lighton30before = subset(current.lightoff, current.lighton$time < 0)
  # dataframe, baseline for per drug/displacement 75*1
  baseline30beforeOff = as.data.frame(with(lighton30before, tapply(mean, genotype, mean)))
  colnames(baseline30beforeOff) = 'mean'
  # create mean normalized variable
  workingData$mean_baseline_normalized = workingData$mean
  
  # normalize Mean variable according to drug types
  for (ij in rownames(baseline)) {
    workingData$mean_baseline_normalized[workingData$genotype == ij] = workingData$mean[workingData$genotype == ij] - baseline30beforeOff[ij, 'mean'] + grandmean['mean'] 
  }
  rm(lighton30before)
  
  ## 4. Integrated normalization ----
  # Integrated normlalization for 4 variables
  lm.integrated = lm(mean_light_normalized ~ factor(rep), data = workingData)
  workingData$mean_int_normalized = lm.integrated$residuals
  
  # normalize mean variable according to drug types
  for (ij in rownames(baseline30beforeOff)) {
    workingData$mean_int_normalized[workingData$genotype == ij] = workingData$mean_int_normalized[workingData$genotype == ij] - baseline30beforeOff[ij, 'mean'] + grandmean['mean']
  }
  
  rm(baseline30beforeOff, lm.integrated)
  return(workingData)
}

# 1.2 functions outputting figures

# 1.2.1 function to draw a plot with ribbons
drawPlot = function(meanDF, sdDF, title, limits) {
  melted.df = melt(meanDF)
  names(melted.df)[1:3] = c('genotype', 'time', 'average')
  melted.df$sd = as.vector( sdDF / 10)
  
  
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
drawVar = function(dataset, variable, chosen, title, limits) {
  ## Calculate mean and standard deviation for mean variable 
  meanDrugTime = tapply(as.numeric(unlist(dataset[,variable])), list(dataset$genotype, 
                                                                     dataset$time), mean)
  sdDrugTime = tapply(as.numeric(unlist(dataset[,variable])), list(dataset$genotype, 
                                                                   dataset$time), sd)
  # indices of selected drugs
  chosen.idx = sort(match(chosen, row.names(meanDrugTime)))
  
  # selected variable
  drawPlot(meanDrugTime[chosen.idx,], sdDrugTime[chosen.idx,], 
           title = title, limits=limits)
}

# 1.2.3 function plot.diy, outputting figures
plot.diy = function(workingData, plot.range, chosen, OnOff = 'Light-Off'){
  fig = list() # a list storing figures
  # raw Mean variable
  fig$mean = drawVar(workingData, 'mean', limits=plot.range, 
                     chosen = chosen, title = paste(OnOff, '(Mean Original)'))
  fig$mean_light_normalized = drawVar(workingData, 'mean_light_normalized', 
                                      limits=plot.range, chosen = chosen, 
                                      title = paste(OnOff, '(Mean Light Normalization)'))
  fig$mean_batch_normalized = drawVar(workingData, 'mean_batch_normalized', 
                                      limits=plot.range, chosen = chosen, 
                                      title = paste(OnOff, '(Mean Batch Normalization)'))
  fig$mean_baseline_normalized = drawVar(workingData, 'mean_baseline_normalized', 
                                         limits=plot.range, chosen = chosen, 
                                         title = paste(OnOff, '(Mean Baseline Normalization)'))
  fig$mean_int_normalized = drawVar(workingData, 'mean_int_normalized', 
                                    limits=plot.range, chosen = chosen, 
                                    title = paste(OnOff, '(Mean Integrated Normalization)'))
  return(fig)
}

# 1.3 function test.diy, conducting the Hotelling test
# func_wide = function(data){
#     dt.wide = matrix(0, nrow= min(summary(factor(data[,1]))), ncol = (2*t_interval) )
#     for(i in range_time){
#         temp = data[ data[,1] == i, ]
#         dt.wide[,i+t_interval+1] = temp[1:min(summary(factor(data[,1]))),2]
#     }
#     return(dt.wide)
# }
func_wide = function(data){
  matrix(data[,2], ncol = t_interval*2)
}
test.diy = function(workingData, data.old, OnOff = 'Light-Off'){
  
  before.time = 1:t_interval
  after.time = 1:t_interval + t_interval
  
  # list of drugs for each folder, Q344X and Rho included at the end
  DrugList = sort(unique(workingData$genotype), decreasing = F)
  
  # p-value results for Mean variable
  result = array(0, dim = c(length(DrugList), 8, 2), dimnames = list(DrugList, c(paste('Before', OnOff), paste('BeforeT', OnOff), paste('BeforeDF1', OnOff),paste('BeforeDF2', OnOff),paste('After', OnOff),paste('AfterT', OnOff),paste('AfterDF1', OnOff),paste('AfterDF2', OnOff)), c('Q344X', 'Rho')))
  pvalue.mat = array(0, dim = c(length(DrugList), 2, 2), dimnames = list(DrugList, c(paste('Before', OnOff),paste('After', OnOff)), c('Q344X', 'Rho')))
  # pvalue.mat = array(0, dim = c(length(DrugList), 4, 2), dimnames = list(DrugList, c('Original Before Lightoff', 'Normalized Before Lightoff', 'Original After Lightoff', 'Normalized After Lightoff'), c('Q344X', 'Rho')))
  
  for(ijk in c('Q344X', 'Rho')){
    ji = 'mean'
    jii = paste0(ji, '_int_normalized') # normalized
    
    # Normalized
    control1 = data.old[data.old$genotype == ijk, c('time', jii)]
    # control1 is first Q344X and then Rho in the ijk loop
    dt.n1 = func_wide(as.data.frame(control1))
    for (ij in DrugList) {
      drug.tmp = workingData[workingData$genotype == ij, c('time', jii)]
      dt.n2 = func_wide(as.data.frame(drug.tmp))
      
      # before lightoff
      mytest = hotelling.test(as.matrix(dt.n1[,before.time]), 
                              as.matrix(dt.n2[,before.time]), shrinkage = T)
      num = round(mytest$pval, digits = 48)
      pvalue.mat[ij, paste0('Before ', OnOff), ijk] = mytest$pval
      result[ij, paste0('Before ', OnOff), ijk] = signif(mytest$pval, digits = 6)
      result[ij, paste0('BeforeT ', OnOff), ijk] = mytest$stats$statistic
      result[ij, paste0('BeforeDF1 ', OnOff), ijk] = mytest$stats$df[1]
      result[ij, paste0('BeforeDF2 ', OnOff), ijk] = mytest$stats$df[2]
      # after lightoff
      mytest = hotelling.test(as.matrix(dt.n1[,after.time]), 
                              as.matrix(dt.n2[,after.time]), shrinkage = T)
      pvalue.mat[ij, paste0('After ', OnOff), ijk] = mytest$pval
      result[ij, paste0('After ', OnOff), ijk] = mytest$pval
      result[ij, paste0('AfterT ', OnOff), ijk] = mytest$stats$statistic
      result[ij, paste0('AfterDF1 ', OnOff), ijk] = mytest$stats$df[1]
      result[ij, paste0('AfterDF2 ', OnOff), ijk] = mytest$stats$df[2]
      # cat(ij, '\t')
    }
  }
  result = as.data.frame(result)
  # FDR on Q344X
  pvalue.mat2 = apply(pvalue.mat[,,'Q344X'],2, p.adjust, method = 'BH')
  #pp = as.data.frame(pvalue.mat2)
  # FDR on Rho
  pvalue.mat2 = cbind(pvalue.mat2, apply(pvalue.mat[,,'Rho'],2, p.adjust, method = 'BH'))
  
  colnames(pvalue.mat2) = c(paste(c('Before', 'After'), OnOff, '(Q344X)'), paste(c('Before', 'After'), OnOff, '(Rho)'))
  # print out
  result[,1] = pvalue.mat2[,1]
  result[,5] = pvalue.mat2[,2]
  result[,9] = pvalue.mat2[,3]
  result[,13] = pvalue.mat2[,4]
  rm(mytest, drug.tmp, dt.n1, dt.n2)
  return(result)
}

test.diy2 = function(workingData, data.old, OnOff = 'Light-Off'){
  
  before.time = 1:t_interval
  after.time = 1:t_interval + t_interval
  
  # list of drugs for each folder, Q344X and Rho included at the end
  
  # p-value results for Mean variable
  pvalue.mat = rep(0,2)
  names(pvalue.mat) = c(paste('Before', OnOff), paste('After', OnOff))
  
  
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
  
  # after lightoff
  mytest = hotelling.test(as.matrix(dt.n1[,after.time]), 
                          as.matrix(dt.n2[,after.time]), shrinkage = T)
  pvalue.mat[paste0('After ', OnOff)] = mytest$pval
  
  colnames(pvalue.mat) = paste(c('Before', 'After'), OnOff)
  # print out
  rm(mytest, drug.tmp, dt.n1, dt.n2)
  return(pvalue.mat)
}


# 1.4 function to show p-value table via kable
pval.table = function(pvalue.mat2, OnOff = 'Light-Off'){
  kable(pvalue.mat2, caption = paste('P-values of Tests Between Drugs and Q344X/Rho Before and After', OnOff, '(Mean)'))
}
# 1.5 conduct a two-sample t-test on average mean variable of two replicates
test.consistency = function(data1 = folder1$off30data, data2 = folder2$off30data, drugs.common = drugs.common){
  drugs.common = sort(drugs.common, decreasing = F)
  data1.tmp = subset(data1, data1$genotype %in% drugs.common)
  data2.tmp = subset(data2, data2$genotype %in% drugs.common)
  
  # p-value results for two replicate comparison
  pvalues = sapply(drugs.common, function(ijk) {
    drug1 = data1.tmp[data1.tmp$genotype == ijk, c('time', 'mean_int_normalized')] 
    drug2 = data2.tmp[data2.tmp$genotype == ijk, c('time', 'mean_int_normalized')] 
    # reshape the data
    dt.n1 = func_wide(as.data.frame(drug1))
    dt.n2 = func_wide(as.data.frame(drug2))
    # do a high-dimensional two-sample t-test via testMean in HDtest package for p>=2
    # do standard two-sample t-test for p=1
    # t.test(apply(dt.n1, 2, mean), apply(dt.n2, 2, mean))$p.value
    result = testMean(dt.n1, dt.n2, method = "HD", filter = F)
    testMean(dt.n1, dt.n2, method = "HD", filter = F)$Student$p.value
  })
  statistics = sapply(drugs.common, function(ijk) {
    drug1 = data1.tmp[data1.tmp$genotype == ijk, c('time', 'mean_int_normalized')] 
    drug2 = data2.tmp[data2.tmp$genotype == ijk, c('time', 'mean_int_normalized')] 
    # reshape the data
    dt.n1 = func_wide(as.data.frame(drug1))
    dt.n2 = func_wide(as.data.frame(drug2))
    # do a high-dimensional two-sample t-test via testMean in HDtest package for p>=2
    # do standard two-sample t-test for p=1
    # t.test(apply(dt.n1, 2, mean), apply(dt.n2, 2, mean))$p.value
    testMean(dt.n1, dt.n2, method = "HD", filter = F)$Student$statistics
  })
  # pvalue and FDR
  result.df = data.frame( pval = pvalues, fdr = p.adjust(pvalues, method = 'BH'),Stat = statistics)
  # return the above table in the ascending order
  return(result.df[order(result.df$fdr, decreasing = T),])
}