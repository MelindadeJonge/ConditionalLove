# Author: Melinda de Jonge
# Date: 14-11-2018
associationChange <- function(sample){
  model_path = file.path('output',paste('ConditionalModel',sample,sep=''))
  figure_path = file.path('figures')
  
  #Combine associations at 5th, 25th, 75th and 95th percentile for all 3 chains
  Quantile1_1 = read.csv(file.path(model_path,'Chain1','Correlations','QuantileCorrelations1.csv'),header=FALSE)
  Quantile1_2 = read.csv(file.path(model_path,'Chain2','Correlations','QuantileCorrelations1.csv'),header=FALSE)
  Quantile1_3 = read.csv(file.path(model_path,'Chain3','Correlations','QuantileCorrelations1.csv'),header=FALSE)
  Quantile2_1 = read.csv(file.path(model_path,'Chain1','Correlations','QuantileCorrelations2.csv'),header=FALSE)
  Quantile2_2 = read.csv(file.path(model_path,'Chain2','Correlations','QuantileCorrelations2.csv'),header=FALSE)
  Quantile2_3 = read.csv(file.path(model_path,'Chain3','Correlations','QuantileCorrelations2.csv'),header=FALSE)
  Quantile4_1 = read.csv(file.path(model_path,'Chain1','Correlations','QuantileCorrelations4.csv'),header=FALSE)
  Quantile4_2 = read.csv(file.path(model_path,'Chain2','Correlations','QuantileCorrelations4.csv'),header=FALSE)
  Quantile4_3 = read.csv(file.path(model_path,'Chain3','Correlations','QuantileCorrelations4.csv'),header=FALSE)
  Quantile5_1 = read.csv(file.path(model_path,'Chain1','Correlations','QuantileCorrelations5.csv'),header=FALSE)
  Quantile5_2 = read.csv(file.path(model_path,'Chain2','Correlations','QuantileCorrelations5.csv'),header=FALSE)
  Quantile5_3 = read.csv(file.path(model_path,'Chain3','Correlations','QuantileCorrelations5.csv'),header=FALSE)
  Quantile1 = cbind2(Quantile1_1,cbind2(Quantile1_2,Quantile1_3))
  Quantile2 = cbind2(Quantile2_1,cbind2(Quantile2_2,Quantile2_3))
  Quantile4 = cbind2(Quantile4_1,cbind2(Quantile4_2,Quantile4_3))
  Quantile5 = cbind2(Quantile5_1,cbind2(Quantile5_2,Quantile5_3))
  
  # Check if associations at 5th and 95th are significantly different
  Pairwise25th75th = data.frame(matrix(,ncol=5,nrow=nrow(Quantile2)))
  colnames(Pairwise25th75th) = c('Mean difference','#N medium<low','#N medium>low','Probability medium>low')
  for (i in 1:nrow(Pairwise25th75th)){
    Pairwise25th75th[i,1] = mean(as.numeric(Quantile4[i,]))-mean(as.numeric(Quantile2[i,]))
    diff = as.numeric(Quantile4[i,]) > as.numeric(Quantile2[i,])
    Pairwise25th75th[i,2] = sum(diff==FALSE)
    Pairwise25th75th[i,4] = sum(diff==TRUE)
    Pairwise25th75th[i,5] = 1 - pnorm(0,mean=(mean(as.numeric(Quantile4[i,]))-mean(as.numeric(Quantile2[i,]))),sd=sqrt((sd(as.numeric(Quantile2[i,]))^2+sd(as.numeric(Quantile4[i,]))^2)))
  }
  write.csv(Pairwise25th75th,file=file.path(model_path,'associationChange25th75th.csv'),row.names=FALSE)
  
  # Check if associations at 5th and 95th are significantly different
  Pairwise595 = data.frame(matrix(,ncol=5,nrow=nrow(Quantile1)))
  colnames(Pairwise595) = c('Mean difference','#N medium<low','#N medium>low','Probability medium>low')
  for (i in 1:nrow(Pairwise595)){
    Pairwise595[i,1] = mean(as.numeric(Quantile5[i,]))-mean(as.numeric(Quantile1[i,]))
    diff = as.numeric(Quantile5[i,]) > as.numeric(Quantile1[i,])
    Pairwise595[i,2] = sum(diff==FALSE)
    Pairwise595[i,3] = sum(diff==TRUE)
    Pairwise595[i,4] = 1 - pnorm(0,mean=(mean(as.numeric(Quantile5[i,]))-mean(as.numeric(Quantile1[i,]))),sd=sqrt((sd(as.numeric(Quantile1[i,]))^2+sd(as.numeric(Quantile5[i,]))^2)))
  }
  write.csv(Pairwise595,file=file.path(model_path,'associationChange5th95th.csv'),row.names=FALSE)
}