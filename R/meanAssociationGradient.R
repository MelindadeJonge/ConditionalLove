meanAssociationGradient <- function(sample){
  model_path = file.path('output',paste('ConditionalModel',sample,sep=''))
  figure_path = file.path('figures')
  
  gradient = read.csv(file.path(model_path,'Chain1','Correlations','XrGrad.csv'),stringsAsFactors = FALSE,header=FALSE)
  gradlength = length(gradient[,2])
  combis = read.csv(file.path(model_path,'Chain1','Correlations','combinations.csv'),stringsAsFactors = FALSE,header=FALSE)
  ninteractions = nrow(combis)
  nspecies= max(combis)
  means=matrix(,nrow=ninteractions,ncol=gradlength)
  CIhigh=matrix(,nrow=ninteractions,ncol=gradlength)
  CIlow=matrix(,nrow=ninteractions,ncol=gradlength)
  
  #Calculate the mean an 95% CI for each interaction at each gradient intervall
  for (i in 1:gradlength){
    X1 = read.csv(file.path(model_path,'Chain1','Correlations',paste('Correlations',i,'.csv',sep='')),stringsAsFactors = FALSE,header=FALSE)
    X2 = read.csv(file.path(model_path,'Chain1','Correlations',paste('Correlations',i,'.csv',sep='')),stringsAsFactors = FALSE,header=FALSE)
    X3 = read.csv(file.path(model_path,'Chain1','Correlations',paste('Correlations',i,'.csv',sep='')),stringsAsFactors = FALSE,header=FALSE)
    X = cbind(X1,cbind(X2,X3))
    X = t(apply(X,1,sort))
    X = as.matrix(X)
    means[,i] = apply(X,1,mean)
    CIhigh[,i] = X[,round(ncol(X)/100*97.5)]
    CIlow[,i] = X[,round(ncol(X)/100*2.5)]
  }
  write.csv(means,file=file.path(model_path,'MeanCorrelationsSample.csv'),row.names=FALSE)
  write.csv(CIhigh,file=file.path(model_path,'CIhighCorrelationsSample.csv'),row.names=FALSE)
  write.csv(CIlow,file=file.path(model_path,'CIlowCorrelationsSample.csv'),row.names=FALSE)
}
