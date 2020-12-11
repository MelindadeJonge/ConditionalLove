# Author: Melinda de Jonge
# Date: 14-11-2018
checkResiduals <- function(sample){
  library(DHARMa)
  library(ggplot2)
  library(viridis)
  library(pracma)
  source('R/multiplot.R')
  
  ## -------- Ellenberg -----------------------------------------------------------------------------------------------------------------
  model_folder = file.path('output',paste('ConditionalModel',sample,sep=''))
  SpeciesList = read.csv(file.path('data','SpeciesList.csv'))
  Evalues = SpeciesList$Ellenberg_M
  group=Evalues
  group[Evalues >= 5 & !is.na(Evalues)] = "Sensitive"
  group[Evalues <= 3 & !is.na(Evalues)] = "Tolerant"
  group[Evalues > 3 & Evalues < 5 & !is.na(Evalues)] = "Intermediate"
  
  ## -------- Scaled Residuals -----------------------------------------------------------------------------------------------------------------
  Y = read.csv(file.path('data','Occurrences.csv'),header = FALSE,stringsAsFactors=FALSE)
  X = read.csv(file.path('data','Covariates.csv'),stringsAsFactors = FALSE)

  plot_combinations = read.csv(file.path('data','UniquePlotCombinations.csv'),header = FALSE)
  Y = Y[(Y[,1] %in% plot_combinations[,sample]),2:ncol(Y)]
  X = X[(X[,1] %in% plot_combinations[,sample]),2:ncol(X)]
  ny=nrow(Y)
  scaledRes = matrix(NA,nrow=nrow(Y),ncol=ncol(Y))
  for (i in 1:ncol(Y)){
    Ypred = as.matrix(read.csv(file.path(model_folder,'Predictions',paste('Species',i,'.csv',sep='')),header=FALSE))
    print(i)
    MedPred = apply(Ypred,1,median)
    D = createDHARMa(simulatedResponse = Ypred[(ny+1):(2*ny),],
                     observedResponse = Y[,i], 
                     fittedPredictedResponse = MedPred[(ny+1):(2*ny)],
                     integerResponse=TRUE)
    scaledRes[,i] = D$scaledResiduals
  }
  scaledRes = data.frame(scaledRes, X$CWD)
  colnames(scaledRes) = c(as.character(SpeciesList$SpeciesName),'CWD')
  write.csv(scaledRes,file=file.path(model_folder,'scaledResiduals.csv'),row.names=FALSE)
  
  ## ---------Specieswise residual correlation CWD --------------------------------------------------------------------------------------------------------------
  cres = cor(scaledRes$CWD,scaledRes[,1:161])
  cres = data.frame(correlations=cres[1,],
                    group=group)
  for(i in 1:161){
    fit = lm(scaledRes[,i]~scaledRes$CWD)
    cres$Beta[i] = fit$coefficients[2]
  }

  # -------- Plotting Scaled residuals specieswise -----------------------------------------------------------------------------------------------------------------
  science_theme = theme(panel.grid.major = element_line(size = 0,colour="White"),
                        axis.title.y = element_text(vjust=1.25,size=8,family="Arial"),
                        axis.title.x = element_text(vjust=-0.3,size=8,family="Arial"),
                        axis.line = element_line(size = 0.7, color = "black"),
                        axis.text.y = element_text(size=8),
                        axis.text.x = element_text(size=8),
                        legend.title=element_blank(),
                        legend.position = 'right',
                        legend.direct = 'vertical',
                        legend.key.size = unit(0.1,'cm'),
                        strip.text = element_text(size=12,family="Arial"),
                        plot.title = element_text(size=8))
  tryCatch({dir.create(file.path('figures','Residuals'))})
  for(i in 1:nrow(cres)){
    X = data.frame(R = scaledRes[,i],
                   CWD = scaledRes$CWD)
    gg = ggplot(data = X,aes(x=CWD,y=R)) +
      geom_point(size=0.25) +
      geom_smooth(method = lm,size=0.4) +
      xlab('Climate Water Deficit (mm)') +
      ylab('scaled residuals') +
      geom_hline(yintercept=0.5,linetype='dashed',size=0.25)+
      ggtitle(paste(colnames(scaledRes)[i], ' (slope: ', round(cres$Beta[i],6), ')',sep='')) +
      theme_bw() +
      science_theme
    ggsave(gg,
           file=file.path('figures','Residuals',paste('Species',i,'.png',sep='')),
           width = 2.5, height = 2.5)
  }
  cres = cres[!is.na(cres$group),]
  corrhist = ggplot(data = cres,aes(x=Beta,fill=group)) +
    geom_histogram(data=subset(cres,group == 'Sensitive'),alpha=1/2,binwidth=0.00002)+
    geom_histogram(data=subset(cres,group == 'Tolerant'),alpha=1/2,binwidth=0.00002)+
    geom_histogram(data=subset(cres,group == 'Intermediate'),alpha=1/2,binwidth=0.00002)+
    xlab('Slope') +
    scale_fill_viridis(discrete = TRUE) +
    xlim(c(-0.0005,0.0005)) +
    theme_bw() +
    science_theme
  corrhist
  
  ggsave(corrhist,
         file=file.path('figures','SlopesResiduals.png'),
         width = 5, height = 2.5)
}
