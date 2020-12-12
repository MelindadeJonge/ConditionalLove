# Author: Melinda de Jonge
# Date: 12-12-2020
EffectiveSampleSize <- function(sample){
  library(mcmcse)
  library(ggplot2)
  
  model_folder = file.path('output',paste('ConditionalModel',sample,sep=''))
  
  #Effective sample size of Beta parameters
  Beta1 = read.csv(file.path(model_folder,'Chain1','Beta.csv'),header=FALSE)
  Beta2 = read.csv(file.path(model_folder,'Chain1','Beta.csv'),header=FALSE)
  Beta3 = read.csv(file.path(model_folder,'Chain1','Beta.csv'),header=FALSE)
  essBeta = data.frame(chain1 = rep(NA,ncol(Beta1)),
                       chain2 = rep(NA,ncol(Beta1)),
                       chain3 = rep(NA,ncol(Beta1)),
                       spNumber = rep(1:161,each=13))
  for(i in 1:nrow(essBeta)){
    try(essBeta$chain1[i] <- ess(Beta1[,i]),silent = F)
    try(essBeta$chain2[i] <- ess(Beta2[,i]),silent = F)
    try(essBeta$chain3[i] <- ess(Beta3[,i]),silent = F)
    print(paste('Beta',i))
  }
  essBeta$Total = rowSums(essBeta[,1:3],na.rm=T)
  write.csv(essBeta,file='output/essBeta.csv',row.names=FALSE)
  
  #Effective sample size Omega at median CWD
  SpeciesList = read.csv(file.path('data','SpeciesList.csv'))
  keep_omega = as.vector(lower.tri(matrix(nrow=nrow(SpeciesList),ncol=nrow(SpeciesList)),diag=TRUE))  #Used to remove duplicate omega parameters
  Omega1 = read.csv(file.path(model_folder,'Chain1','Omega3.csv'),header=FALSE)
  Omega1 = Omega1[,keep_omega]
  Omega2 = read.csv(file.path(model_folder,'Chain2','Omega3.csv'),header=FALSE)
  Omega2 = Omega2[,keep_omega]
  Omega3 = read.csv(file.path(model_folder,'Chain3','Omega3.csv'),header=FALSE)
  Omega3 = Omega3[,keep_omega]
  
  essOmega = data.frame(chain1 = rep(NA,ncol(Omega3)),
                            chain2 = rep(NA,ncol(Omega3)),
                            chain3 = rep(NA,ncol(Omega3)))
  for(i in 1:nrow(essOmega)){
    try(essOmega$chain1[i] <- ess(Omega1[,i]),silent = F)
    try(essOmega$chain2[i] <- ess(Omega2[,i]),silent = F)
    try(essOmega$chain3[i] <- ess(Omega3[,i]),silent = F)
    print(paste('Omega',i))
  }
  essOmega$Total = rowSums(essOmega[,1:3],na.rm=T)
  write.csv(essOmega,file='output/essOmega.csv',row.names=FALSE)

  #Histograms
  science_theme = theme(panel.grid.major = element_line(size = 0,colour="White"),
                        axis.title.y = element_text(vjust=1.25,size=8,family="Arial"),
                        axis.title.x = element_text(vjust=-0.3,size=8,family="Arial"),
                        axis.line = element_line(size = 0.7, color = "black"), 
                        axis.text.y = element_text(size=8),
                        axis.text.x = element_text(size=8), 
                        legend.title=element_blank(),
                        legend.position = 'none',
                        legend.direct = 'vertical',
                        strip.text = element_text(size=8,family="Arial"))
  
  betaHist = ggplot(data = essBeta,aes(x=Total)) +
    geom_histogram(fill = 'blue', alpha = 1/2) +
    xlab('effective sample size') +
    theme_bw() +
    xlim(c(0,nrow(Beta1))) + 
    science_theme
  betaHist
  
  ggsave(betaHist,
         file=paste0('figures/essBeta.png',sep=''),
         width = 2.5, height = 2.5)
  
  omegaHist = ggplot(data = essOmega,aes(x=Total)) +
    geom_histogram(fill = 'blue',alpha = 1/2) +
    xlab('effective sample size') +
    xlim(c(0,nrow(Omega1))) + 
    theme_bw() +
    science_theme
  
  ggsave(omegaHist,
         file=paste0('figures/essOmega.png',sep=''),
         width = 2.5, height = 2.5)
}
