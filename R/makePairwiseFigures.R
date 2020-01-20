# Author: Melinda de Jonge
# Date: 14-11-2018
makePairwiseFigures <- function(sample){
  library(ggplot2)
  library(viridis)
  library(pracma)
  model_path = file.path('output',paste('ConditionalModel',sample,sep=''))
  figure_path = file.path('figures')
  data_path = file.path('data')
  Pairwise25th75th = read.csv(file.path(model_path,'associationChange25th75th.csv'),stringsAsFactors=FALSE)
  
  means = read.csv(file.path(model_path,'MeanCorrelationsSample.csv'),stringsAsFactors = FALSE)
  combinations = read.csv(file.path(model_path,'Chain1','Correlations','combinations.csv'),stringsAsFactors = FALSE,header=FALSE)
  gradient = read.csv(file.path(model_path,'Chain1','Correlations','XrGrad.csv'),stringsAsFactors = FALSE,header=FALSE)
  gradient = gradient[,2]*187.3861+544.73
  
  # Get ellenberg values
  SpeciesList = read.csv(file.path(data_path,'SpeciesList.csv'))
  Evalues = SpeciesList$Ellenberg_M
  group=Evalues
  group[Evalues >= 5 & !is.na(Evalues)] = "Drought intolerant"
  group[Evalues <= 3 & !is.na(Evalues)] = "Drought tolerant"
  speciesnumbersDry = which(group=="Drought tolerant" & !is.na(group))
  speciesnumbersWet = which(group=="Drought intolerant" & !is.na(group))

  ## -----------------------------------------------------------------------------------------------------------------------
  # MAKE GROUPS
  Data = data.frame(association=matrix(unlist(means),ncol=1),
                    gradient = rep(gradient,each=(nrow(means))),
                    sp1 = rep(combinations[,1],length(gradient)),
                    sp2 = rep(combinations[,2],length(gradient)),
                    Pr1 = rep(Pairwise25th75th$Probability.medium.low/3,length(gradient)))
  Data$type = NA
  Data$type[Data$Pr1>0.95] = 'Positive'
  Data$type[Data$Pr1<0.05] = 'Negative'
  Data$type[is.na(Data$type)] = 'No response'
  
  DataDS = Data[(Data$sp1 %in% speciesnumbersWet | Data$sp2 %in% speciesnumbersWet),]
  DataDT = Data[(Data$sp1 %in% speciesnumbersDry | Data$sp2 %in% speciesnumbersDry),]
  DataDS$group = 'Drought-sensitive'
  DataDT$group = 'Drought-tolerant'
  Data=rbind(DataDS,DataDT)
  Data$combo = paste(as.character(Data$sp1),'vs',as.character(Data$sp2),sep="")
  
  ## -----------------------------------------------------------------------------------------------------------------------
  # Plot lines
  #PLOTTING
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
  
  Data = transform(Data,type=factor(type,levels=c('Positive',
                                                  'Negative',
                                                  'No response')))
  
  g1 = ggplot(data=Data,aes(x=gradient,y=association,group=combo))+
    geom_line(alpha=0.1)+
    facet_grid(type~group)+
    xlab('Climate Water Deficit (mm)') +
    ylab('Species association') +
    theme_bw() +
    science_theme 
  
  ggsave(file=file.path(figure_path,'Fig3_pairwise.tiff'),
         device = c('tiff'),
         plot = g1,
         dpi = 600,
         units = c('mm'),
         width = 82, height = 123)
}
