CAPlots <- function(sample){
  library(ggplot2)
  library(viridis)
  library(pracma)
  source(file.path('R','multiplot.R'))
  
  model_path = file.path('output',paste('ConditionalModel',sample,sep=''))
  figure_path = file.path('figures')
  data_path = file.path('data')

  # DATA PREPARATION
  # Load in nessary data, this means average and 95% confidence intervalls of correlation coefficiencts for each pairwise interaction as a function of CWD for the CWD dependent model
  means = read.csv(file.path(model_path,'MeanCorrelationsSample.csv'),stringsAsFactors = FALSE)
  CIhigh = read.csv(file.path(model_path,'CIhighCorrelationsSample.csv'))
  CIlow =   read.csv(file.path(model_path,'CIlowCorrelationsSample.csv'))
  combinations = read.csv(file.path(model_path,'Chain1','Correlations','combinations.csv'),stringsAsFactors = FALSE,header=FALSE)
  gradient = read.csv(file.path(model_path,'Chain1','Correlations','XrGrad.csv'),stringsAsFactors = FALSE,header=FALSE)

  # Replace all non-signifiant correlations with NAs
  means[(CIhigh*CIlow)<=0] = NA  
  # Transform the gradient back to true CWD values (by deviding by the original standard deviation and adding the original mean)
  gradient = gradient[,2]*187.3861+544.73
  gradlength = length(gradient)
  
  ## -----------------------------------------------------------------------------------------------------------------------
  # COMBINE WITH ELLENBERG
  #Combine with Ellenberg
  SpeciesList = read.csv(file.path(data_path,'SpeciesList.csv'))
  Evalues = SpeciesList$Ellenberg_M
  group=Evalues
  group[Evalues >= 5 & !is.na(Evalues)] = "Sensitive"
  group[Evalues <= 3 & !is.na(Evalues)] = "Tolerant"
  group[Evalues > 3 & Evalues < 5 & !is.na(Evalues)] = "Intermediate"
  speciesnumbersDry = which(group=="Tolerant" & !is.na(group))
  speciesnumbersWet = which(group=="Sensitive" & !is.na(group))
  speciesnumbersInt = which(group=="Intermediate" & !is.na(group))

  ## -----------------------------------------------------------------------------------------------------------------------'
  # CA overall
  SpInteractionsWetAll = means[(combinations[,1] %in% speciesnumbersWet | combinations[,2] %in% speciesnumbersWet),]
  Npos = apply((sign(SpInteractionsWetAll) == 1),2,sum,na.rm=TRUE)
  Nneg = apply((sign(SpInteractionsWetAll) == -1),2,sum,na.rm=TRUE)
  TotalSig = Npos + Nneg
  RatioWetAll= (Npos - Nneg) / TotalSig
  
  SpInteractionsDryAll = means[(combinations[,1] %in% speciesnumbersDry | combinations[,2] %in% speciesnumbersDry),]
  Npos = apply((sign(SpInteractionsDryAll) == 1),2,sum,na.rm=TRUE)
  Nneg = apply((sign(SpInteractionsDryAll) == -1),2,sum,na.rm=TRUE)
  TotalSig = Npos + Nneg
  RatioDryAll = (Npos - Nneg) / TotalSig
  
  ## -----------------------------------------------------------------------------------------------------------------------
  # CA between different groups
  SpInteractionsWetWet = means[(combinations[,1] %in% speciesnumbersWet & combinations[,2] %in% speciesnumbersWet),]
  Npos = apply((sign(SpInteractionsWetWet) == 1),2,sum,na.rm=TRUE)
  Nneg = apply((sign(SpInteractionsWetWet) == -1),2,sum,na.rm=TRUE)
  TotalSig = Npos + Nneg
  RatioWetWet= (Npos - Nneg) / TotalSig
  
  SpInteractionsDryDry = means[(combinations[,1] %in% speciesnumbersDry & combinations[,2] %in% speciesnumbersDry),]
  Npos = apply((sign(SpInteractionsDryDry) == 1),2,sum,na.rm=TRUE)
  Nneg = apply((sign(SpInteractionsDryDry) == -1),2,sum,na.rm=TRUE)
  TotalSig = Npos + Nneg
  RatioDryDry = (Npos - Nneg) / TotalSig
  
  SpInteractionsDryWet = means[((combinations[,1] %in% speciesnumbersDry & combinations[,2] %in% speciesnumbersWet) | (combinations[,2] %in% speciesnumbersDry & combinations[,1] %in% speciesnumbersWet)),]
  Npos = apply((sign(SpInteractionsDryWet) == 1),2,sum,na.rm=TRUE)
  Nneg = apply((sign(SpInteractionsDryWet) == -1),2,sum,na.rm=TRUE)
  TotalSig = Npos + Nneg
  RatioDryWet = (Npos - Nneg) / TotalSig
  
  SpInteractionsDryInt = means[((combinations[,1] %in% speciesnumbersDry & combinations[,2] %in% speciesnumbersInt) | (combinations[,2] %in% speciesnumbersDry & combinations[,1] %in% speciesnumbersInt)),]
  Npos = apply((sign(SpInteractionsDryInt) == 1),2,sum,na.rm=TRUE)
  Nneg = apply((sign(SpInteractionsDryInt) == -1),2,sum,na.rm=TRUE)
  TotalSig = Npos + Nneg
  RatioDryInt = (Npos - Nneg) / TotalSig
  
  SpInteractionsWetInt = means[((combinations[,1] %in% speciesnumbersWet & combinations[,2] %in% speciesnumbersInt) | (combinations[,2] %in% speciesnumbersWet & combinations[,1] %in% speciesnumbersInt)),]
  Npos = apply((sign(SpInteractionsWetInt) == 1),2,sum,na.rm=TRUE)
  Nneg = apply((sign(SpInteractionsWetInt) == -1),2,sum,na.rm=TRUE)
  TotalSig = Npos + Nneg
  RatioWetInt = (Npos - Nneg) / TotalSig

  ## -----------------------------------------------------------------------------------------------------------------------
  # Gradient plots
  
  GradientRat = data.frame(avgint = c(RatioWetWet,RatioWetInt,
                                      RatioDryWet,RatioDryInt,RatioDryDry),
                           InteractingGroup = c(rep('DS vs DS',length(RatioWetWet)),
                                                rep('DS vs intermediate',length(RatioWetInt)),
                                                rep('DT vs. DS',length(RatioDryWet)),
                                                rep('DT vs intermediate',length(RatioDryInt)),
                                                rep('DT vs DT',length(RatioDryDry))),
                           gradient = rep(gradient,5),
                           InteractingGroupF = c(rep(3,length(RatioWetWet)),
                                                 rep(4,length(RatioWetInt)),
                                                 rep(5,length(RatioDryWet)),
                                                 rep(6,length(RatioDryInt)),
                                                 rep(7,length(RatioDryDry))))
  GradientRatOverall = data.frame(avgint = c(RatioDryAll,RatioWetAll),
                                  InteractingGroup = c(rep('Drought-tolerant (DT) vs. all',length(RatioDryAll)),
                                                       rep('Drought-sensitive (DS) vs. all',length(RatioWetAll))),
                                  gradient = rep(gradient,2),
                                  InteractingGroupF = c(rep(1,length(RatioDryAll)),
                                                        rep(2,length(RatioWetAll))))
  
  science_themea = theme(panel.grid.major = element_line(size = 0,colour="White"),
                         axis.title.y = element_text(vjust=1.25,size=8,family="Arial"),
                         axis.title.x = element_text(vjust=-0.3,size=8,family="Arial"),
                         axis.line = element_line(size = 0.7, color = "black"), 
                         axis.text.y = element_text(size=8),
                         axis.text.x = element_text(size=8), 
                         legend.title=element_blank(),
                         legend.position = c(0.62,0.9),
                         legend.direct = 'vertical',
                         legend.background = element_blank(),
                         legend.text = element_text(size=6),
                         legend.key.height = unit(1,'mm'),
                         legend.key.width = unit(5,'mm'),
                         plot.title = element_text(size=8))
  
  science_themeb = theme(panel.grid.major = element_line(size = 0,colour="White"),
                         axis.title.y = element_text(vjust=1.25,size=8,family="Arial"),
                         axis.title.x = element_text(vjust=-0.3,size=8,family="Arial"),
                         axis.line = element_line(size = 0.7, color = "black"), 
                         axis.text.y = element_text(size=8),
                         axis.text.x = element_text(size=8), 
                         legend.title=element_blank(),
                         legend.position = c(0.24,0.20),
                         legend.direct = 'vertical',
                         legend.background = element_blank(),
                         legend.text = element_text(size=6,family="Arial"),
                         legend.key.height = unit(1,'mm'),
                         legend.key.width = unit(5,'mm'),
                         plot.title = element_text(size=8))
  
  science_themec = theme(panel.grid.major = element_line(size = 0,colour="White"),
                         axis.title.y = element_text(vjust=1.25,size=8,family="Arial"),
                         axis.title.x = element_text(vjust=-0.3,size=8,family="Arial"),
                         axis.line = element_line(size = 0.7, color = "black"), 
                         axis.text.y = element_text(size=8),
                         axis.text.x = element_text(size=8), 
                         legend.title=element_blank(),
                         legend.position = c(0.44,0.15),
                         legend.direct = 'vertical',
                         legend.background = element_blank(),
                         legend.text = element_text(size=6,family="Arial"),
                         legend.key.height = unit(1,'mm'),
                         legend.key.width = unit(5,'mm'),
                         plot.title = element_text(size=8))
  
  CAPlot2 = ggplot(data=GradientRat, aes(x=gradient,y=avgint,group=as.factor(InteractingGroupF),colour=as.factor(InteractingGroupF),linetype=as.factor(InteractingGroupF),size=as.factor(InteractingGroupF))) +
    geom_line() + 
    scale_y_continuous(limits=c(-1,1)) + 
    scale_colour_manual(values=c('#F98C0AFF','#F98C0AFF','#BB3754FF','#000004FF','#000004FF'),
                        labels = c('DS vs. DS',
                                   'DS vs. rest',
                                   'DS vs. DT',
                                   'DT vs. rest',
                                   'DT vs. DT')) +
    scale_linetype_manual(values=c('dashed','dotted','dotdash','dotted','dashed'),
                          labels = c('DS vs. DS',
                                     'DS vs. rest',
                                     'DS vs. DT',
                                     'DT vs. rest',
                                     'DT vs. DT')) +
    scale_size_manual(values=c(0.5,0.5,0.5,0.5,0.5),
                      labels = c('DS vs. DS',
                                 'DS vs. rest',
                                 'DS vs. DT',
                                 'DT vs. rest',
                                 'DT vs. DT')) +
    ggtitle('c)') + 
    xlab(expression(Climatic ~ water ~ deficit ~ (mm))) +
    ylab(expression(Community ~ Association))+
    theme_bw() +
    science_themeb 
  
  CAPlot1 = ggplot(data=GradientRatOverall, aes(x=gradient,y=avgint,group=InteractingGroup,colour=InteractingGroup)) +
    geom_line(size=0.7) +
    scale_y_continuous(limits=c(-1,1)) +
    ggtitle('b)') +
    xlab(expression(Climatic ~ water ~ deficit ~ (mm))) +
    ylab(expression(Community ~ Association))+
    theme_bw() +
    science_themec +
    scale_colour_manual(values=c('#F98C0AFF','#000004FF')) 
  
  ggsave(file=file.path(figure_path,'Fig4_CA.tiff'),
         device = c('tiff'),
         plot = multiplot(CAPlot1,CAPlot2,cols=2),
         dpi = 600,
         units = c('mm'),
           width = 110, height = 55)
}