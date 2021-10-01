library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(ggplotify)
library(viridis)

data <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek3/united_RT_X0.txt",sep='\t',header = TRUE)
data <- data.table(data)
setnames(data,c("cancer","structure","isAPOBEC","sample","trgIn","cntIn","trgOut","cntOut","trgOutOld","cntOutOld","sign","RTbin"))
data[sign < 5 | sign > 95, signBinary:=1]
data[is.na(signBinary), signBinary := 0]
data$signBinary <- as.factor(data$signBinary)
data[, ratioIn := (cntIn/trgIn)]
data[, ratioOut := (cntOut/trgOut)]
data[, ratio := log((cntIn/trgIn)/(cntOut/trgOut))]

dataAll <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek3/united_X0.txt",sep='\t',header = TRUE)
dataAll <- data.table(dataAll)
setnames(dataAll,c("cancer","structure","isAPOBEC","sample","trgIn","cntIn","trgOut","cntOut","trgOutOld","cntOutOld","sign"))
dataAll[sign < 5 | sign > 95, signBinary:=1]
dataAll[is.na(signBinary), signBinary := 0]
dataAll$signBinary <- as.factor(dataAll$signBinary)
dataAll[, ratioIn := (cntIn/trgIn)]
dataAll[, ratioOut := (cntOut/trgOut)]
dataAll[, ratio := log((cntIn/trgIn)/(cntOut/trgOut))]

activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))

data <- merge(data,activity,by="sample",all.x = TRUE)
dataAll <- merge(dataAll,activity,by="sample",all.x = TRUE)
dataAll[,RTbin := 7]

dataPlot <- rbind(data,dataAll) 

structures <- unique(data$structure)
cancers <- unique(data$cancer)

dataPlot <- dataPlot[isAPOBEC == 1]
dataPlot$RTbin <- as.factor(dataPlot$RTbin)

for(s in structures)
{
  for(c in cancers)
  {
     dt <- dataPlot[cancer==c & structure == s]
     
     #samplesInfinite <- unique(dt[is.infinite(ratio),sample])
     #dt <- dt[!(sample %in% samplesInfinite)]
     dt <- dt[!is.infinite(ratio)]
    
     p <- ggplot(dt,aes(x=enrichment,y=ratio,color=RTbin)) + geom_point(size=0.5) +
       scale_shape_manual(values=c(15, 16)) +
       geom_hline(yintercept=0,color=rgb(243,94,90,maxColorValue = 255)) +
       stat_smooth(method = "lm", se=F,formula = y ~ log(x)) +
       scale_color_manual(values=c(rgb(204,123,177,maxColorValue = 255),
                                   rgb(244,155,73,maxColorValue = 255),
                                   rgb(157,209,184,maxColorValue = 255),
                                   rgb(72,193,241,maxColorValue = 255),
                                   rgb(246,234,92,maxColorValue = 255),
                                   rgb(107,181,58,maxColorValue = 255),
                                   rgb(186,182,220,maxColorValue = 255),
                                   rgb(233,72,126,maxColorValue = 255))) +
       #ggtitle(paste0("Structure=",s," , Cancer=",c," , isAPOBEC=",a)) +
       theme(panel.background = element_blank(),
             plot.title = element_text(size=8),
             axis.title = element_blank(),
             axis.line = element_line(color="black"),
             panel.grid.major = element_line(size = rel(0.5), colour='grey92'))
     
     ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig5/bbHideInfinite/",s,"_",c,".tiff"),units="mm",width=200,height=120,dpi=300)

     # samples <- unique(dt$sample)
     # for(p in samples){
     # 
     #   dt2 <- dt[sample == p]
     #   dt2[, ratioIn := cntIn/trgIn]
     #   dt2[, ratioOut := cntOut/trgOut]
     # 
     #   dt2 <- dt2[,.(RTbin,ratioIn,ratioOut)]
     #   dt2melt <- melt(dt2,id.vars = "RTbin")
     # 
     #   ggplot(dt2melt,aes(x=RTbin,y=value,fill=variable)) + geom_bar(stat="identity",position = "dodge")
     #   ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig5/bbb/",s,"_",c,"_",p,".tiff"))
     # 
     # }
     
  }
}
     
     
          

plots <- list()

i <- 1
for(s in structures)
{
  for(c in cancers)
  {
    for(a in c(0,1))
    {
      for(rt in 0:6)
        {
        
         dt <- data[cancer==c & isAPOBEC == a & structure == s & RTbin == rt]
      
          p <- ggplot(dt,aes(x=enrichment,y=ratio)) + geom_point(color=rgb(38,120,178,maxColorValue = 255),aes(shape=signBinary)) +
            scale_shape_manual(values=c(15, 16)) +
            geom_hline(yintercept=0,color=rgb(243,94,90,maxColorValue = 255)) +
            stat_smooth(method = "lm", col = rgb(140,210,185,maxColorValue = 255),se=F,formula = y ~ log(x)) +
            #scale_color_manual(values=c(rgb(38,120,178,maxColorValue = 255),rgb(140,210,185,maxColorValue = 255))) +
            #ggtitle(paste0("Structure=",s," , Cancer=",c," , isAPOBEC=",a)) +
            theme(panel.background = element_blank(),
                  plot.title = element_text(size=8),
                  axis.title = element_blank(),
                  axis.line = element_line(color="black"),
                  legend.position = "none",
                  panel.grid.major = element_line(size = rel(0.5), colour='grey92'))

          ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig5/b/",s,"_",c,"_",a,"_",rt,".tiff"),units="mm",width=100,height=60,dpi=300)
      
      }
    }
  }
}


