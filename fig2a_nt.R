library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(ggplotify)
library(viridis)

data <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek3/united_X0.txt",sep='\t',header = TRUE)
data <- data.table(data)
setnames(data,c("cancer","structure","isAPOBEC","sample","trgIn","cntIn","trgOut","cntOut","trgOutOld","cntOutOld","sign"))
data[sign < 5 | sign > 95, signBinary:=1]
data[is.na(signBinary), signBinary := 0]
data$signBinary <- as.factor(data$signBinary)
data[, ratio := log((cntIn/trgIn)/(cntOut/trgOut))]

activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))

data <- merge(data,activity,by="sample",all.x = TRUE)

structures <- unique(data$structure)
cancers <- unique(data$cancer)

plots <- list()

i <- 1
for(s in structures)
{
  for(c in cancers)
  {
    for(a in c(0,1))
    {
      dt <- data[cancer==c & isAPOBEC == a & structure == s]
      
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

      ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig2/a/",s,"_",c,"_",a,".tiff"),units="mm",width=100,height=60,dpi=300)
      
    }
  }
}


