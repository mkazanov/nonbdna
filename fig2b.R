library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(ggplotify)
library(cowplot)
library(reshape2)

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
#      dt <- rbind(data[isAPOBEC == 1 & structure == s & enrichment > 2.0],
#                  data[isAPOBEC == 0 & structure == s & enrichment <= 2.0])
        dt <- rbind(data[isAPOBEC == 1 & structure == s],
                    data[isAPOBEC == 1 & structure == s])
        dt[enrichment > 2.0, isEnriched := 1]
        dt[enrichment <= 2.0, isEnriched := 0]
        
        dt$isAPOBEC <- as.factor(dt$isAPOBEC)
        dt$isEnriched <- as.factor(dt$isEnriched)
      
      dts <- dt[,.(sample,cancer,isEnriched,ratio)]
      dtm <- melt(dts, id.vars = c("sample","cancer","isEnriched"))
      #dts <- dt[,.(sample,cancer,isAPOBEC,ratio)]
      #dtm <- melt(dts, id.vars = c("sample","cancer","isAPOBEC"))
      
      ggplot(dtm, aes(x=cancer,y=value,fill=isEnriched,color=isEnriched)) +
  #    ggplot(dtm, aes(x=cancer,y=value,fill=isAPOBEC,color=isAPOBEC)) +
      geom_boxplot(position=position_dodge(1))  +
        scale_fill_manual(name= "Clarity", values = c("#c78ba9", "#edd585")) +
        scale_color_manual(name = "Clarity", values = c("#a3567c", "#e9c44a")) +
        geom_hline(yintercept = 0) +
        theme(panel.background = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.line.x = element_blank(),
              axis.ticks.x=element_blank(),
              axis.line.y = element_line(color="black"),
              panel.grid.major.y = element_line(size = rel(0.5), colour='grey92'),
              legend.position = "none") +
        geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5),color="grey92")
      
      
      
      ggsave(filename=paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig2/b/fig2b_",s,".tiff"),units="mm",width=100,height=50,dpi=300)
    
}


