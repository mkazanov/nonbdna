library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(ggplotify)

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
      
      p <- ggplot(dt,aes(x=enrichment,y=ratio)) + geom_point(aes(color=signBinary)) +
        ggtitle(paste0("Structure=",s," , Cancer=",c," , isAPOBEC=",a)) +
        theme_light() +
        theme(plot.title = element_text(size=8),
              axis.title.x = element_blank())
      plots[[i]] <- as.grob(p)  
      i <- i + 1
    }
  }
}

lay <- rbind(c(1,2),
             c(3,4),
             c(5,6),
             c(7,8),
             c(9,10))

plts <- marrangeGrob(grobs=plots, nrow=6, ncol=2, top=NULL, layout_matrix = lay) 
ggexport(plts,filename="/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig3/resultsNEW.pdf")

