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

dtNT <- unique(data[,.(structure,isAPOBEC,trgIn,trgOut)])
dtNTmerged <- merge(dtNT[isAPOBEC == 1],dtNT[isAPOBEC == 0],by="structure")
dtNTmerged[,TCperNTin := trgIn.x/trgIn.y]
dtNTmerged[,TCperNTout := trgOut.x/trgOut.y]
dtNTinout <- dtNTmerged[,.(structure,TCperNTin,TCperNTout)]

activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))

data <- merge(data,activity,by="sample",all.x = TRUE)
data <- merge(data,dtNTinout,by="structure",all.x = TRUE)

nonAnonS <- copy(data)
nonAnonS <- nonAnonS[isAPOBEC==0]
nonAnonS[,trg := trgOut]
nonAnonS[,cnt := cntOut]
nonAnonS[,type := "nonAnonS"]
nonAnonS[,density := cnt/trg]

nonAinS <- copy(data)
nonAinS <- nonAinS[isAPOBEC==0]
nonAinS[,trg := trgIn]
nonAinS[,cnt := cntIn]
nonAinS[,type := "nonAinS"]
nonAinS[,density := cnt/trg]

isAnonS <- copy(data)
isAnonS <- isAnonS[isAPOBEC==1]
isAnonS[,trg := trgOut]
isAnonS[,cnt := cntOut]
isAnonS[,type := "isAnonS"]
isAnonS[,density := cnt/trg]

isAnonSnt <- copy(isAnonS)
isAnonSnt[, type := "NTisAnonS"]
isAnonSnt[, density := density * TCperNTout]

isAinS <- copy(data)
isAinS <- isAinS[isAPOBEC==1]
isAinS[,trg := trgIn]
isAinS[,cnt := cntIn]
isAinS[,type := "isAinS"]
isAinS[,density := cnt/trg]

isAinSnt <- copy(isAinS)
isAinSnt[, type := "NTisAinS"]
isAinSnt[, density := density * TCperNTin]

dt <- rbind(nonAnonS,nonAinS,isAnonS,isAinS,isAnonSnt,isAinSnt)
dtnt <- rbind(isAnonSnt,isAinSnt)

structures <- unique(data$structure)
cancers <- unique(data$cancer)

plots <- list()

i <- 1
for(s in structures)
{
  for(c in cancers)
  {
      dt2 <- dt[cancer==c & structure == s]
      dt2[,sampleEnrich := paste0(sample,'__',round(enrichment,2))]
      sampleLevels <- unique(dt2[order(enrichment),sampleEnrich])
      dt2$sampleEnrich <- factor(dt2$sampleEnrich,levels=sampleLevels)
      
      p <- ggplot(dt2,aes(x=sampleEnrich,y=density,fill=type)) + geom_bar(stat="identity",position="dodge") +
        theme(panel.background = element_blank(),
              axis.text.x = element_blank(),
              axis.title = element_blank(),
              axis.line = element_line(color="black"),
              legend.position = "bottom") +
        scale_fill_manual(values=c(rgb(203,73,123,maxColorValue = 255), 
                                   rgb(161,201,171,maxColorValue = 255), 
                                   rgb(242,234,102,maxColorValue = 255),
                                   rgb(139,197,229,maxColorValue = 255),
                                   "grey",
                                   "cyan"))
      
        ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig2/c_nt/",s,"_",c,"_.jpg"),plot=p,width=220,height=50,units="mm",dpi=300)

        dt3 <- dtnt[cancer==c & structure == s]
        dt3[,sampleEnrich := paste0(sample,'__',round(enrichment,2))]
        sampleLevels <- unique(dt3[order(enrichment),sampleEnrich])
        dt3$sampleEnrich <- factor(dt3$sampleEnrich,levels=sampleLevels)
        
        
        p <- ggplot(dt3,aes(x=sampleEnrich,y=density,fill=type)) + geom_bar(stat="identity",position="dodge") +
          theme(panel.background = element_blank(),
                axis.text.x = element_blank(),
                axis.title = element_blank(),
                axis.line = element_line(color="black"),
                legend.position = "bottom") +
          scale_fill_manual(values=c(rgb(203,73,123,maxColorValue = 255), 
                                     rgb(161,201,171,maxColorValue = 255)))
        
        ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig2/c_nt/NT_",s,"_",c,"_.jpg"),plot=p,width=220,height=50,units="mm",dpi=300)
        
        
        
              
      plots[[i]] <- as.grob(p)  
      i <- i + 1
    
  }
}

#lay <- rbind(c(1,2),
#             c(3,4),
#             c(5,6),
#             c(7,8),
#             c(9,10))

#plts <- marrangeGrob(grobs=plots, nrow=5, ncol=1, top=NULL) 
#ggexport(plts,filename="/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek2/resultsNEWmy.pdf")

