library(data.table)
library(reshape2)
library(ggplot2)

activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))

dataBDNA <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek3/united_X0.txt", sep='\t')
dataBDNA <- data.table(dataBDNA)
setnames(dataBDNA,c("Cancer","structure","isAPOBEC","sample","trgIn","cntIn","trgOut","cntOut","trgOutOld","cntOutOld","sign"))
densBDNA <- dataBDNA[structure == "gq"]
nomTCperNTin <- unique(densBDNA[isAPOBEC == 1, trgIn])
denTCperNTin <- unique(densBDNA[isAPOBEC == 0, trgIn])
TCperNTin <- nomTCperNTin / denTCperNTin
nomTCperNTout <- unique(densBDNA[isAPOBEC == 1, trgOut])
denTCperNTout <- unique(densBDNA[isAPOBEC == 0, trgOut])
TCperNTout <- nomTCperNTout / denTCperNTout


dtboxp <- data.table()
dtpctNT <- data.table()

for(canc in c("BLCA","BRCA","HNSC","CESC","LUAD","LUSC")){

  data <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek_GQ/GQ_",canc,".txt"), sep='\t')
  data <- data.table(data)
  
  nomTCperNTinStack <- unique(data[,TargStack])
  nomTCperNTinStrand <- unique(data[,TargStrnd])
  TCperNTinStack <- nomTCperNTinStack / denTCperNTin
  TCperNTinStrand <- nomTCperNTinStrand / denTCperNTin
  
  data[, densityStack := APOstack/TargStack]
  data[, densityStrand := APOstrnd/TargStrnd]
  
  dt <- merge(data,activity,by="sample",all.x = TRUE)
  
  dt[,sampleEnrich := paste0(sample,'__',round(enrichment,2))]
  sampleLevels <- unique(dt[order(enrichment),sampleEnrich])
  dt$sampleEnrich <- factor(dt$sampleEnrich,levels=sampleLevels)
  
  dt2 <- dt[,.(densityStack,densityStrand,sampleEnrich)]
  dtmelt <- melt(dt2,id.vars = c("sampleEnrich"))
    
  dtBDNA <- dataBDNA[Cancer == canc & structure == "gq" & isAPOBEC == 1]
  dtBDNA[, densityGQin := cntIn/trgIn]
  dtBDNA[, densityGQout := cntOut/trgOut]
  dtBDNA <- merge(dtBDNA,activity,by="sample",all.x = TRUE)
  dtBDNA[,sampleEnrich := paste0(sample,'__',round(enrichment,2))]
  sampleLevels <- unique(dtBDNA[order(enrichment),sampleEnrich])
  dtBDNA$sampleEnrich <- factor(dtBDNA$sampleEnrich,levels=sampleLevels)
  
  dtBDNAin <- dtBDNA[,.(sampleEnrich,"variable"="densityIn","value"=densityGQin)]
  dtBDNAout <- dtBDNA[,.(sampleEnrich,"variable"="densityOut","value"=densityGQout)]
  
  finaldt <- rbind(dtmelt,dtBDNAin,dtBDNAout)
  finaldt <- data.table(finaldt)
  finaldtNT <- copy(finaldt)
  finaldtNT[variable == "densityOut", value := value * TCperNTout]
  finaldtNT[variable == "densityIn", value := value * TCperNTin]
  finaldtNT[variable == "densityStack", value := value * TCperNTinStack]
  finaldtNT[variable == "densityStrand", value := value * TCperNTinStrand]
  finaldtNT[,variable := paste0("NT",variable)]
  
  finaldt <- rbind(finaldt,finaldtNT)
  
  ggplot(finaldt,aes(x=sampleEnrich,y=value,fill=variable)) + geom_bar(stat="identity",position="dodge") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90))
  
  ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek_GQ/pics/gq_",canc,".tiff"),units="mm",dpi=300,width=500,height=200)
  
  dt3 <- dt[,.(sample,enrichment,densityStack,densityStrand)]
  dtBDNA2 <- dtBDNA[,.(sample,densityGQout)]
  dt4 <- merge(dt3,dtBDNA2,by="sample")
  dt4[,logStackRatio := log(densityStack/densityGQout)]
  dt4[,logStrandRatio := log(densityStrand/densityGQout)]
  final2 <- dt4[,.(sample,enrichment,logStackRatio,logStrandRatio)]

  final2melt <- melt(final2,id.vars = c("sample","enrichment"))
  final2melt <- data.table(final2melt)
  
  ggplot(final2melt,aes(x=enrichment,y=value,color=variable)) + geom_point() +
    #geom_hline(yintercept=0,color=rgb(243,94,90,maxColorValue = 255)) +
    scale_color_manual(values=c(rgb(203,73,123,maxColorValue = 255), 
                                rgb(161,201,171,maxColorValue = 255))) +#, 
                                #rgb(139,197,229,maxColorValue = 255))) +
    stat_smooth(formula = y ~ log(x), method="lm",se=F) +
    #ggtitle(paste0("Structure=",s," , Cancer=",c," , isAPOBEC=",a)) +
    scale_x_continuous(breaks = seq(0, 4, by = 0.5)) +
    theme(panel.background = element_blank(),
          plot.title = element_text(size=8),
          axis.title = element_blank(),
          axis.line = element_line(color="black"),
          legend.position = "none",
          panel.grid.major = element_line(size = rel(0.5), colour='grey92'))

  ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek_GQ/pics/gq2_",canc,".tiff"),units="mm",dpi=300,width=100,height=66)

  final2melt[, cancer := canc]
  dtboxp <- rbind(dtboxp,final2melt)
  tmp <- finaldtNT[variable %in% c("NTdensityStack","NTdensityStrand")]
  tmp[, cancer := canc]
  dtpctNT <- rbind(dtpctNT,tmp)
}

dtboxpHigh <- dtboxp[enrichment >= 2.0]
dtboxpLow <- dtboxp[enrichment < 2.0]

ggplot(dtboxpHigh, aes(x=cancer,y=value,fill=variable,color=variable)) +
  #    ggplot(dtm, aes(x=cancer,y=value,fill=isAPOBEC,color=isAPOBEC)) +
  geom_boxplot(position=position_dodge(1))  +
  scale_fill_manual(name= "Clarity", values = c(rgb(118,184,174,maxColorValue = 255), rgb(229,115,161,maxColorValue = 255))) +
  scale_color_manual(name = "Clarity", values = c(rgb(119,172,166,maxColorValue = 255), rgb(203,73,123,maxColorValue = 255))) +
  geom_hline(yintercept = 0) +
  theme(panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(color="black"),
        panel.grid.major.y = element_line(size = rel(0.5), colour='grey92'))#,
       # legend.position = "none") +
  geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5),color="grey92")

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek_GQ/pics/gqHigh.tiff"),units="mm",dpi=300,width=125,height=80)

ggplot(dtboxpLow, aes(x=cancer,y=value,fill=variable,color=variable)) +
  #    ggplot(dtm, aes(x=cancer,y=value,fill=isAPOBEC,color=isAPOBEC)) +
  geom_boxplot(position=position_dodge(1))  +
  scale_fill_manual(name= "Clarity", values = c(rgb(118,184,174,maxColorValue = 255), rgb(229,115,161,maxColorValue = 255))) +
  scale_color_manual(name = "Clarity", values = c(rgb(119,172,166,maxColorValue = 255), rgb(203,73,123,maxColorValue = 255))) +
  geom_hline(yintercept = 0) +
  theme(panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(color="black"),
        panel.grid.major.y = element_line(size = rel(0.5), colour='grey92'),
        legend.position = "none") +
  geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5),color="grey92")

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek_GQ/pics/gqLow.tiff"),units="mm",dpi=300,width=125,height=80)

dtpctNT[, enrich := as.numeric(tstrsplit(sampleEnrich,'__')[[2]])]
dtpctNT[, total := sum(value), by=sampleEnrich]
dtpctNT[, pct := round(value/total,3)]
dtpctNT_copy <- copy(dtpctNT)
dtpctNT <- dtpctNT[variable == "NTdensityStack"]

p <- ggplot(dtpctNT,aes(x=enrich,y=pct)) + geom_point(color=rgb(38,120,178,maxColorValue = 255)) +
  stat_smooth(method = "lm", col = rgb(140,210,185,maxColorValue = 255),se=F,formula = y ~ log(x)) +
  #scale_color_manual(values=c(rgb(38,120,178,maxColorValue = 255),rgb(140,210,185,maxColorValue = 255))) +
  #ggtitle(paste0("Structure=",s," , Cancer=",c," , isAPOBEC=",a)) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size=8),
        axis.title = element_blank(),
        axis.line = element_line(color="black"),
        legend.position = "none",
        panel.grid.major = element_line(size = rel(0.5), colour='grey92'))

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek_GQ/pics/gqPctNT.tiff"),units="mm",width=100,height=60,dpi=300)

for(canc in c("BLCA","BRCA","HNSC","CESC","LUAD","LUSC")){
 dt4 <- dtpctNT_copy[cancer == canc]
 
 ggplot(dt4, aes(x=sampleEnrich,y=pct,fill=variable)) + geom_bar(position="fill", stat="identity") +
   scale_fill_manual(values = c(rgb(118,184,174,maxColorValue = 255), rgb(229,115,161,maxColorValue = 255))) +
  # scale_color_manual(values = c(rgb(119,172,166,maxColorValue = 255), rgb(203,73,123,maxColorValue = 255))) +
   theme(axis.text.x = element_blank(),
         panel.background = element_blank(),
         legend.position = "none",
         axis.title = element_blank(),
         axis.ticks.x = element_blank())
 
 ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek_GQ/pics/gq3_",canc,".tiff"),units="mm",width=150,height=75,dpi=300)
 
}

