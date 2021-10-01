library(data.table)
library(ggplot2)
library(ggpubr)

savetable <- data.table()

for(cancer in c("BLCA","BRCA","HNSC","CESC","LUAD","LUSC")){

dataLoopC <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem/irloop_",cancer,".txt"), sep='\t')
dataLoopC <- data.table(dataLoopC)

dataNotLoopC <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem/irstem_",cancer,".txt"), sep='\t')
dataNotLoopC <- data.table(dataNotLoopC)

dataBDNA <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek3/united_X0.txt", sep='\t')
dataBDNA <- data.table(dataBDNA)
setnames(dataBDNA,c("Cancer","structure","isAPOBEC","sample","trgIn","cntIn","trgOut","cntOut","trgOutOld","cntOutOld","sign"))
sizeIR <- unique(dataBDNA[structure == "ir" & isAPOBEC == 0, trgIn])
sizeBDNA <- unique(dataBDNA[structure == "ir" & isAPOBEC == 0, trgOut])

dataBDNA <- dataBDNA[Cancer == cancer & structure == "ir" & isAPOBEC == 1]


activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))

resScatter <- data.table()

samples <- unique(dataLoopC$sample)
for(s in samples){
  dt1 <- dataLoopC[sample == s]
  dt2 <- dataNotLoopC[sample == s]
  dt3 <- dataBDNA[sample == s]
  enrich <- activity[sample == s, enrichment]
  
  totalTrgs1 <- dt1[pC == "T", sum(Cnt)]
  totalMuts1 <- dt1[pC == "T" & mut %in% c("G","T"),sum(Cnt)]
  
  totalTrgs1end3notT <- dt1[, sum(Cnt)]
  totalMuts1end3notT <- dt1[mut %in% c("G","T"),sum(Cnt)]
  
  totalTrgAPBloop <- dt2[LOOP==1 & pC == "T" & nuc == "C", sum(Cnt)]
  totalMutsAPBloop <- dt2[LOOP==1 & pC == "T" & nuc == "C" & mut %in% c("G","T"), sum(Cnt)]
  
  totalTrgAPBstem <- dt2[LOOP==0 & pC == "T" & nuc == "C", sum(Cnt)]
  totalMutsAPBstem <- dt2[LOOP==0 & pC == "T" & nuc == "C" & mut %in% c("G","T"), sum(Cnt)]

  totalTrgBDNA <- dt3$trgOut
  totalMutsBDNA <- dt3$cntOut
  
  totalTrgIR <- dt3$trgIn
  totalMutsIR <- dt3$cntIn
  
  FULLIRvsBDNAdensityNT <- log(((totalMuts1 + totalMutsAPBloop + totalMutsAPBstem)/sizeIR) / (totalMutsBDNA/sizeBDNA))
  LOOPSTEMvsBDNAdensityNT <- log(((totalMutsAPBloop + totalMutsAPBstem)/sizeIR) / (totalMutsBDNA/sizeBDNA))
  expectCLOOPvsBDNAdensityNT <- log(((totalTrgs1*(totalMutsAPBloop/totalTrgAPBloop+totalMutsAPBstem/totalTrgAPBstem)/2 + totalMutsAPBloop + totalMutsAPBstem)/sizeIR) / (totalMutsBDNA/sizeBDNA))
  
  resScatter <- rbind(resScatter, data.table("sample"=s, "type"="fullIR", "logratio"=FULLIRvsBDNAdensityNT, "enrichment"=enrich))
  #resScatter <- rbind(resScatter, data.table("sample"=s, "type"="LoopStem", "logratio"=LOOPSTEMvsBDNAdensityNT, "enrichment"=enrich))
  resScatter <- rbind(resScatter, data.table("sample"=s, "type"="expected", "logratio"=expectCLOOPvsBDNAdensityNT, "enrichment"=enrich))
 
  savetable <- rbind(savetable, data.table("cancer"=cancer, "sample"=s, 
                                          "totalMutsIRstemloop"=totalMuts1 + totalMutsAPBloop + totalMutsAPBstem,
                                          "totalTrgIRstemloop"=totalTrgs1+totalTrgAPBloop+totalTrgAPBstem,
                                          "totalMutsIR"=totalMutsIR,
                                          "totalTrgIR"=totalTrgIR,
                                          "sizeIR"=sizeIR,
                                          "totalMutBDNA"=totalMutsBDNA,
                                          "totalTrgBDNA"=totalTrgBDNA,
                                          "sizeBDNA"=sizeBDNA,
                                          "density2a"=log((totalMutsIR/totalTrgIR)/(totalMutsBDNA/totalTrgBDNA)),
                                          "density2aIR"=log(((totalMuts1 + totalMutsAPBloop + totalMutsAPBstem)/(totalTrgs1+totalTrgAPBloop+totalTrgAPBstem))/(totalMutsBDNA/totalTrgBDNA)),
                                          "density3a"=log(((totalMuts1 + totalMutsAPBloop + totalMutsAPBstem)/(sizeIR))/(totalMutsBDNA/sizeBDNA)),
                                          "enrichment"=enrich))
  
}

model <- lm(logratio ~ enrichment, resScatter[type=="fullIR"])
print(cancer)
print(summary(model))

ggplot(resScatter,aes(x=enrichment,y=logratio,color=type,linetype=type)) + geom_point() +
  #geom_hline(yintercept=0,color=rgb(243,94,90,maxColorValue = 255)) +
  scale_color_manual(values=c(rgb(229,115,161,maxColorValue = 255),###rgb(203,73,123,maxColorValue = 255), 
 #                             rgb(161,201,171,maxColorValue = 255), 
                              rgb(139,197,229,maxColorValue = 255))) +
  stat_smooth(formula = y ~ log(x), method="lm",se=F) +
  #ggtitle(paste0("Structure=",s," , Cancer=",c," , isAPOBEC=",a)) +
  scale_x_continuous(breaks = seq(0, 4, by = 0.5)) +
  scale_linetype_manual(values = c("dashed","solid")) +
    theme(panel.background = element_blank(),
        plot.title = element_text(size=8),
        axis.title = element_blank(),
        axis.line = element_line(color="black"),
        legend.position = "none",
        panel.grid.major = element_line(size = rel(0.5), colour='grey92'))

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig4old3/IRloopstem_",cancer,".tiff"),units="mm",dpi=300,width=100,height=66)
}

setorderv(savetable,c("cancer","enrichment"))
write.table(savetable,"/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig4old3/savetable.csv",row.names = F, quote = F, sep='\t')
