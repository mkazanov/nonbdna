library(data.table)
library(ggplot2)
library(ggpubr)

cancer <- "BRCA"

dataLoopC <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem/irloop_",cancer,".txt"), sep='\t')
dataLoopC <- data.table(dataLoopC)

dataNotLoopC <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem/irstem_",cancer,".txt"), sep='\t')
dataNotLoopC <- data.table(dataNotLoopC)

dataBDNA <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek3/united_X0.txt", sep='\t')
dataBDNA <- data.table(dataBDNA)
setnames(dataBDNA,c("Cancer","structure","isAPOBEC","sample","trgIn","cntIn","trgOut","cntOut","trgOutOld","cntOutOld","sign"))
dataBDNA <- dataBDNA[Cancer == cancer & structure == "ir" & isAPOBEC == 1]

activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))

results <- data.table()
resScatter <- data.table()
resCount <- data.table()
resCountTrg <- data.table()

samples <- unique(dataLoopC$sample)
for(s in samples){
  dt1 <- dataLoopC[sample == s]
  dt2 <- dataNotLoopC[sample == s]
  dt3 <- dataBDNA[sample == s]
  enrich <- activity[sample == s, enrichment]
  
  totalTrgs1 <- dt1[,sum(Cnt)]
  totalMuts1 <- dt1[mut %in% c("G","T"),sum(Cnt)]
  
  totalTrgAPBloop <- dt2[LOOP==1 & pC == "T" & nuc == "C", sum(Cnt)]
  totalMutsAPBloop <- dt2[LOOP==1 & pC == "T" & nuc == "C" & mut %in% c("G","T"), sum(Cnt)]
  
  totalTrgOTHloop <- dt2[LOOP==1 & ((nuc == "C" & pC != "T") | nuc == "T"), sum(Cnt)]
  totalMutsOTHloop <- dt2[LOOP==1 & ((nuc == "C" & pC != "T") | nuc == "T") & mut != '-', sum(Cnt)]
  
  totalTrgAPBstem <- dt2[LOOP==0 & pC == "T" & nuc == "C", sum(Cnt)]
  totalMutsAPBstem <- dt2[LOOP==0 & pC == "T" & nuc == "C" & mut %in% c("G","T"), sum(Cnt)]

  totalTrgOTHstem <- dt2[LOOP==0 & ((nuc == "C" & pC != "T") | nuc == "T"), sum(Cnt)]
  totalMutsOTHstem <- dt2[LOOP==0 & ((nuc == "C" & pC != "T") | nuc == "T") & mut != '-', sum(Cnt)]
  
  totalTrgBDNA <- dt3$trgOut
  totalMutsBDNA <- dt3$cntOut
  
  loopCdensity <- totalMuts1 / totalTrgs1
  loopAPBdensity <- totalMutsAPBloop / totalTrgAPBloop
  loopOTHdensity <- totalMutsOTHloop / totalTrgOTHloop
  stemAPBdensity <- totalMutsAPBstem / totalTrgAPBstem
  stemOTHdensity <- totalMutsOTHstem / totalTrgOTHstem
  DBNAdensity <- totalMutsBDNA / totalTrgBDNA
  
  results <- rbind(results, data.table("sample"=s, "type"="1LoopC", "density"=loopCdensity))
  results <- rbind(results, data.table("sample"=s, "type"="2loopAPB", "density"=loopAPBdensity))
  results <- rbind(results, data.table("sample"=s, "type"="4bdnaAPB", "density"=DBNAdensity))
  results <- rbind(results, data.table("sample"=s, "type"="5loopOTH", "density"=loopOTHdensity))
  results <- rbind(results, data.table("sample"=s, "type"="3stemAPB", "density"=stemAPBdensity))
  results <- rbind(results, data.table("sample"=s, "type"="6stemOTH", "density"=stemOTHdensity))
  
  resScatter <- rbind(resScatter, data.table("sample"=s, "type"="1LoopC", "logratio"=loopCdensity-DBNAdensity, "enrichment"=enrich))
  resScatter <- rbind(resScatter, data.table("sample"=s, "type"="2LoopA", "logratio"=loopAPBdensity-DBNAdensity, "enrichment"=enrich))
  resScatter <- rbind(resScatter, data.table("sample"=s, "type"="3StemA", "logratio"=stemAPBdensity-DBNAdensity, "enrichment"=enrich))
 
  resCount <- rbind(resCount, data.table("sample"=s, "type"="1LoopC", "count"=totalMuts1, "enrichment"=enrich))  
  resCount <- rbind(resCount, data.table("sample"=s, "type"="2loopAPB", "count"=totalMutsAPBloop, "enrichment"=enrich))
  resCount <- rbind(resCount, data.table("sample"=s, "type"="3stemAPB", "count"=totalMutsAPBstem, "enrichment"=enrich))
  
  resCountTrg <- rbind(resCountTrg, data.table("sample"=s, "type"="1LoopC", "count"=totalTrgs1, "enrichment"=enrich))  
  resCountTrg <- rbind(resCountTrg, data.table("sample"=s, "type"="2loopAPB", "count"=totalTrgAPBloop, "enrichment"=enrich))
  resCountTrg <- rbind(resCountTrg, data.table("sample"=s, "type"="3stemAPB", "count"=totalTrgAPBstem, "enrichment"=enrich))
  
}

results <- merge(results,activity,by="sample",all.x = TRUE)

results[,sampleEnrich := paste0(sample,'__',round(enrichment,2))]
sampleLevels <- unique(results[order(enrichment),sampleEnrich])
results$sampleEnrich <- factor(results$sampleEnrich,levels=sampleLevels)

ggplot(results,aes(x=sampleEnrich,y=density,fill=type)) + geom_bar(stat="identity",position="dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem/pics/IRloopstem_",cancer,".pdf"),width=10,height=10)

ggplot(resScatter,aes(x=enrichment,y=logratio,color=type)) + geom_point() +
  #geom_hline(yintercept=0,color=rgb(243,94,90,maxColorValue = 255)) +
  scale_color_manual(values=c(rgb(203,73,123,maxColorValue = 255), 
                              rgb(161,201,171,maxColorValue = 255), 
                              rgb(139,197,229,maxColorValue = 255))) +
  stat_smooth(formula = y ~ exp(x), method="lm",se=F) +
  #ggtitle(paste0("Structure=",s," , Cancer=",c," , isAPOBEC=",a)) +
  scale_x_continuous(breaks = seq(0, 4, by = 0.5)) +
    theme(panel.background = element_blank(),
        plot.title = element_text(size=8),
        axis.title = element_blank(),
        axis.line = element_line(color="black"),
        legend.position = "none",
        panel.grid.major = element_line(size = rel(0.5), colour='grey92'))

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig3/IRloopstem_",cancer,".tiff"),units="mm",dpi=300,width=100,height=66)

resCount[,sampleEnrich := paste0(sample,'__',round(enrichment,2))]
sampleLevels <- unique(resCount[order(enrichment),sampleEnrich])
resCount$sampleEnrich <- factor(resCount$sampleEnrich,levels=sampleLevels)

ggplot(resCount,aes(x=sampleEnrich,y=count,fill=type)) + geom_bar(position="fill", stat="identity") +
  theme(axis.text.x=element_text(angle = 90, hjust = 0))

resCountTrg[,sampleEnrich := paste0(sample,'__',round(enrichment,2))]
sampleLevels <- unique(resCountTrg[order(enrichment),sampleEnrich])
resCountTrg$sampleEnrich <- factor(resCountTrg$sampleEnrich,levels=sampleLevels)

ggplot(resCountTrg,aes(x=sampleEnrich,y=count,fill=type)) + geom_bar(position="fill", stat="identity") +
  theme(axis.text.x=element_text(angle = 90, hjust = 0))
