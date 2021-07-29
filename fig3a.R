library(data.table)
library(ggplot2)
library(ggpubr)

results <- data.table()

for(cancer in c("BLCA","BRCA","HNSC","CESC","LUAD","LUSC")){

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

resultsNT <- data.table()

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

  totalTrgAPBstem <- dt2[LOOP==0 & pC == "T" & nuc == "C", sum(Cnt)]
  totalMutsAPBstem <- dt2[LOOP==0 & pC == "T" & nuc == "C" & mut %in% c("G","T"), sum(Cnt)]
  
  totalTrgBDNA <- dt3$trgOut
  totalMutsBDNA <- dt3$cntOut
  
  loopCdensity <- totalMuts1 / totalTrgs1
  loopAPBdensity <- totalMutsAPBloop / totalTrgAPBloop
  stemAPBdensity <- totalMutsAPBstem / totalTrgAPBstem
  DBNAdensity <- totalMutsBDNA / totalTrgBDNA
  
  results <- rbind(results, data.table("cancer"=cancer, "sample"=s, "type"="LoopC", "density"=log(loopCdensity/DBNAdensity), "enrich"=enrich))
  results <- rbind(results, data.table("cancer"=cancer, "sample"=s, "type"="loopAPB", "density"=log(loopAPBdensity/DBNAdensity), "enrich"=enrich))
  results <- rbind(results, data.table("cancer"=cancer, "sample"=s, "type"="stemAPB", "density"=log(stemAPBdensity/DBNAdensity), "enrich"=enrich))

  resultsNT <- rbind(resultsNT, data.table("sample"=s, "enrich"=enrich, "type"="LoopC", "cnt"=totalMuts1))
  resultsNT <- rbind(resultsNT, data.table("sample"=s, "enrich"=enrich, "type"="loopAPB", "cnt"=totalMutsAPBloop))
  resultsNT <- rbind(resultsNT, data.table("sample"=s, "enrich"=enrich, "type"="stemAPB", "cnt"=totalMutsAPBstem))
  
}  
  
resultsNT[, sampleEnrich := paste0(sample,'__',round(enrich,2))]
sampleLevels <- unique(resultsNT[order(enrich),sampleEnrich])
resultsNT$sampleEnrich <- factor(resultsNT$sampleEnrich,levels=sampleLevels)

ggplot(resultsNT, aes(x=sampleEnrich,y=cnt,fill=type,labels=cnt)) + geom_bar(position="fill", stat="identity") +
  geom_text(aes(label = resultsNT$cnt),size = 3, position = position_fill(vjust = .5)) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 90))

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem/pics/ir_samples_",cancer,".tiff"),units="mm",width=300,height=150,dpi=300)

}

results <- results[!is.na(enrich)]
resultsHigh <- results[enrich >= 2.0]
resultsLow <- results[enrich < 2.0]

ggplot(resultsHigh, aes(x=cancer,y=density,fill=type,color=type)) +
  #    ggplot(dtm, aes(x=cancer,y=value,fill=isAPOBEC,color=isAPOBEC)) +
  geom_boxplot(position=position_dodge(1))  +
  scale_fill_manual(name= "Clarity", values = c("#c78ba9", "#edd585", "yellow")) +
  scale_color_manual(name = "Clarity", values = c("#a3567c", "#e9c44a", "yellow")) +
  geom_hline(yintercept = 0) +
  theme(panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(color="black"),
        panel.grid.major.y = element_line(size = rel(0.5), colour='grey92')#,
        #legend.position = "none"
        ) +
  geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5),color="grey92")

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem/pics/ir_high.tiff"),units="mm",dpi=300,width=125,height=80)


ggplot(resultsLow, aes(x=cancer,y=density,fill=type,color=type)) +
  #    ggplot(dtm, aes(x=cancer,y=value,fill=isAPOBEC,color=isAPOBEC)) +
  geom_boxplot(position=position_dodge(1))  +
  scale_fill_manual(name= "Clarity", values = c("#c78ba9", "#edd585", "yellow")) +
  scale_color_manual(name = "Clarity", values = c("#a3567c", "#e9c44a", "yellow")) +
  geom_hline(yintercept = 0) +
  theme(panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y = element_line(color="black"),
        panel.grid.major.y = element_line(size = rel(0.5), colour='grey92')#,
        #legend.position = "none"
        ) +
  geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5),color="grey92")

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem/pics/ir_low.tiff"),units="mm",dpi=300,width=125,height=80)
