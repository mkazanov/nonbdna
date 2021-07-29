library(data.table)
library(ggplot2)

cancer <- "LUSC"

dataLoopC <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem/irloop_",cancer,".txt"), sep='\t')
dataLoopC <- data.table(dataLoopC)

dataNotLoopC <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem/irstem_",cancer,".txt"), sep='\t')
dataNotLoopC <- data.table(dataNotLoopC)

activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))

results <- data.table()

samples <- unique(dataLoopC$sample)
for(s in samples){
  dt1 <- dataLoopC[sample == s]
  dt2 <- dataNotLoopC[sample == s]
  
  totalTrgs1 <- dt1[,sum(Cnt)]
  totalMuts1 <- dt1[mut %in% c("G","T"),sum(Cnt)]
  
  totalTrg2LOOP <- dt2[LOOP==1 & nuc %in% c("C","T"), sum(Cnt)]
  totalMuts2LOOP <- dt2[LOOP==1 & nuc %in% c("C","T") & mut != '-', sum(Cnt)]
  
  totalTrg2STEM <- dt2[LOOP==0 & nuc %in% c("C","T"), sum(Cnt)]
  totalMuts2STEM <- dt2[LOOP==0 & nuc %in% c("C","T") & mut != '-', sum(Cnt)]

  loopCdensity <- totalMuts1 / totalTrgs1
  loopNotCdensity <- totalMuts2LOOP / totalTrg2LOOP
  stemDensity <- totalMuts2STEM / totalTrg2STEM
    
  results <- rbind(results, data.table("sample"=s, "type"="LoopC", "density"=loopCdensity))
  results <- rbind(results, data.table("sample"=s, "type"="loopNotC", "density"=loopNotCdensity))
  results <- rbind(results, data.table("sample"=s, "type"="stem", "density"=stemDensity))
}

results <- merge(results,activity,by="sample",all.x = TRUE)

results[,sampleEnrich := paste0(sample,'__',round(enrichment,2))]
sampleLevels <- unique(results[order(enrichment),sampleEnrich])
results$sampleEnrich <- factor(results$sampleEnrich,levels=sampleLevels)

ggplot(results,aes(x=sampleEnrich,y=density,fill=type)) + geom_bar(stat="identity",position="dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem/pics/IRloopstem_",cancer,".pdf"),width=10,height=10)
