library(data.table)
library(ggplot2)

cancer <- "LUSC"

data <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/LoopCheck/garvard/garv_",cancer,".txt"), sep='\t')
data <- data.table(data)

activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))

results <- data.table()

samples <- unique(data$sample)
for(s in samples){
  dt <- data[sample == s]
  
  # APOBEC no loop density
  aCG <- dt[Lng == 0 & pC == "T" & mut %in% c("G","T"), sum(Cnt)]
  aTrg <- dt[Lng == 0 & pC == "T" & mut == "-", Cnt]
  aDensity <- aCG / (aTrg + aCG)

  # APOBEC loop density
  loopA <- dt[Lng > 5 & pC == "A" & mut %in% c("G","T"), sum(Cnt)]
  loopC <- dt[Lng > 5 & pC == "C" & mut %in% c("G","T"), sum(Cnt)]
  loopG <- dt[Lng > 5 & pC == "G" & mut %in% c("G","T"), sum(Cnt)]
  loopT <- dt[Lng > 5 & pC == "T" & mut %in% c("G","T"), sum(Cnt)]

  loopAtrg <- dt[Lng > 5 & pC == "A" & mut == "-", sum(Cnt)]
  loopCtrg <- dt[Lng > 5 & pC == "C" & mut == "-", sum(Cnt)]
  loopGtrg <- dt[Lng > 5 & pC == "G" & mut == "-", sum(Cnt)]
  loopTtrg <- dt[Lng > 5 & pC == "T" & mut == "-", sum(Cnt)]
  
  loopAdensity <- loopA / (loopA + loopAtrg)
  loopCdensity <- loopC / (loopC + loopCtrg)
  loopGdensity <- loopG / (loopG + loopGtrg)
  loopTdensity <- loopT / (loopT + loopTtrg)
  
  results <- rbind(results, data.table("sample"=s, "type"="LoopA", "density"=loopAdensity))
  results <- rbind(results, data.table("sample"=s, "type"="LoopC", "density"=loopCdensity))
  results <- rbind(results, data.table("sample"=s, "type"="LoopG", "density"=loopGdensity))
  results <- rbind(results, data.table("sample"=s, "type"="LoopT", "density"=loopTdensity))
  results <- rbind(results, data.table("sample"=s, "type"="noloop", "density"=aDensity))

}

results <- merge(results,activity,by="sample",all.x = TRUE)

results[,sampleEnrich := paste0(sample,'__',round(enrichment,2))]
sampleLevels <- unique(results[order(enrichment),sampleEnrich])
results$sampleEnrich <- factor(results$sampleEnrich,levels=sampleLevels)

ggplot(results,aes(x=sampleEnrich,y=density,fill=type)) + geom_bar(stat="identity",position="dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/LoopCheck/garvard/garv_",cancer,".pdf"),width=10,height=10)
