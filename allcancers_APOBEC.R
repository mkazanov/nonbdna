library(data.table)
library(ggplot2)

ROOTDIR <- "/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/all_cancers_APOBEC/"

results <- data.table()
for(NBtype in c('apr','dr','gq','ir','mr','str','z')){
  data <- read.csv(paste0(ROOTDIR,NBtype,"_mu.txt"),sep='\t')
  data <- data.table(data)
  data[, CancerType := substr(Cancer,1,4)]
  #data[, CancerType := Cancer]
  data[, logDensity := log10((APOin/TCxIN)/(APOout/TCxOUT))]
  data <- data[order(CancerType,logDensity)]
  data[, ind := 1:.N, by=CancerType]
  data[, groupN := .N, by=CancerType]
  data[, xx := (ind-1)/(groupN-1)]
  
  dt <- data[is.finite(logDensity),.(meanLogDensity = mean(logDensity)), by=Cancer]
  meanOfMeans <- dt[,mean(meanLogDensity)]  
  p <- wilcox.test(dt$meanLogDensity, mu = 0, alternative = "two.sided")
  
  results <- rbind(results,data.table("NBtype"=NBtype,
                                      "logMeanOfMeans"=meanOfMeans,
                                      "foldEnrichment"=10^meanOfMeans,
                                      "p-value"=p$p.value))
  
  ggplot(data, aes(x=xx,y=logDensity)) + geom_point(size=0.1) + facet_grid(~CancerType) + theme_bw() + geom_hline(yintercept=0,color="red") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 10),
          axis.title.y = element_blank(),
          panel.spacing = unit(0, "lines"))
  
  ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/figS1/panel_",NBtype,".tiff"),units="mm",dpi=300,width=500,height=50)
}

write.table(results,"/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/figS1/figS1_stat.txt",row.names = F, quote = F, sep='\t')
