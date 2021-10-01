library(data.table)
library(ggplot2)
library(ggpubr)

results <- data.table()

dataBDNA <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/all_cancers_APOBEC/ir_mu.txt", sep='\t')
dataBDNA <- data.table(dataBDNA)



for(cancer in c(#"BLCA-US",
                "BOCA-UK",
                #"BRCA-EU",
                #"BRCA-UK",
                #"BRCA-US",
                "BTCA-SG",
                #"CESC-US",
                "CLLE-ES",
                "CMDI-UK",
                "COAD-US",
                "DLBC-US",
                "EOPC-DE",
                "ESAD-UK",
                "GACA-CN",
                "GBM-US",
                #"HNSC-US",
                "KICH-US",
                "KIRC-US",
                "KIRP-US",
                "LAML-KR",
                "LGG-US",
                "LICA-FR",
                "LIHC-US",
                "LINC-JP",
                "LIRI-JP",
                #"LUAD-US",
                #"LUSC-US",
                "MALY-DE",
                "MELA-AU",
                "ORCA-IN",
                "OV-AU",
                "OV-US",
                "PACA-AU",
                "PACA-CA",
                "PAEN-AU",
                "PAEN-IT",
                "PBCA-DE",
                "PRAD-CA",
                "PRAD-UK",
                "PRAD-US",
                "READ-US",
                "RECA-EU",
                "SARC-US",
                "SKCM-US",
                "STAD-US",
                "THCA-US",
                "UCEC-US")){

  dataLoopC <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem_allcancers/irloop_",cancer,".txt"), sep='\t')
  dataLoopC <- data.table(dataLoopC)
  
  dataNotLoopC <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem_allcancers//irstem_",cancer,".txt"), sep='\t')
  dataNotLoopC <- data.table(dataNotLoopC)
  
  samples <- unique(dataLoopC$sample)
  for(s in samples){
    dt1 <- dataLoopC[sample == s]
    dt2 <- dataNotLoopC[sample == s]
    dt3 <- dataBDNA[sample == s]
    
    totalTrgs1 <- dt1[pC == "T",sum(Cnt)]
    totalMuts1 <- dt1[pC == "T" & mut %in% c("G","T"),sum(Cnt)]
    
    totalTrgAPBloop <- dt2[LOOP==1 & pC == "T" & nuc == "C", sum(Cnt)]
    totalMutsAPBloop <- dt2[LOOP==1 & pC == "T" & nuc == "C" & mut %in% c("G","T"), sum(Cnt)]
    
    totalTrgAPBstem <- dt2[LOOP==0 & pC == "T" & nuc == "C", sum(Cnt)]
    totalMutsAPBstem <- dt2[LOOP==0 & pC == "T" & nuc == "C" & mut %in% c("G","T"), sum(Cnt)]
    
    totalTrgBDNA <- dt3$TCxOUT
    totalMutsBDNA <- dt3$APOout
    
    loopCdensity <- totalMuts1 / totalTrgs1
    loopAPBdensity <- totalMutsAPBloop / totalTrgAPBloop
    stemAPBdensity <- totalMutsAPBstem / totalTrgAPBstem
    DBNAdensity <- totalMutsBDNA / totalTrgBDNA

    results <- rbind(results, data.table("cancer"=cancer, "sample"=s, "type"="1LoopC", "density"=log(loopCdensity/DBNAdensity)))
    results <- rbind(results, data.table("cancer"=cancer, "sample"=s, "type"="2loopAPB", "density"=log(loopAPBdensity/DBNAdensity)))
    results <- rbind(results, data.table("cancer"=cancer, "sample"=s, "type"="3stemAPB", "density"=log(stemAPBdensity/DBNAdensity)))
    
  }  
}
  
  ggplot(results, aes(x=cancer,y=density,fill=type,color=type)) +
    #    ggplot(dtm, aes(x=cancer,y=value,fill=isAPOBEC,color=isAPOBEC)) +
    geom_boxplot(position=position_dodge(1))  +
    scale_fill_manual(name= "Clarity", values = c(rgb(253,190,147,maxColorValue = 255), 
                                                  rgb(233,148,151,maxColorValue = 255), 
                                                  rgb(146,187,216,maxColorValue = 255))) +
    scale_color_manual(name = "Clarity", values = c(rgb(253,127,40,maxColorValue = 255), 
                                                    rgb(212,42,47,maxColorValue = 255), 
                                                    rgb(38,120,178,maxColorValue = 255))) +
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
    geom_vline(xintercept=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,
                            20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5,30.5,31.5,32.5,33.5,34.5,35.5,36.5,37.5,
                            38.5,39.5,40.5,41.5,42.5,43.5,44.5,45.5,46.5,47.5),color="grey92")
  
  ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/IRloopstem_allcancers/0_pics/ir_allcancers.tiff"),units="mm",dpi=300,width=1025,height=80)
  
