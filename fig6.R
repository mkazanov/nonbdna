library(ggplot2)
library(data.table)

std <- function(x) sd(x)/sqrt(length(x))

dtsize <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek3/united_X0.txt", sep='\t')
dtsize <- data.table(dtsize)
dtsize <- dtsize[APO == 0 & struID == "z",.(szIN,szOUall)]
dts <- unique(dtsize[,])

for(cancer in c("MELA-AU","SKCM-US")){
  data <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/ZDNA/Z_",cancer,".txt"),sep='\t')
  data <- data.table(data)
  dt <- unique(data[,.(tTCin,tTCout,tCTin,tCTout,tCCin,tCCout)])
  dt2 <- cbind(dt,dts)
  dt2[, tTCinRatio := tTCin / szIN]
  dt2[, tTCoutRatio := tTCout / szOUall]
  dt2[, tCTinRatio := tCTin / szIN]
  dt2[, tCToutRatio := tCTout / szOUall]
  dt2[, tCCinRatio := tCCin / szIN]
  dt2[, tCCoutRatio := tCCout / szOUall]
  dt2 <- dt2[,.(tTCinRatio,tTCoutRatio,tCCinRatio,tCCoutRatio)]
  dt2melt <- melt(dt2)
  dt2melt[,motif:=substr(variable,2,3)]
  dt2melt[,region:=substr(variable,4,5)]
  
  # update due to Gena's specific counting of targets
  dt2melt[motif=='TC' & region == "in", value := 0.007927096]
  dt2melt[motif=='TC' & region == "ou", value := 0.108816240]
  
  ggplot(dt2melt,aes(x=motif,y=value,fill=region,color=region)) + geom_bar(stat="identity",position="dodge") +
    theme(panel.background = element_blank(),
          axis.title = element_blank(),
          axis.line = element_line(color="black"),
          axis.text.x = element_blank(),
          legend.position = "none") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_fill_manual(values = c(rgb(233,148,151,maxColorValue = 255), 
                                                  rgb(146,187,216,maxColorValue = 255))) +
    scale_color_manual(values = c(rgb(212,42,47,maxColorValue = 255), 
                                                    rgb(38,120,178,maxColorValue = 255)))
  
  ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig6/fig6b.tiff"),dpi=300,units="mm",width=70,height=50)
  
    
  # dm <- data[,.(mTCin=sum(mTCin),
  #               mTCout=sum(mTCout),
  #               mCTin=sum(mCTin),
  #               mCTout=sum(mCTout),
  #               mcCin=sum(mcCin),
  #               mcCout=sum(mcCout),
  #               mCcin=sum(mCcin),
  #               mCcout=sum(mCcout))]
  # 
  # dm2 <- cbind(dm,dt)
  # dm2[, mTCinRatio := mTCin / tTCin]
  # dm2[, mTCoutRatio := mTCout / tTCout]
  # dm2[, mCTinRatio := mCTin / tCTin]
  # dm2[, mCToutRatio := mCTout / tCTout]
  # dm2[, mcCinRatio := mcCin / tCCin]
  # dm2[, mcCoutRatio := mcCout / tCCout]
  # dm2[, mCcinRatio := mCcin / tCCin]
  # dm2[, mCcoutRatio := mCcout / tCCout]
  
  data[, dTCin := mTCin / tTCin]
  data[, dTCout := mTCout / tTCout]
  data[, dCCin := (mcCin+mCcin) / tCCin]
  data[, dCCout := (mcCout+mCcout) / tCCout ]
  data[, dCTin := mCTin / tCTin]
  data[, dCTout := mCTout / tCTout ]
  data[, rTCin := dTCin/(dTCin+dTCout)] 
  data[, rTCout := dTCout/(dTCin+dTCout)]
  data[, rCCin := dCCin/(dCCin+dCCout)] 
  data[, rCCout := dCCout/(dCCin+dCCout)]
  data[, rCTin := dCTin/(dCTin+dCTout)] 
  data[, rCTout := dCTout/(dCTin+dCTout)]
  
  dm <- data[(mcCin + mCcin) >= 5,.(rTCin,rTCout,rCCin,rCCout)]
  
  print(wilcox.test(dm$rTCin,dm$rTCout))
  print(wilcox.test(dm$rCCin,dm$rCCout))
  
  dmMeans <- dm[,.("rTCin"=mean(rTCin),"rTCout"=mean(rTCout),"rCCin"=mean(rCCin),"rCCout"=mean(rCCout))]
  dmSE <- dm[,.("rTCin"=std(rTCin),"rTCout"=std(rTCout),"rCCin"=std(rCCin),"rCCout"=std(rCCout))]
  dmSTD <- dm[,.("rTCin"=sd(rTCin),"rTCout"=sd(rTCout),"rCCin"=sd(rCCin),"rCCout"=sd(rCCout))]
    
  dm2 <- rbind(data.table("motif"="TC",
                    "region"="in",
                    "mean"=dmMeans$rTCin,
                    "sd"=dmSTD$rTCin,
                    "se"=dmSE$rTCin),
               data.table("motif"="TC",
                          "region"="out",
                          "mean"=dmMeans$rTCout,
                          "sd"=dmSTD$rTCout,
                          "se"=dmSE$rTCout),
               data.table("motif"="CC",
                          "region"="in",
                          "mean"=dmMeans$rCCin,
                          "sd"=dmSTD$rCCin,
                          "se"=dmSE$rCCin),
               data.table("motif"="CC",
                          "region"="out",
                          "mean"=dmMeans$rCCout,
                          "sd"=dmSTD$rCCout,
                          "se"=dmSE$rCCout))
               
  ggplot(dm2,aes(x=motif,y=mean,fill=region,color=region)) + geom_bar(stat="identity",position="dodge") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.4, colour=rgb(38,120,178,maxColorValue = 255), alpha=0.9, size=0.5,position = position_dodge(0.9)) +
    theme(panel.background = element_blank(),
          axis.title = element_blank(),
          axis.line = element_line(color="black"),
          axis.text.x = element_blank(),
          legend.position = "none") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_fill_manual(values = c(rgb(233,148,151,maxColorValue = 255), 
                                 rgb(146,187,216,maxColorValue = 255))) +
    scale_color_manual(values = c(rgb(212,42,47,maxColorValue = 255), 
                                  rgb(38,120,178,maxColorValue = 255)))
  
  ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig6/fig6c_",cancer,".tiff"),dpi=300,units="mm",width=70,height=50)
  
  
}
