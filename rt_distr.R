library(data.table)

data <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek3/united_RT_X0.txt",sep='\t',header = TRUE)
data <- data.table(data)

setnames(data,c("cancer","structure","isAPOBEC","sample","trgIn","cntIn","trgOut","cntOut","trgOutOld","cntOutOld","sign","RTbin"))

data[,RTbin2 := 6 - RTbin]

dt <- unique(data[isAPOBEC==1,.(structure,RTbin2,trgIn,trgOut)])
dt[,fraction := trgIn/(trgIn+trgOut)]
dt$RTbin2 <- as.factor(dt$RTbin2)

for(s in unique(dt$structure)){
 dt2 <- dt[structure == s]  
 
 dt2group <- dt2[,.(m=mean(fraction),sd=sd(fraction)),by=.(structure,RTbin2)]
 
 ggplot(dt2group,aes(x=RTbin2,y=m)) + geom_bar(stat="identity", width=0.7,
                                                           fill=rgb(146,187,216,maxColorValue = 255),
                                                           color=rgb(38,120,178,maxColorValue = 255)) +
   geom_errorbar(aes(ymin=m-sd, ymax=m+sd),width=0.2, colour=rgb(38,120,178,maxColorValue = 255), alpha=0.9, size=0.5) +
    theme(panel.background = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank(),
          axis.line = element_line(color="black"),
          legend.position = "bottom") +
   scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
 #   scale_fill_manual(values=c(rgb(203,73,123,maxColorValue = 255)))#, 
#                               rgb(161,201,171,maxColorValue = 255), 
#                               rgb(242,234,102,maxColorValue = 255),
 #                              rgb(139,197,229,maxColorValue = 255)))    
 ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig5/fig5a_TC_",s,".tiff"),dpi=300,units="mm",width=70,height=70)
}

dt <- unique(data[isAPOBEC==0,.(structure,RTbin2,trgIn,trgOut)])
dt[,fraction := trgIn/(trgIn+trgOut)]
dt$RTbin2 <- as.factor(dt$RTbin2)

for(s in unique(dt$structure)){
  dt2 <- dt[structure == s]  
  ggplot(dt2,aes(x=RTbin2,y=fraction)) + geom_bar(stat="identity") +
    theme(panel.background = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.line = element_line(color="black"),
          legend.position = "bottom") #+
  #    scale_fill_manual(values=c(rgb(203,73,123,maxColorValue = 255), 
  #                               rgb(161,201,171,maxColorValue = 255), 
  #                               rgb(242,234,102,maxColorValue = 255),
  #                              rgb(139,197,229,maxColorValue = 255)))    
  ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig5/fig5a_all_",s,".tiff"))
}
