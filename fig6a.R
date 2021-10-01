library(data.table)

data <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek3/united_X0.txt", sep='\t')
data <- data.table(data)
setnames(data,c("cancer","structure","isAPOBEC","sample","trgIn","cntIn","trgOut","cntOut","trgOutOld","cntOutOld","sign"))

dt <- unique(data[,.(structure,isAPOBEC,trgIn,trgOut,trgOutOld)])

dt1 <- dt[isAPOBEC==1]
dt2 <- dt[isAPOBEC==0]

d <- merge(dt1,dt2,by="structure")
d[, ratio := trgIn.x / trgIn.y]

bdnaRatio <- unique(dt[isAPOBEC==1,trgOut]) / unique(dt[isAPOBEC==0,trgOut])

final <- rbind(d[,.(structure,ratio)],data.table("structure"="B-DNA","ratio"=bdnaRatio))

final[structure=="apr", FullName := "A-phased repeats"]
final[structure=="dr", FullName := "Direct repeats"]
final[structure=="gq", FullName := "G-quadruplex"]
final[structure=="ir", FullName := "Inverted repeats"]
final[structure=="mr", FullName := "Mirror repeats"]
final[structure=="str", FullName := "Short tandem repeats"]
final[structure=="z", FullName := "Z-DNA"]
final[structure=="B-DNA", FullName := "B-DNA"]

final$FullName <- factor(final$FullName, levels=c("G-quadruplex",
                                                  "A-phased repeats",
                                                  "Inverted repeats",
                                                  "Mirror repeats",
                                                  "Direct repeats",
                                                  "Short tandem repeats",
                                                  "Z-DNA",
                                                  "B-DNA"))

ggplot(final,aes(x=FullName,y=ratio)) + geom_bar(stat="identity", width=0.7,
                                              fill=rgb(146,187,216,maxColorValue = 255),
                                              color=rgb(38,120,178,maxColorValue = 255)) +
 # geom_errorbar(aes(ymin=m-sd, ymax=m+sd),width=0.2, colour=rgb(38,120,178,maxColorValue = 255), alpha=0.9, size=0.5) +
  theme(panel.background = element_blank(),
        axis.text.x = element_blank(),
       # axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(color="black"),
        legend.position = "bottom") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
#   scale_fill_manual(values=c(rgb(203,73,123,maxColorValue = 255)))#, 
#                               rgb(161,201,171,maxColorValue = 255), 
#                               rgb(242,234,102,maxColorValue = 255),
#                              rgb(139,197,229,maxColorValue = 255)))    

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig6/fig6a.tiff"),dpi=300,units="mm",width=70,height=40)







