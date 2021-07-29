library(data.table)
library(ggplot2)

TCNStack <- 109115
TCNStrand <- 1039608
TCNBDNA <- 311504008
GQinSize <- 10104513
GQoutSize <- 2862661023

vars <- c("TCNdensityStack","TCNdensityStrand","TCNdensityBDNA")
vals <- c(TCNStack/GQinSize,TCNStrand/GQinSize,TCNBDNA/(2*GQoutSize))
dt <- data.table("var"=factor(vars, levels=vars),"val"=vals)

ggplot(dt,aes(x=var,y=val,fill=var)) + geom_bar(stat="identity") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
#  scale_fill_manual(values=c(rgb(146,187,216,maxColorValue = 255), 
#                             rgb(233,148,151,maxColorValue = 255))) +
  scale_fill_manual(values=c(rgb(203,73,123,maxColorValue = 255), 
                             rgb(161,201,171,maxColorValue = 255),
                             rgb(139,197,229,maxColorValue = 255))) +
  theme(panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig4/fig4a.tiff"),units="mm",width=50,height=60,dpi=300)

