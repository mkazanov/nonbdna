library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(ggplotify)
library(cowplot)

data <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek3/united_X0.txt",sep='\t',header = TRUE)
data <- data.table(data)
setnames(data,c("cancer","structure","isAPOBEC","sample","trgIn","cntIn","trgOut","cntOut","trgOutOld","cntOutOld","sign"))
data[sign < 5 | sign > 95, signBinary:=1]
data[is.na(signBinary), signBinary := 0]
data$signBinary <- as.factor(data$signBinary)
data[, ratio := log((cntIn/trgIn)/(cntOut/trgOut))]

activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))

data <- merge(data,activity,by="sample",all.x = TRUE)

structures <- unique(data$structure)
cancers <- unique(data$cancer)

plots <- list()

i <- 1
for(s in structures)
{
  for(c in cancers)
  {
      dt <- rbind(data[cancer==c & isAPOBEC == 1 & structure == s & enrichment > 2.0],
                  data[cancer==c & isAPOBEC == 0 & structure == s & enrichment <= 2.0])
      dt$isAPOBEC <- as.factor(dt$isAPOBEC)
      dt$signBinary <- as.numeric(dt$signBinary)
      dt$signBinary <- - dt$signBinary
      dt$signBinary <- as.factor(dt$signBinary)
      
      p <- ggscatter(dt,x="enrichment",y="ratio", fill="isAPOBEC",color="isAPOBEC",shape="signBinary",add="reg.line",conf.int=TRUE,
                     xticks.by=0.5) +
#        ggtitle(paste0("Structure=",s," , Cancer=",c)) +
#        geom_violin(alpha=0.4, position = position_dodge(width = .75),color="black") +
#        geom_boxplot(notch = TRUE,  outlier.size = -1, color="black", alpha = 0.7)+
        ggpubr::fill_palette("jco") +
        theme_light() +
        theme(plot.title = element_text(size=8),
              axis.title = element_blank(),
              legend.position = "none")

      ydens <- axis_canvas(p,axis="y",coord_flip = TRUE) +
        geom_density(data = dt, aes(x = ratio, fill = isAPOBEC),
                     alpha = 0.7, size = 0.2)+
        coord_flip()+
        ggpubr::fill_palette("jco")
      
      yplot <- ggboxplot(dt, x = "isAPOBEC", y = "ratio",
                         color = "isAPOBEC", fill = "isAPOBEC", palette = "jco",
                         alpha = 0.5)
      
      p2 <- insert_yaxis_grob(p, yplot, grid::unit(.2, "null"), position = "right")

      ggsave(plot=p2,filename=paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/pics/fig3/",s,"_",c,".tiff"),units="mm",width=140,height=100,dpi=300)
    
  }
}


