library(data.table)
library(ggplot2)
library(reshape2)

clusterType <- list()
clusterType[1] <- "CorG_coord"
clusterType[2] <- "coord_with_terminal_CG"
clusterType[3] <- "CG_single_switch"
clusterType[4] <- "CG_multiple_switch"

i <- 4

data <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Gordenin/nonbdna_clusters_",clusterType[[i]],".txt"))
data <- data.table(data)
data[,X:=NULL]

activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))
activity[,cancer := substr(project,1,4)]

data <- merge(data,activity,by=c("cancer","sample"))
data[,apr_fraction := apr/cnt]
data[,dr_fraction := dr/cnt]
data[,gq_fraction := gq/cnt]
data[,ir_fraction := ir/cnt]
data[,mr_fraction := mr/cnt]
data[,str_fraction := str/cnt]
data[,z_fraction := z/cnt]


for(canc in c("BLCA","BRCA","HNSC","CESC","LUAD","LUSC")){
  dt <- data[cancer == canc]
  dt <- dt[,.(enrichment,apr_fraction,dr_fraction,gq_fraction,ir_fraction,mr_fraction,str_fraction,z_fraction)]
  dtmelt <- melt(dt,id.vars = "enrichment")
  ggplot(dtmelt, aes(x=enrichment,y=value,fill=variable,color=variable)) + geom_point() 
  ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Gordenin/",clusterType[[i]],"/point_",canc,".tiff"))
}
