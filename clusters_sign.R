library(data.table)
library(ggplot2)
library(reshape2)

clusterType <- list()
clusterType[1] <- "CorG_coord"
clusterType[2] <- "coord_with_terminal_CG"
clusterType[3] <- "CG_single_switch"
clusterType[4] <- "CG_multiple_switch"
clusterType[5] <- "noargs"

i <- 4

data <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Gordenin/nonbdna_clusters_",clusterType[[i]],".txt"))
data <- data.table(data)
data[,X:=NULL]

activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))
activity[,cancer := substr(project,1,4)]

data <- merge(data,activity,by=c("cancer","sample"))

dataBDNA <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Denek3/united_X0.txt", sep='\t')
dataBDNA <- data.table(dataBDNA)
setnames(dataBDNA,c("Cancer","structure","isAPOBEC","sample","trgIn","cntIn","trgOut","cntOut","trgOutOld","cntOutOld","sign"))
sizesIn <- unique(dataBDNA[isAPOBEC == 0, .(structure,trgIn,trgOut)])
sizeOut <- unique(sizesIn$trgOut)

data[,apr_fraction := aprSize/sizesIn[structure == "apr",trgIn]]
data[,dr_fraction := drSize/sizesIn[structure == "dr",trgIn]]
data[,gq_fraction := gqSize/sizesIn[structure == "gq",trgIn]]
data[,ir_fraction := irSize/sizesIn[structure == "ir",trgIn]]
data[,mr_fraction := mrSize/sizesIn[structure == "mr",trgIn]]
data[,str_fraction := strSize/sizesIn[structure == "str",trgIn]]
data[,z_fraction := zSize/sizesIn[structure == "z",trgIn]]
data[,bdna_fraction := (totalsize - (aprSize + drSize + gqSize + irSize + mrSize + strSize + zSize)) / sizeOut]

data[, apr_ratio := log(apr_fraction/bdna_fraction)]
data[, dr_ratio := log(dr_fraction/bdna_fraction)]
data[, gq_ratio := log(gq_fraction/bdna_fraction)]
data[, ir_ratio := log(ir_fraction/bdna_fraction)]
data[, mr_ratio := log(mr_fraction/bdna_fraction)]
data[, str_ratio := log(str_fraction/bdna_fraction)]
data[, z_ratio := log(z_fraction/bdna_fraction)]

dataAHigh <- data[enrichment >= 2.0]
dataALow <- data[enrichment < 2.0] 


results <- data.table()
for(canc in c("BLCA","BRCA","HNSC","CESC","LUAD","LUSC")){
 for(str in c("apr","dr","gq","ir","mr","str","z")){
  x <- dataAHigh[,get(paste0(str,"_fraction"))]
  y <- dataAHigh[,bdna_fraction]
  res <- wilcox.test(x,y)
  results <- rbind(results,data.table("cancer"=canc,"structure"=str,"pvalue"=res$p.value))
 }
}

write.csv(results,"/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/Gordenin/cluster_significance.txt",row.names = F, quote = F)
