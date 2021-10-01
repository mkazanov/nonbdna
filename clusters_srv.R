library(data.table)
library(openxlsx)

args = commandArgs(trailingOnly=TRUE)

chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
          "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")

wb <- loadWorkbook("/data1/DATA/NONBDNA_CLUSTERS/TabS5_journal.pbio.3000464.s011.xlsx")
dataClusters <- read.xlsx(wb,sheet=2)
dataClusters <- data.table(dataClusters)
dataClusters[cancer_type == "Bladder-TCC ", cancer := "BLCA"]
dataClusters[cancer_type == "Breast-AdenoCA ", cancer := "BRCA"]
dataClusters[cancer_type == "Head-SCC ", cancer := "HNSC"]
dataClusters[cancer_type == "Cervix-SCC ", cancer := "CESC"]
dataClusters[cancer_type == "Lung-AdenoCA ", cancer := "LUAD"]
dataClusters[cancer_type == "Lung-SCC ", cancer := "LUSC"]

if (length(args)!=0) {
  arg1 <- args[[1]]
  dataClusters <- dataClusters[CG_cluster_classification_more_3mutations == arg1]
} else {
  arg1 <- "noargs"
}

activity <- read.csv("/data1/DATA/NONBDNA_CLUSTERS/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))
activity[,cancer := substr(project,1,4)]

strs <- list()
for(str in c("apr","dr","gq","ir","mr","str","z"))
{
 dt <- data.table()
 for(chr in chrs){
  dttmp <- read.csv(paste0("/data1/DATA/NONBDNA_CLUSTERS/human_hg19.tsv/",str,"/tsv/",chr,"_",toupper(str),".tsv"), sep='\t')
  dt <- rbind(dt, data.table(dttmp))
 }
 strs[str] <- list(dt)
}

results <- data.table()
for(canc in c("BLCA","BRCA","HNSC","CESC","LUAD","LUSC")){
  print(canc)
  dc <- dataClusters[cancer == canc]
  da <- activity[cancer == canc]
  samples <- unique(da$sample)
  for(s in samples){
    print(s)
    dcs <- dc[Tumor_Sample_Barcode == s]
    if(nrow(dcs) == 0){
      print(paste0(canc,': ',s))
    }
    
    cntList <- list()
    sizeList <- list()
    totalsize <- dcs[,sum(Cluster_Length)]
    totalsize_check <- 0
    for(str in c("apr","dr","gq","ir","mr","str","z")){
      dt <- data.table(strs[[str]])
      cnt <- 0
      sizeIn <- 0
      sizeOut <- 0
      for(i in 1:nrow(dcs)){
        chr <- dcs[i]$Chromosome
        spos <- dcs[i]$Cluster_start
        epos <- dcs[i]$Cluster_end
        dtt <- copy(dt) 
        dtt[Sequence_name == chr, maxStart := ifelse(spos > Start, spos, Start)]
        dtt[Sequence_name == chr, minEnd := ifelse(epos > Stop, Stop, epos)]
        #dtintersect1 <- dt[Sequence_name == chr & ((spos >= Start & Start >= epos) | (spos >= Stop & Stop >= epos))]
        dtintersect2 <- dtt[!is.na(maxStart) & ((minEnd - maxStart + 1) > 0)]
        
        if(nrow(dtintersect2) > 0){
          curSize <- dtintersect2[,sum(minEnd - maxStart + 1)]
        } else {
          curSize <- 0
        }
        
        sizeIn <- sizeIn + curSize
      }
      sizeList[str] <- sizeIn
    }
    
    results <- rbind(results, data.table("cancer"=canc,
                                         "sample"=s,
                                         "apr"=cntList[['apr']],
                                         "aprSize"=sizeList[['apr']],
                                         "dr"=cntList[['dr']],
                                         "drSize"=sizeList[['dr']],
                                         "gq"=cntList[['gq']],
                                         "gqSize"=sizeList[['gq']],
                                         "ir"=cntList[['ir']],
                                         "irSize"=sizeList[['ir']],
                                         "mr"=cntList[['mr']],
                                         "mrSize"=sizeList[['mr']],
                                         "str"=cntList[['str']],
                                         "strSize"=sizeList[['str']],
                                         "z"=cntList[['z']],
                                         "zSize"=sizeList[['z']],
                                         "cnt"=nrow(dcs),
                                         "totalsize"=totalsize))
  }
}

write.csv(results,paste0("nonbdna_clusters_",arg1,".txt"))
