library(openxlsx)
library(reshape2)

wb <- loadWorkbook("/Users/mar/BIO/PROJECTS/STUDENTS/DENISOVA/FinalResultsBDNA/ss_bg_08.xlsx")
wb2 <- loadWorkbook("/Users/mar/BIO/PROJECTS/STUDENTS/DENISOVA/FinalResultsBDNA/structure_BDNA.xlsx")

data <- read.xlsx(wb, sheet=1)
data <- data.table(data)
data2 <- read.xlsx(wb2, sheet=1)
data2 <- data.table(data2)

dt <- data[X1 == "All"]
dt2 <- data2[X1 == "all"]

dt[,X1 := NULL]
dt2[,X1 := NULL]

dtm <- melt(dt)
dt2m <- melt(dt2)

final <- merge(dtm,dt2m,by="variable")
final <- data.table(final)
final[,ratio := value.x/value.y]

ggplot(final,aes(x=variable,y=ratio)) + geom_bar(stat="identity",fill="lightblue") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle=90),
        axis.title.x = element_blank(),
        axis.line = element_line(color="black"),
        legend.position = "bottom") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

  



