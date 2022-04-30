library(tidyverse)
library(RColorBrewer)
library(ggsignif)
library(resample)

cibox<-function(x){
  b<-resample::bootstrap(x,R=30000,statistic=mean);
  m = b$observed;
  data.frame(lower=CI.percentile(b)[1],
             upper=CI.percentile(b)[2],
             ymin=CI.percentile(b)[1],
             ymax=CI.percentile(b)[2],
             middle=m)
}

permtest <- function(x,y,permutations=20000,statistic=mean) {
  b <- resample::permutationTest2(data=x, data2 = y,
                                  R=permutations,statistic=statistic)
  list( p.value = b$stats$PValue )
}

permmean <- function(x) {
  b = resample::bootstrap(x,R=5000,statistic=mean)$observed;
  data.frame(y = b,
             ymin = b,
             ymax = b)
}


theme_manish <- function() {
  theme_bw(base_size=10, base_line_size = 0.5, base_rect_size = 0.5) +
    theme(axis.text = element_text(size = rel(0.8), colour = "black", angle = 0)) +
    theme(axis.title = element_text(size = rel(1), colour = "black", angle = 0)) +
    theme(axis.ticks = element_line(colour = "black", size = rel(0.25))) + 
    theme(panel.border  = element_rect(colour = "black", size = rel(0.25))) + 
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(panel.background = element_blank()) +
    theme(plot.background = element_blank()) +
    theme(legend.position="none") +
    theme(strip.background = element_rect(size=rel(0.25)) ) + 
    # theme(legend.key.width=unit(0.3,"cm")) +
    theme(legend.text=element_text(size=6)) +
    theme(legend.title=element_text(size=8)) 
}

setwd("2021_CyTOF_CSVs")

# read all 

#totalPhosph<-bind_rows(all_files)

totalCellCt<-totalPhosph %>% group_by(Subject) %>% summarise(TotalCells=n())

cellCts<-totalPhosph %>% count(Subject,Cell)

cellNums<-full_join(cellCts,totalCellCt)

cellNums<-cellNums %>% mutate(percentTotal=n/TotalCells*100)

cellNums %>% arrange(desc(percentTotal))

cellNums<-cellNums %>% mutate(Group=case_when(
  grepl("Control", Subject) ~ "Control",
  !grepl("Control", Subject) ~ "Patient"
))

# myeloid<-c("Basophil","Monocyte","Eosinophil","Neutrophil","Dendritic")
# 
# lymphocyte<-c("B_Cell","T Cell","Natural Killer","Plasmablast","Dendritic")
# 
# cellNums<-cellNums %>% mutate(CellGroup=case_when(
#   grepl("B_Cell", Subject) ~ "B Cell",
#   grepl("Control", Subject) ~ "Myeloid"
# ))
# 
# unique(cellNums$Cell)

cellCtPermTest<-function(x) {
  pval<-permtest((x %>% filter(Group=="Control"))$percentTotal,(x %>% 
                                                       filter(Group=="Patient"))$percentTotal,
           permutations=20000)
  return(tibble(pval=pval[[1]],Cell=unique(x$Cell)))
}

cellCtPVals<-cellNums %>% split(cellNums$Cell) %>% map_dfr(cellCtPermTest)

write.csv(cellNums,"2021_12_CellCounts_FromPhospho.csv")
write.csv(cellCtPVals,"2021_12_CellCounts_pvals.csv")

