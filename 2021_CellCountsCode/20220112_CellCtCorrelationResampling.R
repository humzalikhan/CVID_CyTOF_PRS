library(tidyverse)
library(readxl)
library(ComplexHeatmap)

IDDA<-read_csv("NEW_IDDAScore_HK_Calculated_2022_03.csv")

# counts----

Counts<-read_csv("/Users/humzakhan/Desktop/CyTOF\ CSVs/2021_12_CellCounts_FromPhospho.csv") %>% select(-...1)

IDDACounts<-Counts %>% full_join(IDDA) %>% select(IDDAScore, Cell, Subject, n, TotalCells, percentTotal, Group)

IDDACountspt<-IDDACounts %>% dplyr::filter(Group=="Patient") 
# make 0s for things 
IDDACountspt<-IDDACountspt %>% split(IDDACountspt$Cell)

correlateResample <-function(x,R=12) {
  x2<-x %>% slice_sample(n=R,replace=T)
  notIn<-setdiff(x,x2)
  y<-cor(x2$percentTotal,x2$IDDAScore,use="complete.obs")
  return(tibble(Cell=unique(x$Cell),pearsonCoeff=y,excluded=paste(notIn$Subject,collapse=",")))
}

corResampleDf<-as_tibble(NULL)
for (i in c(1:1000)) {
  cellCor<-map_dfr(IDDACountspt,~correlateResample(.,R=11))
  cellCor$Iter<-i
  corResampleDf<-rbind(cellCor,corResampleDf)
}

corResampleDfMeans<-corResampleDf %>% split(corResampleDf$Cell) %>% map(~(mean(.$pearsonCoeff)))

#corResampleDfMeansPerm<-corResampleDf %>% split(corResampleDf$Cell) %>% map(~(permmean(.$pearsonCoeff)))

corResampleDfMeans<-corResampleDfMeans %>% as_tibble() %>% 
  t() %>% 
  as.data.frame() 

corResampleDfMeans %>% as.data.frame() %>% arrange(desc(V1))

png("Figure2d_CellCtCorr.png",width=120,height=100,units="mm",res=3000)
Heatmap(as.matrix(corResampleDfMeans),cluster_rows = T,
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8))
setwd("/Users/humzakhan/Desktop/2021_Lab/CVID\ Project\ Organization/Paper/Plots")
# 
# # 
# # corResampleDfMeansPerm<-corResampleDfMeansPerm %>% as_tibble() %>% 
# #   t() %>% 
# #   as.data.frame() %>% 
# #   rownames_to_column() %>% 
# #   as_tibble() %>% 
# #   filter(!grepl("ymin",rowname),!grepl("ymax",rowname))
# 
# corResampleDfMeans %>% filter(abs(V1)>.35)
# 
# ggplot(corResampleDf %>% filter(grepl("Monocyte",Cell)),aes(x=pearsonCoeff))+
#   geom_histogram()+
#   theme_bw()+
#   facet_grid(.~Cell)
# 
# ggplot(corResampleDf %>% filter(grepl("B_Cell",Cell)),aes(x=pearsonCoeff))+
#   geom_histogram()+
#   theme_bw()+
#   facet_grid(.~Cell)
# 
# ggplot(corResampleDf %>% filter(grepl("CD4 T Cell",Cell)),aes(x=pearsonCoeff))+
#   geom_histogram()+
#   theme_bw()+
#   facet_grid(.~Cell)
# 
# ggplot(corResampleDf %>% filter(grepl("CD8 T Cell",Cell)),aes(x=pearsonCoeff))+
#   geom_histogram()+
#   theme_bw()+
#   facet_grid(.~Cell)
# 
# 
# ggplot(corResampleDf %>% filter(grepl("Dendritic Cell",Cell)),aes(x=pearsonCoeff))+
#   geom_histogram()+
#   theme_bw()+
#   facet_grid(.~Cell)
# 
# ggplot(corResampleDf %>% filter(grepl("Dendritic Cell",Cell)),aes(x=pearsonCoeff))+
#   geom_histogram()+
#   theme_bw()+
#   facet_grid(.~Cell)


