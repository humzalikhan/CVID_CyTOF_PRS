library(tidyverse)
library(readxl)
library(sunburstR)
library(ComplexHeatmap)
library(DataExplorer)
library(ggpubr)
library(MetBrewer)

setwd("/Users/humzakhan/Desktop/CyTOF CSVs")

# idda -----
IDDA<-read_xlsx("/Users/humzakhan/Desktop/2021_Lab/CVID\ Project\ Organization/NEW_IDDAScore_UA_2022_03.xlsx",n_max=200) %>% 
  dplyr::filter(!is.na(B)) %>% mutate(A=as.numeric(A)) %>% mutate(M=100*M,O=100*O)

IDDA<-IDDA %>% mutate(IDDAScore=(((A+B+C+D+E+F+G+H+I+J+K+L)/(M/150))+ifelse(N<40,N*.1,4)+ifelse(O<10,O*.8,8)+
                                   P+Q+R+S)) %>% rename(Subject='UCLA ID2')
# counts----

Counts<-read_csv("2021_12_CellCounts_FromPhospho.csv") %>% select(-...1)

pivotCounts<-Counts %>% select(Subject,Cell,percentTotal,Group) %>% pivot_wider(names_from = Cell, values_from = percentTotal)
pivotCounts[is.na(pivotCounts$CD21loIgD_Memory_B_Cells),]$CD21loIgD_Memory_B_Cells<-0
pivotCounts[is.na(pivotCounts$Plasmablast),]$Plasmablast<-0

cCounts<-pivotCounts %>% dplyr::filter(Group=="Control")
cCountsIDDA<-cCounts
cCountsIDDA$IDDA<-NA

pCounts<-pivotCounts %>% dplyr::filter(Group=="Patient")
pCountsIDDA<-IDDA %>% select(Subject,IDDAScore) %>% full_join(pCounts)
pCountsIDDA<-pCountsIDDA %>% mutate(IDDACat=case_when(
  IDDAScore>12 ~ "IDDAhi",
  IDDAScore<12 ~ "IDDAlo"
))

plot_correlation(cCounts[,3:29],title="Control Cell Count Correlation Matrix")
plot_correlation(pCounts[,3:29],title="Patient Cell Count Correlation Matrix")

plot_correlation(pCountsIDDA[,c(2,4:30)],
                 title="Patient Cell Count-IDDA Score Correlation")

cor(pCountsIDDA[,c(2,4:30)]) %>% as_tibble()
# monocyte pos correlated with IDDA Score. b cells neg corr

iddaPCountCorr<-cor(pCountsIDDA[,c(2,4:30)]) 

cols<-rev(c("#FF2A00","#009BFF"))

Heatmap(as.matrix(iddaPCountCorr[2:nrow(iddaPCountCorr),1]),cluster_rows = T,
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8),column_title="IDDA Score-Cell Count Correlation")

pCountsIDDA2<-pCountsIDDA %>% dplyr::filter(Subject!="UCLA102")

iddaPCountCorr2<-cor(pCountsIDDA2[,c(2,4:30)]) 

Heatmap(as.matrix(iddaPCountCorr2[2:nrow(iddaPCountCorr2),1]),cluster_rows = T,
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8),column_title="IDDA Score-Cell Count Correlation without UCLA102")

pCountsIDDA2<-pCountsIDDA %>% dplyr::filter(Subject!="UCLA102")

iddaCorr<-iddaPCountCorr[,1] %>% as.data.frame() %>% 
  tibble::rownames_to_column("Feature") %>% 
  as_tibble() %>% 
  rename(r='.')

iddaCorr<-iddaCorr[order(iddaCorr$r,decreasing = T),] %>% view()

iddaPCountCorr2<-iddaPCountCorr2[,1] %>% as.data.frame() %>% 
  tibble::rownames_to_column("Feature") %>% 
  as_tibble() %>% 
  rename(r='.')

iddaPCountCorr2<-iddaPCountCorr2[order(iddaPCountCorr2$r,decreasing = T),] %>% view()


# manually do correlation matrix
# cor(cCounts[,3:29]) %>% as.data.frame() %>% view()
# cor(pCounts[,3:29]) %>% as_tibble()
# 
# cor_args_list<-pCounts[,3:29]
# cor_args<-NULL
# plot_data <- reshape2::melt(cor(pCounts[,3:29]))
# 
# ggplot(plot_data,aes(x = Var1, y = Var2, fill = value))+
#   geom_tile()

ggplot(pivotCounts,aes(x=`Effector Memory CD4 T Cell`, CD21hiSwitched_Memory_B_Cells,color=Group)) +
  geom_point()+
  facet_grid(.~Group)+
  theme_bw()


IDDACounts<-Counts %>% full_join(IDDA %>% rename(Subject='UCLA ID2')) %>% select(IDDAScore, Cell, Subject, n, TotalCells, percentTotal, Group)

IDDACountspt<-IDDACounts %>% dplyr::filter(Group=="Patient") 

IDDACountspt<-IDDACountspt %>% split(IDDACountspt$Cell)

correlate<-function(x) {
  x<-x
  y<-cor(x$percentTotal,x$IDDAScore,use="complete.obs")
  return(tibble(Cell=unique(x$Cell),pearsonCoeff=y))
}



corr<-map_dfr(IDDACountspt,correlate)

corr2<-map_dfr(IDDACountspt,correlate2)

corr[order(corr2$pearsonCoeff,decreasing = T),]

corr2<-corr2 %>% as.data.frame()

row.names(corr2)<-corr2$Cell

Heatmap((corr2[,2]),cluster_rows = T,
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8),)

Heatmap(as.matrix(corr2[2:nrow(corr2),1]),cluster_rows = T,
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8),column_title="IDDA Score-Cell Count Correlation")



