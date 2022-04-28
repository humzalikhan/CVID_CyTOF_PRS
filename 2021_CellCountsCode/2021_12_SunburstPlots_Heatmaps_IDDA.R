library(tidyverse)
library(readxl)
library(sunburstR)
library(ComplexHeatmap)
library(DataExplorer)
library(ggpubr)

setwd("/Users/humzakhan/Desktop/CyTOF CSVs")
Counts<-read_csv("2021_12_CellCounts_FromPhospho.csv") %>% select(-...1)
Cells<-unique(Counts$Cell) %>% as_tibble()


IDDA<-read_csv("/Users/humzakhan/Desktop/2021_Lab/CVID\ Project\ Organization/2021_IDDAScore_HK_Calculated.csv") %>% 
  select('UCLA ID2',IDDAScore) %>%
  mutate(IDDACat=case_when(
    IDDAScore>11.3 ~ "IDDAhi",
    IDDAScore<11.3 ~ "IDDAlo"
  )) %>% rename(Subject='UCLA ID2')

Counts<-Counts %>% full_join(IDDA)
Counts[is.na(Counts$IDDACat),]$IDDACat<-"Control"

myeloid<-"Basophil|Monocyte|Dendritic|Eosinophil|Neutrophil"
lymphoid<-"B_Cell|T Cell|Killer|Plasma"

Cells<-Cells %>% mutate(InLists = case_when(
  grepl(myeloid,value) ~ 'Myeloid',
  grepl(lymphoid,value) ~ 'Lymphoid',
)) 


Cells<-Cells %>% mutate(Subtype = case_when(
  grepl("Monocyte",value) ~ 'Monocyte',
  grepl("B_Cell",value) ~ 'B Cell',
  grepl("blast",value) ~ 'B Cell',
  grepl("T Cell",value) ~ 'T Cell',
  grepl("T Cell",value) ~ 'T Cell',
  grepl("phil",value) ~ 'Granulocyte',
  grepl("Dendritic",value) ~ 'Dendritic Cell',
  grepl("Natural",value) ~ 'Natural Killer Cell'
)) %>% rename(Cell=value)

Cells[Cells$Cell=="Natural Killer T Cell",]$Subtype<-"Natural Killer Cell"

Cells<-Cells %>% mutate(CD21 = case_when(
  grepl("CD21lo",Cell) ~ 'CD21lo B Cell',
  grepl("CD21hi",Cell) ~ 'CD21hi B Cell',
))

annotatedCounts<-Cells %>% full_join(Counts)

annotatedCounts2<-annotatedCounts %>% group_by(IDDACat,Cell,InLists,Subtype) %>% 
  summarise(mean=mean(percentTotal)) %>% 
  mutate(path = paste(IDDACat,InLists,Subtype,Cell, sep = "-"))

sunburst(data = data.frame(xtabs(mean~path, annotatedCounts2 %>% dplyr::filter(IDDACat=="Control",
                                                                        Cell!="Neutrophil"))), legend = T,
         withD3=F,
         percent=T,count=T)

sunburst(data = data.frame(xtabs(mean~path, annotatedCounts2 %>% dplyr::filter(IDDACat=="IDDAhi",
                                                                        Cell!="Neutrophil"))), legend = T,
         withD3=F,
         percent=T,count=T)

sunburst(data = data.frame(xtabs(mean~path, annotatedCounts2 %>% dplyr::filter(IDDACat=="IDDAlo",
                                                                        Cell!="Neutrophil"))), legend = T,
         withD3=F,
         percent=T,count=T)

# heatmaps 

CellCountsSpread <- Counts %>% select(-n,-TotalCells) %>%  pivot_wider(names_from = Cell,values_from=percentTotal)

CellCountsSpread<-CellCountsSpread %>% as.data.frame()

rownames(CellCountsSpread) = CellCountsSpread$IDDACat

groupDF = data.frame("IDDACat" = CellCountsSpread$IDDACat)
groupDF = groupDF[order(groupDF$IDDACat),]
#rownames(groupDF) = rownames(FoldChangeCStimClean) 

CellCountsSpread<-CellCountsSpread[order(CellCountsSpread$IDDACat),]
CellCountsSpread[is.na(CellCountsSpread)]<-0
#View(FoldChangeCStimClean1)

Heatmap(as.matrix(scale(CellCountsSpread[,5:31])),cluster_rows =F,
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8),
        column_title = "Cell Counts") + 
  rowAnnotation(Group=groupDF,simple_anno_size=ht_opt$simple_anno_size/3,
                col=list(Group=c("Control" = "grey", "IDDAhi" = "red","IDDAlo"="lightblue2")))+
  rowAnnotation(rn = anno_text(rownames(CellCountsSpread),gp=gpar(fontsize=8)))

# 

CellCountsSpread <- Counts %>% select(-n,-TotalCells) %>% 
  group_by(IDDACat,Cell) %>% 
  summarise(mean=mean(percentTotal,na.rm=T)) %>% 
  pivot_wider(names_from = Cell,values_from=mean)

CellCountsSpread<-CellCountsSpread %>% as.data.frame()

rownames(CellCountsSpread) = CellCountsSpread$IDDACat

groupDF = data.frame("IDDACat" = CellCountsSpread$IDDACat)
groupDF = groupDF[order(groupDF$IDDACat),]
#rownames(groupDF) = rownames(FoldChangeCStimClean) 

CellCountsSpread<-CellCountsSpread[order(CellCountsSpread$IDDACat),]
CellCountsSpread[is.na(CellCountsSpread)]<-0
#View(FoldChangeCStimClean1)

Heatmap(as.matrix(scale(CellCountsSpread[,2:ncol(CellCountsSpread)])),cluster_rows =F,
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8),
        column_title = "Cell Counts") + 
  rowAnnotation(Group=groupDF,simple_anno_size=ht_opt$simple_anno_size/3,
                col=list(Group=c("Control" = "grey", "IDDAhi" = "red","IDDAlo"="lightblue2")))+
  rowAnnotation(rn = anno_text(rownames(CellCountsSpread),gp=gpar(fontsize=8)))

