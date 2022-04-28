library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(ggpubr)
library(MetBrewer)

setwd("/Users/humzakhan/Desktop/CyTOF CSVs")

# idda -----
IDDA<-read_xlsx("/Users/humzakhan/Desktop/2021_Lab/CVID\ Project\ Organization/NEW_IDDAScore_UA_2022_03.xlsx",n_max=200) %>% 
  dplyr::filter(!is.na(B)) %>% mutate(A=as.numeric(A)) %>% mutate(M=100*M,O=100*O)

IDDA<-IDDA %>% mutate(IDDAScore=(((A+B+C+D+E+F+G+H+I+J+K+L)/(M/150))+ifelse(N<40,N*.1,4)+ifelse(O<10,O*.8,8)+
                                   P+Q+R+S)) %>% rename(Subject='UCLA ID2')
# 
# IDDA2<-read_xlsx("/Users/humzakhan/Desktop/2021_Lab/CVID\ Project\ Organization/NEW_IDDAScore_UA_2022_03.xlsx",n_max=200) %>% 
#   dplyr::filter(!is.na(B)) %>% mutate(A=as.numeric(A)) %>% mutate(M=100*M,O=100*O)
# 
# # check usf2 N
# IDDA2<-IDDA2 %>% mutate(IDDAScore=(((A+B+C+D+E+F+G+H+I+J+K+L)/(M/150))+ifelse(N<40,N*.1,4)+ifelse(O<10,O*.8,8)+
#                                    P+Q+R+S))
# 
# ggplot(IDDA,aes(x=IDDAScore,fill=`UCLA ID`))+geom_histogram()
# ggplot(IDDA2,aes(x=IDDAScore,fill=`UCLA ID`))+geom_histogram()
# 
# mean(IDDA$IDDAScore)
# IDDA %>% dplyr::filter(IDDAScore>11.333)
# 
# IDDA2 %>% dplyr::filter(IDDAScore>13)

# counts----

Counts<-read_csv("2021_12_CellCounts_FromPhospho.csv") %>% select(-...1)

# cell count IDDA strat --------
countIDDA<-IDDA %>% full_join(Counts) %>% select(percentTotal,Group,Cell,Subject,IDDAScore) %>% 
  mutate(IDDACat=case_when(
    IDDAScore>12 ~ "IDDAhi",
    IDDAScore<12 ~ "IDDAlo",
    is.na(IDDAScore) ~ "Control"
  ))

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
  b = resample::bootstrap(x,R=20000,statistic=mean)$observed;
  data.frame(y = b,
             ymin = b,
             ymax = b)
}

compareCountsIDDA<-function(x) {
  iddahi<-x %>% dplyr::filter(IDDACat=="IDDAhi")
  iddalo<-x %>% dplyr::filter(IDDACat=="IDDAlo")
  control<-x %>% dplyr::filter(IDDACat=="Control")
  control.IDDAhi.pval<-permtest(unique(control$percentTotal),unique(iddahi$percentTotal))[[1]]
  #print(pval1)
  IDDAhi.IDDAlo.pval<-permtest(unique(iddahi$percentTotal),unique(iddalo$percentTotal))[[1]]
  control.IDDAlo.pval<-permtest(unique(control$percentTotal),unique(iddalo$percentTotal))[[1]]
  #print(pval2)
  return(tibble(Cell=x$Cell,control.IDDAlo.pval,control.IDDAhi.pval,IDDAhi.IDDAlo.pval)[1,])
}

countIDDASplit<-countIDDA %>% split(list(countIDDA$Cell))

iddaStratCellCtPVals<-map_dfr(countIDDASplit,~compareCountsIDDA(.))

iddaStratCellCtPVals<-bind_rows(iddaStratCellCtPVals)

x<-iddaStratCellCtPVals %>% dplyr::filter(IDDAhi.IDDAlo.pval<.05)

Cells<-x$Cell %>% paste(collapse = "|")

cols<-c("#C5D5E5","#FA461D","#1D88FA")

setwd("/Users/humzakhan/Desktop/2021_Lab/CVID\ Project\ Organization/Paper/Plots")
rel<-countIDDA %>% dplyr::filter(grepl(Cells,Cell))
rel$Cell<-factor(rel$Cell,levels = c("CD21hiSwitched_Memory_B_Cells", "CD21loSwitched_B_Cells", "Regulatory T Cell", "CD21loNaive_B_Cells"))
#rel$Cell<-fct_relevel(rel$Cell, "CD21loNaive_B_Cells", after = Inf)

plot<-ggplot(rel,aes(x=IDDACat,y=Cell,color=IDDACat,size=percentTotal))+
  geom_point()+
  theme_bw()+
  #facet_grid(.~Cell)+
  scale_color_manual(values=cols)+
  # geom_signif(color="black",comparisons = list(c("Control","Patient")),
  #             test=permtest, test.args = (permutations=2000),
  #             step_increase = 0.12, textsize=2.5,
  #             map_signif_level=function(p)sprintf("p = %.2g", p))+
  theme(legend.position="none",axis.text.y=element_blank())
  NULL

ggsave(plot=plot,"Figure2b_CellCts.png",width=60,height=60,units="mm",scale=1,dpi=3000)

# save this 

setwd("/Users/humzakhan/Desktop/2021_Lab/CVID\ Project\ Organization/Paper/Plots")

cols<-c("#C5D5E5","#FA461D","#1D88FA")

plot<-ggplot(countIDDA %>% dplyr::filter(grepl(Cells,Cell)),aes(x=IDDACat,y=percentTotal,color=IDDACat))+
  geom_jitter(height=0)+
  theme_bw()+
  facet_wrap(.~Cell,scales="free_y")+
  stat_summary(alpha=0.5,outlier.shape = NA,outlier.size=0,fun.data=cibox, position=position_dodge(0.95),geom='boxplot') +
  geom_signif(color="black",comparisons = list(c("Control","IDDAhi"),c("Control","IDDAlo"),c("IDDAhi","IDDAlo")),
              test=permtest, test.args = (permutations=2000),
              step_increase = 0.12, textsize=2.5,
              map_signif_level=function(p)sprintf("p = %.2g", p))+
  scale_color_manual(values=cols)+
  theme(legend.position="none")+
  NULL

ggsave(plot=plot,"Figure2c_CellCts_Signif.png",width=120,height=140,units="mm",scale=1,dpi=3000)

