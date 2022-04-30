library(tidyverse)
library(readxl)
library(ggsignif)
library(resample)

setwd("_")
Counts<-read_csv("2021_12_CellCounts_FromPhospho.csv") %>% select(-...1)

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


cellCtPermTest<-function(x) {
  pval<-permtest((x %>% dplyr::filter(Group=="Control"))$percentTotal,(x %>% 
                                                                        dplyr::filter(Group=="Patient"))$percentTotal,
                 permutations=20000)
  return(tibble(pval=pval[[1]],Cell=unique(x$Cell)))
}

cellCtPVals<-Counts %>% split(Counts$Cell) %>% map_dfr(cellCtPermTest)

cellCtPVals<-cellCtPVals %>% mutate(Signif=case_when(
  pval < .05 ~ "Significant",
  pval > .05 ~ "Insignif"
))

Counts<-full_join(cellCtPVals,Counts)
# 
# Counts<-Counts %>% mutate(u116=case_when(
#   grepl("UCLA116",Subject)~"UCLA116",
#   grepl("ControlUCLA116",Subject)~"ControlUCLA116",
#   T~"Not"))
# 
# Counts<-Counts %>% mutate(u116c=case_when(
#     grepl("Control",Subject)&grepl("UCLA116",u116)~"ControlUCLA116",
#     !grepl("Control",Subject)&grepl("UCLA116",u116)~"UCLA116",
#     T~"Not"))

ggplot(Counts %>% dplyr::filter(Signif=="Significant"),aes(x=Group,y=percentTotal,color=Group))+
  scale_color_manual(values = c("grey", "#F7A165"))+
  facet_wrap(~Cell,scales="free_y")+
  geom_jitter()+
  stat_summary(alpha=0.5,outlier.shape = NA,outlier.size=0,fun.data=cibox, position=position_dodge(0.95),geom='boxplot') +
  geom_signif(color="black",comparisons = list(c("Control","Patient")),
              test=permtest, test.args = (permutations=2000),
              step_increase = 0.12, textsize=2,
              map_signif_level=function(p)sprintf("p = %.2g", p))+
  theme_bw(base_size=15)+
  theme(legend.position="none",text=element_text(size=3))+
  NULL

ggsave("FigureS1A_CellCtGlobal.png",width=175,height=250,units="mm",scale=1,dpi=3000)


#IDDA<-read_xlsx("/Users/humzakhan/Desktop/2021_Lab/CVID\ Project\ Organization/NEW_IDDAScore_UA_2022_03.xlsx",n_max=200) %>% 
#  dplyr::filter(!is.na(B)) %>% mutate(A=as.numeric(A)) %>% mutate(M=100*M,O=100*O)

#IDDA<-IDDA %>% mutate(IDDAScore=(((A+B+C+D+E+F+G+H+I+J+K+L)/(M/150))+ifelse(N<40,N*.1,4)+ifelse(O<10,O*.8,8)+
 #                                  P+Q+R+S)) %>% rename(Subject='UCLA ID2')

#ggplot(IDDA,aes(y=IDDAScore,fill=`UCLA ID2`))+geom_histogram()+theme_minimal()
  
