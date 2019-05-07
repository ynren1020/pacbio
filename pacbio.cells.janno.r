############################################2019-05-07#############################
##################pacbio cell line junction data find igf1ras1 junction format#####
##################summarize the numbers of different format########################
#grep chr15 /home/yang4414/tywang/projects/pacbio_IGF1R-AS1/data/mapping/*.janno###
##| grep 988 | grep ? |grep ct-AC >pacbio.cells.igf1ras1.janno#####################
###################################################################################
library(dplyr)
library(tidyr)
library(stringr)


input_file<-"pacbio.cells.igf1ras1.janno"
dat<-read.delim(input_file,header = FALSE,stringsAsFactors = FALSE)
for (i in 1:nrow(dat)){
  dat$cells[i]<-str_split_fixed(dat$V1[i],"/",n=9)[9]
}

dat$cells<-str_replace_all(dat$cells,".janno:chr15","")
dat<-dat%>%separate(cells,c("cells","type"),'[.]')

##filter junctions with score less than 5##
dat<-filter(dat,V5>5)

dat<-dat%>%
  group_by(type)%>%
  arrange(V2,V3)

dat<-dat%>%unite("coord",c("V2","V3"),sep=":")

table(dat$coord,dat$type)

#                    A549 DMS-53 H727 LNCaP VCaP
#98831560:98839791    0      0    1     0    2
#98831563:98839635    0      0    0     0    1
#98831563:98839791    2      2    2     0    2
#98831798:98839791    0      0    1     0    0
#98831922:98839791    0      0    1     0    0
#98831953:98839772    0      0    1     0    0
#98831959:98835340    0      0    1     0    1
#98831959:98838253    1      0    2     0    0
#98831959:98839434    0      0    1     0    1
#98831959:98839635    0      0    1     0    0
#98831959:98839754    0      0    2     0    0
#98831959:98839772    0      0    1     0    0
#98831959:98839791    2      2    2     2    2
#98831959:98839797    0      0    2     0    0
#98831959:98843309    0      0    1     0    0
#98831959:98843675    0      0    1     0    0
#98831964:98835340    1      1    2     0    1
#98831964:98838253    0      0    1     0    0
#98831964:98839754    0      0    1     0    0
#98831964:98839791    2      2    2     1    2
#98840709:98843671    0      0    1     0    0



