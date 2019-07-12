#################2019-06-19#################
##cell1.VCaP.bam relative files#############
##/home/yren/project/pacbio#################
##find transcript's isoforms ###############
############################################
##https://www.rdocumentation.org/packages/GenomicAlignments/versions/1.8.4/topics/cigar-utils
##packages##
library(tidyverse)
library(GenomicAlignments)
##data##
input<-"cell1.VCaP.6cols.chr15.988.txt"
vcap<-read.delim(input,header = FALSE,stringsAsFactors = FALSE)
length(unique(vcap$V1)) #19319
vcap<-vcap[vcap$V4>98800000&vcap$V5>=20,] #mapq>=20, mapq=-10log10(probability of wrongmap) mapq=20, Pwrong=0.01

vcap.reverse<-vcap[vcap$V2==16,]#8071 reverse strand flags=16

length(unique(vcap.reverse$V4)) #365

cigar<-vcap.reverse$V6

##cigar length table##
cigarOpTable(cigar)

#right coordinate## seems not right##
cigarRangesAlongReferenceSpace(cigar,
                               reduce.ranges=TRUE, with.ops=TRUE)[[1]] 


pos2 <- vcap.reverse$V4

cigarRangesAlongReferenceSpace(cigar, pos=pos2, with.ops=TRUE)[[1]]

res1a <- extractAlignmentRangesOnReference(cigar, pos=pos2)
res1b <- cigarRangesAlongReferenceSpace(cigar,
                                        pos=pos2,
                                        ops=setdiff(CIGAR_OPS, "N"),
                                        reduce.ranges=TRUE)
#seqnames<-factor(vcap.reverse$V3)
#extractAlignmentRangesOnReference(cigar, pos=pos2, f=seqnames)
res1b[[1]]

##list to dataframe##
for (i in 1:length(res1b)){
  
}

library(plyr)

names(res1b) <- vcap.reverse$V1
result <- ldply(res1b, data.frame)

result.summ<-result%>%dplyr::group_by(.id)%>%dplyr::summarise(count=n())

result.join<-full_join(result,result.summ,by=".id")

result.join.count2<-result.join[result.join$count==2,]

table(as.factor(result.join$count))

#1     2     3     4     5 
#137 12594  4839    92     5 

##2 exons##
result.join.count2.odd<-result.join.count2 %>% dplyr::filter(row_number() %% 2 == 1) ## Select odd rows
unique(result.join.count2.odd$width)
[1] 1763 1762 1757 1756 1755 1754 1753 1748 1747 1746 1742 1727 1723 1720 1546 1127  979  583  396  385  384  323  416  418  405  451   26
[28]  420  380  398  387  386  424  383  382  286  393  307  381  377  402  500  379  378  156  376  422  375  374  373  372  371  370  369
[55]  368  367  366  365  364  363  362  361  360  359  358  357  356  355  354  351  350  349  348  347  346  345  344  343  341  340  338
[82]  337  336  334  333  331  327  322  320  318  314  313  308  305  304  303  301  300  299  298  296  295  294  289  281  280  279  278
[109]  277  274  263  260  255  254  253  251  248  243  283  231  225  223  220  214  212  211  200  199  192  183  182  181  180  179  174
[136]  172  159  155  153  147  143  142  140  132  131  130  127  124  122  121  119  110  107  106  105  103  102   95   81   78   77   75
[163]   72   55   46   17

unique(result.join.count2.odd$end)[order(unique(result.join.count2.odd$end))]
[1] 98839538 98839674 98839800 98839821 98839835 98839892 98839896 98839898 98839899 98839907 98839910 98839917 98839928 98839930 98839932
[16] 98839936 98839940 98839942 98839963 98840015

table(as.factor(result.join.count2.odd$end))
98839538 98839674 98839800 98839821 98839835 98839892 98839896 98839898 98839899 98839907 98839910 98839917 98839928 98839930 98839932 
1        1        1        1        1        2     6273        2        1        1        1        3        2        1        1 
98839936 98839940 98839942 98839963 98840015 
1        1        1        1        1 

##filter by end 98839896
result.join.count2.odd.9896<-result.join.count2.odd[result.join.count2.odd$end==98839896,]
unique(result.join.count2.odd.9896$width)
[1] 1763 1762 1757 1756 1755 1754 1753 1748 1747 1746 1742 1727 1723 1720 1546 1127  979  583  396  385  384  383  382  381  380  379  378
[28]  377  376  375  374  373  372  371  370  369  368  367  366  365  364  363  362  361  360  359  358  357  356  355  354  351  350  349
[55]  348  347  346  345  344  343  341  340  338  337  336  334  333  331  327  323  322  320  318  314  313  308  305  304  303  301  300
[82]  299  298  296  295  294  289  281  280  279  278  277  274  263  260  255  254  253  251  248  243  231  225  223  220  214  212  211
[109]  200  199  192  183  182  181  180  179  174  172  159  155  153  147  143  142  140  132  131  130  127  124  122  121  119  110  107
[136]  106  105  103  102   95   81   78   77   75   72   55   46   17

table(as.factor(result.join.count2.odd.9896$width))

17   46   55   72   75   77   78   81   95  102  103  105  106  107  110  119  121  122  124  127  130  131  132  140  142  143  147  153 
1    2    1    1    1    1    1    1    1    1    1    1    2    2    1    2    1    7    1    2    4   35    2    1    1    1    1    1 
155  159  172  174  179  180  181  182  183  192  199  200  211  212  214  220  223  225  231  243  248  251  253  254  255  260  263  274 
1    1    1    1    1    1    1    5    1    2    1    1    1    1    2    1    1    2    1    2    1    2    2    1    1    1    1    4 
277  278  279  280  281  289  294  295  296  298  299  300  301  303  304  305  308  313  314  318  320  322  323  327  331  333  334  336 
2    1    4    4   10    2    4    2    3   11    4    4   29    1    1    4    2    1    2    1   10    1    1    1    1    1    5    1 
337  338  340  341  343  344  345  346  347  348  349  350  351  354  355  356  357  358  359  360  361  362  363  364  365  366  367  368 
6    1    7    1   17   18    1    1    1    4   32    3    1   11    3    3    7    1    3    7    6   42   13    2    3    4    4    4 
369  370  371  372  373  374  375  376  377  378  379  380  381  382  383  384  385  396  583  979 1127 1546 1720 1723 1727 1742 1746 1747 
10   17   55   32   27   34   80   41  202  263   89  115 1195  905  124 2470   84    1    1    1    1    1    1    1    1    1    1    2 
1748 1753 1754 1755 1756 1757 1762 1763 
1    3    4    2    4    1   65    1 


##
result.join.count2.even<-result.join.count2 %>% dplyr::filter(row_number() %% 2 == 0) ## Select even rows
table(as.factor(result.join.count2.even$start))

98839655 98843309 98843658 98843665 98843666 98843668 98843669 98843670 98843671 98843673 98843675 98843679 98843701 98843736 98843752 
1       83        1        3        1        1        3        3       59       12     6124        1        2        1        1 
98843753 
1 

##join odd and even for exon2 (count=2) df##
result.join.count2.new<-full_join(result.join.count2.odd.9896,result.join.count2.even,by=".id")
result.join.count2.new<-na.omit(result.join.count2.new)
#based on above frequency make groups##
for (i in 1:nrow(result.join.count2.new)){
  if (result.join.count2.new$width.x[i]>900)result.join.count2.new$group[i]<-"isoform1"
  if (result.join.count2.new$width.x[i]<900&result.join.count2.new$start.y[i]==98843309)result.join.count2.new$group[i]<-"isoform2"
  if (result.join.count2.new$width.x[i]<900&result.join.count2.new$start.y[i]!=98843309)result.join.count2.new$group[i]<-"isoform3"
}

write.table(result.join.count2.new,"cell1_VCaP_pacbio_2exon2_isoform.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
##############################2 exons isoform#########################################
#isoform1 92     98838135-98839896   98843675-98843915                             ###
#isoform2 83     98839516-98839896   98843309-98843404                             ###
#isoform3 6099   98839513-98839896   98843675-98843915                             ###
######################################################################################

result.join.count3<-result.join[result.join$count==3,]

result.join.count3.first<-result.join.count3%>%dplyr::group_by(.id)%>%dplyr::filter(row_number() == 1) ## Select first rows
result.join.count3.second<-result.join.count3%>%dplyr::group_by(.id)%>%dplyr::filter(row_number() == 2) ## Select second rows
result.join.count3.third<-result.join.count3%>%dplyr::group_by(.id)%>%dplyr::filter(row_number() == 3) ##select third rows

result.join.count3.2<-full_join(result.join.count3.first,result.join.count3.second,by=".id")
result.join.count3.rejoin<-full_join(result.join.count3.2,result.join.count3.third,by=".id")

result.join.count3.rejoin$format<-paste(result.join.count3.rejoin$end.x,":",result.join.count3.rejoin$start.y,"-",result.join.count3.rejoin$end.y,":",result.join.count3.rejoin$start)

for (i in 1:nrow(result.join.count3.rejoin)){
  if(result.join.count3.rejoin$end.x[i] < 98828825) result.join.count3.rejoin$group[i]<-"isoform4"
  if(98831559<result.join.count3.rejoin$end.x[i]&result.join.count3.rejoin$end.x[i] < 98831970) result.join.count3.rejoin$group[i]<-"isoform5"
  if(98832070<result.join.count3.rejoin$end.x[i]&result.join.count3.rejoin$end.x[i] < 98836360) result.join.count3.rejoin$group[i]<-"isoform6"
  if(98838270<result.join.count3.rejoin$end.x[i]&result.join.count3.rejoin$end.x[i] < 98838285) result.join.count3.rejoin$group[i]<-"isoform7"
  if(98838520<result.join.count3.rejoin$end.x[i]&result.join.count3.rejoin$end.x[i] < 98839791) result.join.count3.rejoin$group[i]<-"isoform8"
  if(98839895<result.join.count3.rejoin$end.x[i]) result.join.count3.rejoin$group[i]<-"isoform9"
  
}
write.table(result.join.count3.rejoin,"cell1.VCaP.pacbio.3exons.isoform.txt",quote=FALSE,col.names = TRUE,row.names = FALSE,sep = "\t")
###3exons isoforms###
##isoform4: (4)    98828299 : 98839513 - 98839896 : 98843675 (not very uniform)
##isoform5: (1330) 98831563 : 98839791 - 98839896 : 98843675
##isoform6: (7)    98836357 : 98839513 - 98839896 : 98843675
##isoform7: (174)  98838281 : 98839791 - 98839896 : 98843675
##isoform8: (8)    98839549 : 98839644 - 98839896 : 98843675 (not very uniform)
##isoform9: (90)   98839896 : 98843309 - 98843479 : 98843675
#####################

result.join.count4<-filter(result.join,count==4)
result.join.count5<-filter(result.join,count==5)

result.join.count4.first<-result.join.count4%>%dplyr::group_by(.id)%>%dplyr::filter(row_number() == 1) ## Select first rows
result.join.count4.second<-result.join.count4%>%dplyr::group_by(.id)%>%dplyr::filter(row_number() == 2) ## Select second rows
result.join.count4.third<-result.join.count4%>%dplyr::group_by(.id)%>%dplyr::filter(row_number() == 3) ##select third rows
result.join.count4.forth<-result.join.count4%>%dplyr::group_by(.id)%>%dplyr::filter(row_number() == 4) ##select forth rows

result.join.count4.2<-full_join(result.join.count4.first,result.join.count4.second,by=".id")
result.join.count4.3<-full_join(result.join.count4.2,result.join.count4.third,by=".id")
result.join.count4.rejoin<-full_join(result.join.count4.3,result.join.count4.forth,by=".id")

result.join.count4.rejoin$format<-with(result.join.count4.rejoin,paste(end.x,":",start.y,"-",end.y,":",start.x.x,"-",end.x.x,":",start.y.y))

for (i in 1:nrow(result.join.count4.rejoin)){
  if(result.join.count4.rejoin$format[i]=="98831563 : 98839791 - 98839896 : 98841731 - 98841946 : 98843675") result.join.count4.rejoin$group[i]<-"isoform10"
  else if(result.join.count4.rejoin$format[i]=="98831959 : 98835340 - 98835391 : 98839791 - 98839896 : 98843675") result.join.count4.rejoin$group[i]<-"isoform11"
  else result.join.count4.rejoin$group[i]<-NA
}

write.table(result.join.count4.rejoin,"cell1.VCaP.pacbio.4exons.isoform.txt",quote=FALSE,col.names = TRUE,row.names = FALSE,sep = "\t")
#isoform10 (7): 98831028-98831563 98839791-98839896 98841731-98841946 98843675-98843915
#isoform11 (5): 98831672-98831959 98835340-98835391 98839791-98839896 98843675-98843915

##plot##
##cell1 VCaP##
library(ggplot2)
x<-c("isoform 3","isoform 5","isoform 7","isoform 1","isoform 9","isoform 2")
y<-c(6099,1330,174,92,90,83)
coord<-c("chr15:98839513-98843915","chr15:98831016-98843915","chr15:98838135-98843915","chr15:98838135-98843915","chr15:98839513-98843915","chr15:98839516-98843404")
cell1.vcap<-data.frame(isoform=x,count=y)
cell1.vcap$isoform<-factor(cell1.vcap$isoform,levels=c("isoform 2","isoform 9","isoform 1","isoform 7","isoform 5","isoform 3"))
#levels are in the order of bottom to top,smallest to largest##
p<-ggplot(data=cell1.vcap, aes(x=isoform, y=count)) +
  geom_bar(stat="identity")+theme(axis.text.x = element_text(angle =90, hjust = 1))+
  geom_text(aes(label=coord), position=position_dodge(width=0.9), vjust=12.25)
p

# Horizontal bar plot
p1<-p + coord_flip()
p2<-p1+theme(axis.text.x = element_text(angle =45, hjust = 1))+labs(x="cell1.VCaP")






