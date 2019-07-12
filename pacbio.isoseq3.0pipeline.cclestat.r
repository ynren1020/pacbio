##################################2019-07-01###################################
##lbc008-A549;lbc016-H226;lbc024-DMS-53;lbc032-H727;lbc040-VCaP;lbc048-LNCaP###
##no igf1r-as1 isoform: lbc016-H226, lbc048-LNCaP##############################
###############################################################################

isoform_normlize<-function(nreads,totalFL){
  nreads.norm<-nreads/totalFL*10^6
}

#VCaPcell1##
nreads<-c(2,18,24,15,3,3,2,233)
totalFL<-10659
VCaP.cell1<-data.frame(isoform=paste0("PB.53.",1:8),NormCount=isoform_normlize(nreads,totalFL))
#H727cell1##
nreads<-c(18,2,2,3,199,39,2,11,233)
totalFL<-8536
H727.cell1<-data.frame(isoform=paste0("PB.33.",1:9),NormCount=isoform_normlize(nreads,totalFL))
#DMS53#
nreads<-c(3,7,2,24)
totalFL<-12573
DMS53.cell1<-data.frame(isoform=paste0("PB.58.",1:4),NormCount=isoform_normlize(nreads,totalFL))
#A549##
nreads<-c(13,11,35)
totalFL<-10101
A549.cell1<-data.frame(isoform=paste0("PB.45.",1:3),NormCount=isoform_normlize(nreads,totalFL))


##merged cells##
#VCaPcell1##
nreads<-c(2,32,4,51,21,8,3,3,435)
totalFL<-19296
VCaP.merge<-data.frame(isoform=paste0("PB.69.",1:9),NormCount=isoform_normlize(nreads,totalFL))
#H727cell1##
nreads<-c(32,2,6,364,66,3,3,2,20,412)
totalFL<-15425
H727.merge<-data.frame(isoform=paste0("PB.40.",c(1:3,6:12)),NormCount=isoform_normlize(nreads,totalFL))
#DMS53#
nreads<-c(5,9,2,5,45)
totalFL<-22852
DMS53.merge<-data.frame(isoform=paste0("PB.65.",1:5),NormCount=isoform_normlize(nreads,totalFL))
#A549##
nreads<-c(28,2,16,4,2,65)
totalFL<-18098
A549.merge<-data.frame(isoform=paste0("PB.57.",1:6),NormCount=isoform_normlize(nreads,totalFL))

##output data##
write.table(VCaP.merge,"VCaP.merge.isoform.txt",quote=FALSE,col.names = TRUE,row.names = FALSE,sep = "\t")
write.table(H727.merge,"H727.merge.isoform.txt",quote=FALSE,col.names = TRUE,row.names = FALSE,sep = "\t")
write.table(DMS53.merge,"DMS53.merge.isoform.txt",quote=FALSE,col.names = TRUE,row.names = FALSE,sep = "\t")
write.table(A549.merge,"A549.merge.isoform.txt",quote=FALSE,col.names = TRUE,row.names = FALSE,sep = "\t")

###2019-07-12###
##rnaseq reads/isoform quantified by Salmon##
##/home/yren/project/pacbio/rnaseq/transcripts_quant/quant.sf##
tpm<-c(44712.7,0,222866.4,131981.0,206014.8,0,0,0,394425.0)
VCaP.rnaseqSalmon<-data.frame(isoform=paste0("PB.69.",1:9),NormCount=tpm)

tpm<-c(171400.7,0,0,443558.0,96349.2,0,93679.8,0,195012.2,0)
H727.rnaseqSalmon<-data.frame(isoform=paste0("PB.40.",c(1:3,6:12)),NormCount=tpm)

##plot##
##cell1 VCaP##
library(ggplot2)
#x<-c("isoform 3","isoform 5","isoform 7","isoform 1","isoform 9","isoform 2")
#y<-c(6099,1330,174,92,90,83)
#coord<-c("chr15:98839513-98843915","chr15:98831016-98843915","chr15:98838135-98843915","chr15:98838135-98843915","chr15:98839513-98843915","chr15:98839516-98843404")
#cell1.vcap<-data.frame(isoform=x,count=y)
#cell1.vcap$isoform<-factor(cell1.vcap$isoform,levels=c("isoform 2","isoform 9","isoform 1","isoform 7","isoform 5","isoform 3"))
#levels are in the order of bottom to top,smallest to largest##

#VCaP.merge$isoform<-factor(VCaP.merge$isoform,levels=c("PB.69.1","PB.69.7","PB.69.8","PB.69.3","PB.69.6","PB.69.5","PB.69.2","PB.69.4","PB.69.9"))
#p<-ggplot(data=VCaP.merge, aes(x=isoform, y=NormCount)) +
#  geom_bar(stat="identity")+theme(axis.text.x = element_text(angle =90, hjust = 1))
 # Horizontal bar plot
#p1<-p + coord_flip()
#p2<-p1+theme(axis.text.x = element_text(angle =45, hjust = 1))+labs(x="VCaP Transcript Isoforms",y="Transcript Count per Million")
#p2


isoform_plot<-function(dat){
  name<-deparse(substitute(dat))
  cells<-strsplit(name,"[.]")[[1]][1]
  dat$isoform<-factor(dat$isoform,levels=dat$isoform[order(dat$NormCount)])
  p<-ggplot(data=dat,aes(x=isoform,y=NormCount))+
    geom_bar(stat="identity")+theme(axis.text.x = element_text(angle =90, hjust = 1))
  # Horizontal bar plot
  p1<-p + coord_flip()
  p2<-p1+theme(axis.text.x = element_text(angle =45, hjust = 1))+labs(x=paste0(cells,"Transcript Isoforms"),y="Transcript Count per Million")
  p2
  ##pacbio isoform plot##
  ggsave(filename=paste0(cells,".merge.isoform.png"),width=5,height = 5,units = "in",dpi = 300,device = "png")
  ##rnaseq isoform plot##
  #ggsave(filename=paste0(cells,".rnaseqsalmon.isoform.png"),width=5,height = 5,units = "in",dpi = 300,device = "png")
}

isoform_plot(VCaP.merge)
isoform_plot(A549.merge)
isoform_plot(H727.merge)
isoform_plot(DMS53.merge)
isoform_plot(H727.rnaseqSalmon)
