###############2019-06-18###########################
#####SGSeq to predict splicing event################
##VCaP illumina RNAseq data#########################
##igf1r-as1 region bam##############################
##/home/yren/project/pacbio/rnaseq##################

##install and load packages##
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SGSeq")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(SGSeq)
df<-data.frame(sample_name="VCaP",file_bam="VCaP.igf1r-as1.bam",stringsAsFactors = FALSE)
test<-getBamInfo(df)
##create Grange object##
gr <- GRanges(seqnames = "chr15", strand ="-",
             ranges = IRanges(start = 98824244, width = 36449))

#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#txdb.chr15 <- keepSeqlevels(txdb, "chr15")
#seqlevelsStyle(txdb) <- "NCBI"

##make TxDb from GTF##
library(GenomicFeatures)
gtfFile <- "IGF1R-AS1.hg38.gtf"
txdb <- makeTxDbFromGFF(file=gtfFile,
                        dataSource="igf1r-as1",format = "gtf"
                        )

txf_ucsc <- convertToTxFeatures(txdb)
#txf_ucsc <- txf_ucsc[txf_ucsc %over% gr]
sgf_ucsc <- convertToSGFeatures(txf_ucsc)

##solve error in bam header##

#seqlevels_in_bam <-seqlevels(BamFile(test$file_bam[1]))
#sgfdb2 <- keepSeqlevels(sgf_ucsc, seqlevels_in_bam)
##
sgfc_ucsc <- analyzeFeatures(test, features = txf_ucsc)

colData(sgfc_ucsc)
rowRanges(sgfc_ucsc)
head(counts(sgfc_ucsc))
head(FPKM(sgfc_ucsc))

df <- plotFeatures(sgfc_ucsc,geneID = 1)

##Splice graph analysis based on de novo prediction##
##prediction##
sgfc_pred <- analyzeFeatures(test, which = gr)
head(rowRanges(sgfc_pred))
##predicted features can be annotated with respect to known transcripts##
sgfc_pred <- annotate(sgfc_pred, txf_ucsc)
head(rowRanges(sgfc_pred))

df <- plotFeatures(sgfc_pred, color_novel = "red",heightPanels = c(10,7),square = TRUE,geneID = 1)


par(mfrow = c(2, 1), mar = c(1, 3, 1, 1))
plotSpliceGraph(rowRanges(sgfc_pred),  toscale = "none", color_novel = "red",ypos = c(0.3,0.1),geneID = 1)
plotCoverage(sgfc_pred,  toscale = "none",geneID = 1)


##regtools VCaP.igf1r-as1.bam junction##
regtools.out<-data.frame(chrom="chr15",start.reg=c(98803362,98831563,98831964,98831959,98838281,98838281,98839896,98839896,98839896,98843479),
                         end.reg=c(98838731,98839791,98839791,98839791,98839791,98839434,98843675,98843309,98841731,98843675),score=c(1,16,2,6,7,1,298,10,1,12)
)

##prediction output of junction##
pred<-rowRanges(sgfc_pred)
pred<-as.data.frame(pred)

vcapN<-counts(sgfc_pred)
pred<-cbind(pred,vcapN)

pred.junction<-pred[pred$type=="J"&pred$geneID==1,]
regtools.out$start.reg<-as.integer(regtools.out$start.reg)
pred.junction.join<-dplyr::full_join(pred.junction,regtools.out,by=c("start"="start.reg","end"="end.reg"))
pred.junction.join<-pred.junction.join[,-c(7,8,11)]

write.table(pred.junction.join,"VCaP.junction.SGSeq.regtools.match.txt",quote = FALSE,col.names = TRUE,row.names = FALSE,sep = "\t")

