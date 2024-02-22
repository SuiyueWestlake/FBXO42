if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
###切换4.1.1
BiocManager::install("parallel")
library(parallel)
BiocManager::install("BiocGenerics")
library("BiocGenerics")

###切换3.6.3
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

BiocManager::install("ChIPseeker")
library("ChIPseeker")
BiocManager::install("GenomicRanges")
library("GenomicRanges") 
BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
BiocManager::install("ChIPpeakAnno")
#install.packages("ChIPpeakAnno")
library("ChIPpeakAnno")

Sample1 <- readPeakFile("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/20210825CutTagData/Sample1vs7/Sample1vs7_peaks.bed")
Sample2 <- readPeakFile("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/20210825CutTagData/Sample2vs8/Sample2vs8_peaks.bed")
Sample3 <- readPeakFile("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/20210825CutTagData/Sample3vs7/Sample3vs7_peaks.bed")
Sample4 <- readPeakFile("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/20210825CutTagData/Sample4vs8/Sample4vs8_peaks.bed")
Sample5 <- readPeakFile("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/20210825CutTagData/Sample5vs7/Sample5vs7_peaks.bed")
Sample6 <- readPeakFile("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/20210825CutTagData/Sample6vs8/Sample6vs8_peaks.bed")
Sample9 <- readPeakFile("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/20210825CutTagData/Sample9vs7/Sample9vs7_peaks.bed")
Sample10 <- readPeakFile("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/20210825CutTagData/Sample10vs8/Sample10vs8_peaks.bed")

#annotation
anno_sample9 <- annotatePeak(Sample9, tssRegion = c(-2000, 1000), TxDb = txdb, annoDb="org.Hs.eg.db", level="gene")
anno9_new=as.data.frame(anno_sample9)
my_fun <- function(i) {
  strsplit(i," ")[[1]][1]
}
anno9_new$annotation=apply(as.matrix(anno9_new$annotation),1,my_fun)

Table9=as.matrix(table(anno9_new$annotation))
Table9=c(Table9[5:7,],sum(Table9[1:4,]))
names(Table9)[4]="Intergenic"
Table9=Table9[c("Promoter","Exon","Intron","Intergenic")]
par(mfrow=c(1,2))
piepercent9<- paste(round(100*Table9/sum(Table9), 2), "%")
pie(Table9,label=piepercent9,clockwise=T,col = c("#D25565","#F0B775","#2E94B9","#3b9a9c"),main="Sample9")
legend("topright",c("Promoter","Exon","Intron","Intergenic"),cex=0.6,#添加图例
       fill=c("#D25565","#F0B775","#2E94B9","#3b9a9c")) #调色板

###Sample10
anno_sample10 <- annotatePeak(Sample10, tssRegion = c(-2000, 1000), TxDb = txdb, annoDb="org.Hs.eg.db", level="gene")
anno10_new=as.data.frame(anno_sample10)
anno10_new$annotation=apply(as.matrix(anno10_new$annotation),1,my_fun)
Table10=as.matrix(table(anno10_new$annotation))
Table10=c(Table10[5:7,],sum(Table10[1:4,]))
names(Table10)[4]="Intergenic"
Table10=Table10[c("Promoter","Exon","Intron","Intergenic")]
piepercent10<- paste(round(100*Table10/sum(Table10), 2), "%")
pie(Table10,label=piepercent10,clockwise=T,col = c("#D25565","#F0B775","#2E94B9","#3b9a9c"),main="Sample10")
legend("topright",c("Promoter","Exon","Intron","Intergenic"),cex=0.6,#添加图例
       fill=c("#D25565","#F0B775","#2E94B9","#3b9a9c")) #调色板

#Sample5
anno_sample5 <- annotatePeak(Sample5, tssRegion = c(-2000, 1000), TxDb = txdb, annoDb="org.Hs.eg.db", level="gene")
anno5_new=as.data.frame(anno_sample5)
anno5_new$annotation=apply(as.matrix(anno5_new$annotation),1,my_fun)
Table5=as.matrix(table(anno5_new$annotation))
Table5=c(Table5[5:7,],sum(Table5[1:4,]))
names(Table5)[4]="Intergenic"
Table5=Table5[c("Promoter","Exon","Intron","Intergenic")]
piepercent5<- paste(round(100*Table5/sum(Table5), 2), "%")
pie(Table5,label=piepercent5,clockwise=T,col = c("#D25565","#F0B775","#2E94B9","#3b9a9c"),main="Sample5")
legend("topright",c("Promoter","Exon","Intron","Intergenic"),cex=0.6,#添加图例
       fill=c("#D25565","#F0B775","#2E94B9","#3b9a9c")) #调色板

anno_sample6 <- annotatePeak(Sample6, tssRegion = c(-2000, 1000), TxDb = txdb, annoDb="org.Hs.eg.db", level="gene")
anno6_new=as.data.frame(anno_sample6)
anno6_new$annotation=apply(as.matrix(anno6_new$annotation),1,my_fun)
Table6=as.matrix(table(anno6_new$annotation))
Table6=c(Table6[5:7,],sum(Table6[1:4,]))
names(Table6)[4]="Intergenic"
Table6=Table6[c("Promoter","Exon","Intron","Intergenic")]
piepercent6<- paste(round(100*Table6/sum(Table6), 2), "%")
pie(Table6,label=piepercent6,clockwise=T,col = c("#D25565","#F0B775","#2E94B9","#3b9a9c"),main="Sample6")
legend("topright",c("Promoter","Exon","Intron","Intergenic"),cex=0.6,#添加图例
       fill=c("#D25565","#F0B775","#2E94B9","#3b9a9c")) #调色板


###两个样本之间peak overlap分析
Sample1vs9 <- findOverlapsOfPeaks(Sample1,Sample9)
Sample1vs10 <- findOverlapsOfPeaks(Sample1,Sample10)
Sample9vs10 <- findOverlapsOfPeaks(Sample9,Sample10)
Sample5vs6 <- findOverlapsOfPeaks(Sample5,Sample6)
Sample3vs9 <- findOverlapsOfPeaks(Sample3,Sample9)
Sample3vs10 <- findOverlapsOfPeaks(Sample3,Sample10)

Sample1vs2 <- findOverlapsOfPeaks(Sample1,Sample2)
Sample3vs4 <- findOverlapsOfPeaks(Sample3,Sample4)
Sample1vs3 <- findOverlapsOfPeaks(Sample1,Sample3)
Sample2vs4 <- findOverlapsOfPeaks(Sample2,Sample4)

# 绘制venn图
par(side = 3)
makeVennDiagram(Sample1vs9)
makeVennDiagram(Sample1vs10)
makeVennDiagram(Sample9vs10)
makeVennDiagram(Sample5vs6)
makeVennDiagram(Sample3vs9)
makeVennDiagram(Sample3vs10)

makeVennDiagram(Sample1vs2)
makeVennDiagram(Sample3vs4)
makeVennDiagram(Sample1vs3)
makeVennDiagram(Sample2vs4)

###分析overlap peak
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
Sample5vs6
gr <- GRanges(Sample9, format="BED", header=FALSE)
seq <- getAllPeakSequence(gr, upstream=20, downstream=20, genome=Hsapiens)
write2FASTA(seq, "sampleA.peaks.fa")

# 用1号染色体的碱基分布当做背景
freqs <- oligoFrequency(Hsapiens$chr1, MarkovOrder=3)
# oligoLength规定了motif的长度
os <- oligoSummary(seq, oligoLength=7, MarkovOrder=3, quickMotif=TRUE, freqs=freqs)
zscore <- sort(os$zscore)
# 绘制所有6个碱基组合的频率分布图
h <- hist(zscore, breaks=100, xlim=c(-50, 50), main="Histogram of Z-score")
# 频率最大的碱基组合即为motif的结果
text(zscore[length(zscore)], max(h$counts)/10,     labels=names(zscore[length(zscore)]), adj=1)

BiocManager::install("motifStack")
library(motifStack)
pfms <- mapply(function(.ele, id)    new("pfm", mat=.ele, name=paste("SAMPLE motif", id)),    os$motifs, 1:length(os$motifs))
motifStack(pfms[[1]])
plot(pfms[[1]])

# 准备基因组注释信息
BiocManager::install("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)
annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")
# 进行
test_peak=Sample5vs6$peaklist
test_peak=test_peak$`Sample5///Sample6`
test_peak=as.data.frame(test_peak)
test_peak <- toGRanges(test_peak)
overlaps.anno <- annotatePeakInBatch(test_peak,AnnotationData=annoData,output="nearestLocation")
library(org.Hs.eg.db)
overlaps.anno <- addGeneIDs(overlaps.anno, "org.Hs.eg.db", IDs2Add = "entrez_id")
Genes=overlaps.anno$entrez_id
Genes=na.omit(Genes)
write.table(Genes, "/Users/tangsuiyue/Documents/Documents/Project/FBXO42/diff_peak_GO/diff_genes.txt",quote=F,col.names=F,row.names = F)

BiocManager::install("clusterProfiler")
library(clusterProfiler)
ekk <- enrichKEGG(gene = Genes,
                  organism = 'hsa',
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)
ekk <- setReadable(ekk,OrgDb = org.Hs.eg.db,keytype = "ENTREZID")

barplot(ekk,showCategory = 15,title = "EnrichmentGO_BP")
dotplot(ekk,showCategory = 15,title = "EnrichmentGO_BP")


pie1(table(overlaps.anno$insideFeature))
over <- getEnrichedGO(overlaps.anno, orgAnn="org.Hs.eg.db", maxP=.05, minGOterm=10, multiAdjMethod="BH", condense=TRUE)
path<-getEnrichedPATH(overlaps.anno,"org.Hs.eg.db","reactome.db",maxP=.05)
class(over[[1]])
barplot(over[[1]]$bp,showCategory = 15,title = "Enrichment_KEGG")


###以下不用运行
#Pie and Bar plot
plotAnnoPie(anno10_new)

#full annotation overlap
install.packages("ggupset")
library("ggupset")
upsetplot(anno)

###多组比较
peaks <- list(Sample9=Sample9,Sample10=Sample10)
# promotor区间范围可以自己设定
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
#annotatePeak传入annoDb参数,可进行基因ID转换（Entrez，ENSEMBL，SYMBOL，GENENAME）
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList,title="Distribution of transcription factor-binding loci \n relative to TSS")

