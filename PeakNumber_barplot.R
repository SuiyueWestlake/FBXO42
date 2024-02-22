###ATAC seq peaks 
WT=read.csv("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/20220602_ATACseq_Data/Peak_R_analysis/293T_WT_peaks.xls")
FBXO42=read.csv("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/20220602_ATACseq_Data/Peak_R_analysis/293T_FBXO42KO_peaks.xls")
PeakNumber=c(dim(WT)[1],dim(FBXO42)[1])
names(PeakNumber)=c("WT","FBXO42")
barplot(PeakNumber,ylim = c(0,35000))
