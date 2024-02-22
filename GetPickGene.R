###gtf2bed
hg38=read.table("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/peak2gene/gencode.v38.annotation.gff3",sep = "\t")
hg38=hg38[which(hg38[,3]=="gene"),]
hg38=hg38[-which(hg38[,1]=="chrM"),]
length(hg38$V9)
for (i in 1:dim(hg38)[1]) {
  #A=unlist(strsplit(hg38[i,9],split = "gene_name="))[2]
  #hg38[i,9]=unlist(strsplit(A,split = ";"))[1]
  #genelist=c(genelist,genename)
  if(hg38[i,7]=="+"){
    if(hg38[i,4]<1000){
      hg38[i,4]=0
    }
    else{
      hg38[i,4]=hg38[i,4]-1000
    }
  }
  if(hg38[i,7]=="-"){
    hg38[i,5]=hg38[i,5]+1000
  }
}
write.table(hg38,"/Users/tangsuiyue/Documents/Documents/Project/FBXO42/peak2gene/gencode.v38.annotation2.gff3",sep = "\t",col.names = F,row.names = F,quote = F)
###upload gencode.v38.annotation2.gff3 to CPU and run bedtools intersect

peak=read.table("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/peak2gene/slurm-704861.out")
length(unique(peak$V4))
length(unique(peak$V8))
peaks_in_gene=c()
for (i in 1:dim(peak)[1]) {
  A=unlist(strsplit(peak[i,14],split = "gene_name="))[2]
  genename=unlist(strsplit(A,split = ";"))[1]
  peaks_in_gene=c(peaks_in_gene,genename)
}
peak$genename=peaks_in_gene
length(unique(peak$genename))
HES4=hg38[grep(pattern="HES4",hg38[,4]),]
which(peak[,16]=="HEY2")
peak[which(peak[,16]=="HES4"),]

Result_toWeixiang=peak[,c(1,2,3,4,5,16,6,7,8,9,10,11,12,13,14,15)]
write.table(Result_toWeixiang,"/Users/tangsuiyue/Documents/Documents/Project/FBXO42/peak2gene/Result_toWeixiang.tet",sep = "\t",col.names = T,row.names = F,quote = F)
write.csv(Result_toWeixiang,"/Users/tangsuiyue/Documents/Documents/Project/FBXO42/peak2gene/Result2Weixiang.csv",quote = F)
Peak=read.csv("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/peak2gene/Result2Weixiang.csv")


