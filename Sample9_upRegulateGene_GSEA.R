### Cut-Tag diff gene GSEA 
### R 3.3.3
install.packages("dplyr")
library("dplyr")
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("DOSE")
library("DOSE")
require(DOSE)
BiocManager::install("clusterProfiler")
library(clusterProfiler)
require(clusterProfiler)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
require(org.Hs.eg.db)
BiocManager::install("org.Hs.eg.db")
#require(enrichplot)
library("ggplot2")
install.packages(c("Rmisc","lattice","plyr"))
library(Rmisc)
library(lattice)
library(plyr)
my_fun <- function(i) {
  strsplit(i,"-")[[1]][1]
}

Diff_gene_table_Sample9=read.csv("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/diff_peak_GO/Sample9vs10/uniquePeak_in_Sample9_with_geneName.csv")
Diff_gene_table_Sample9=Diff_gene_table_Sample9[,-c(1,7,8,9,10,11,12,13,14,16)]
min(Diff_gene_table_Sample9$V5)
Diff_gene_table_Sample9[,6]=apply(as.matrix(Diff_gene_table_Sample9[,6]),1,my_fun)
length(unique(Diff_gene_table_Sample9[,6]))
Diff_gene_table_Sample9=distinct(Diff_gene_table_Sample9)
length(unique(Diff_gene_table_Sample9[,6]))
Diff_gene_table_Sample9=Diff_gene_table_Sample9[order(Diff_gene_table_Sample9$V5,decreasing = T),]
colnames(Diff_gene_table_Sample9)=c("Chr","Start","End","info","log10 likehood ratio","SYMBOL")
###Unique Gene in Sample to Enrichment analysis
Gene_list=Diff_gene_table_Sample9$SYMBOL
ids <- bitr(Gene_list,fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
head(ids)
genes = ids[,2]
head(genes)
#GO BP MF KEGG
ekk <- enrichKEGG(gene = genes,
                  organism = 'hsa',
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.01)
hh <- as.data.frame(ekk)
hh=hh[c(1,4,18,
        2,14,19,21,
        3,9,13,16,17,
        5,6,7,8,10,11,12,15,20),]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
ggplot(hh,aes(y=order,x=Count))+
  geom_point(aes(size=Count,color=-1*p.adjust,))+# 修改点的大小
  scale_color_gradient(low="#f5e6e8",high = "#b7094c")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()

hh$Type=c(rep("Protein homeostasis", 3), 
          rep("Cell behavior", 4),
          rep("Signaling transduction", 5),
          rep("Infection and disease", 9))
          
COLS <- c("#2E94B9", "#D25565", "#F0B775", "#78c0a8")

ggplot(data=hh, aes(x=Count,y=order, fill=Type)) + 
  geom_point(aes(size=-1*p.adjust,color=Type)) + 
  scale_color_manual(values=COLS)+
  theme_bw() + 
  xlab("Gene Number") + 
  ylab("Pathways") + 
  labs(title = "The Most Enriched KEGG Pathways")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle是坐标轴字体倾斜的角度，可以自己设置

p1 <- ggplot(data=hh[1:3,], aes(x=Count,y=order, fill=Type)) + 
  geom_point(aes(size=-1*p.adjust,color=Type),stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  xlab("Gene Number") + 
  ylab("Pathways") + 
  labs(title = "The Most Enriched KEGG Pathways")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle是坐标轴字体倾斜的角度，可以自己设置
p2 <- ggplot(data=hh[4:7,], aes(x=Count,y=order, fill=Type)) + 
  geom_point(aes(size=-1*p.adjust,color=Type),stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  xlab("Gene Number") + 
  ylab("Pathways") + 
  labs(title = "The Most Enriched KEGG Pathways")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle是坐标轴字体倾斜的角度，可以自己设置
p3 <- ggplot(data=hh[8:12,], aes(x=Count,y=order, fill=Type)) + 
  geom_point(aes(size=-1*p.adjust,color=Type),stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  xlab("Gene Number") + 
  ylab("Pathways") + 
  labs(title = "The Most Enriched KEGG Pathways")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle是坐标轴字体倾斜的角度，可以自己设置
p4 <- ggplot(data=hh[13:21,], aes(x=Count,y=order, fill=Type)) + 
  geom_point(aes(size=-1*p.adjust,color=Type),stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  xlab("Gene Number") + 
  ylab("Pathways") + 
  labs(title = "The Most Enriched KEGG Pathways")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle是坐标轴字体倾斜的角度，可以自己设置

multiplot(p1,p2,p3,p4,cols = 1)

###GSEA analysis of ranking gene list:
my_fun <- function(i) {
  name=strsplit(i,"gene_name=")[[1]][2]
  name=strsplit(name,";")[[1]][1]
  name=strsplit(name,"-")[[1]][1]
}
common_gene_table=read.csv("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/diff_peak_GO/Sample9vs10/commonPeak_inSample9_10.csv")
common_gene_table=common_gene_table[,-c(6,7,8,9,10,11,12,13,15)]
min(common_gene_table$V5)
max(common_gene_table$V5)
common_gene_table[,6]=apply(as.matrix(common_gene_table[,6]),1,my_fun)
length(unique(common_gene_table[,6]))
common_gene_table=distinct(common_gene_table)
length(unique(common_gene_table[,6]))
common_gene_table=common_gene_table[order(common_gene_table$V5,decreasing = T),]
colnames(common_gene_table)=c("Chr","Start","End","info","log10 likehood ratio","SYMBOL")

Diff_gene_table_Sample10=read.csv("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/diff_peak_GO/Sample9vs10/uniquePeak_in_Sample10_with_geneName.csv",stringsAsFactors = F)
Diff_gene_table_Sample10=Diff_gene_table_Sample10[,-c(6,7,8,9,10,11,12,13,15)]
min(Diff_gene_table_Sample10$V5)
Diff_gene_table_Sample10$V5=Diff_gene_table_Sample10$V5*-1
Diff_gene_table_Sample10[,6]=apply(as.matrix(Diff_gene_table_Sample10[,6]),1,my_fun)
length(unique(Diff_gene_table_Sample10[,6]))
Diff_gene_table_Sample10=distinct(Diff_gene_table_Sample10)
length(unique(Diff_gene_table_Sample10[,6]))
Diff_gene_table_Sample10=Diff_gene_table_Sample10[order(Diff_gene_table_Sample10$V5,decreasing = T),]
colnames(Diff_gene_table_Sample10)=c("Chr","Start","End","info","log10 likehood ratio","SYMBOL")

RankingTable=rbind(Diff_gene_table_Sample9,common_gene_table,Diff_gene_table_Sample10)
length(unique(RankingTable$SYMBOL))
RankingTable=distinct(RankingTable)
###Unique Gene in Sample to Enrichment analysis
RankingTable$sumvalue=RankingTable$`log10 likehood ratio`
RankingTable <- RankingTable %>%
  group_by(SYMBOL) %>%
  filter(sumvalue==max(sumvalue))
ids <- bitr(Gene_list,fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
geneList=merge(RankingTable,ids)
geneList=geneList[order(geneList$`log10 likehood ratio`,decreasing = T),]
genes=geneList$`log10 likehood ratio`
names(genes)=geneList$ENTREZID
kegg_result <- gseKEGG(
  genes,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = FALSE,
  by = "fgsea")
gseKEGG_result<-as.data.frame(kegg_result@result)
write.csv(gseKEGG_result,"/Users/tangsuiyue/Documents/Documents/Project/FBXO42/diff_peak_GO/diff_gene_GSEA/GSEA.csv",row.names = F,quote = F)
gseaplot2(kegg_result,"hsa04330",title = "Notch signaling pathway")
gseaplot2(kegg_result,"hsa04350",title = "TGF-beta signaling pathway")
gseaplot2(kegg_result,"hsa05221",title = "Acute myeloid leukemia")
gseaplot2(kegg_result,"hsa03460",title = "Fanconi anemia pathway")


gseaplot2(kegg_result,1,"Glycerolipid metabolism")
dotplot(kegg_result,title = "Enrichment_KEGG")

