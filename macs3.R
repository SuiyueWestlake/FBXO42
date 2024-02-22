pkgs <- c(
  "ChIPseeker", "GenomicFeatures", "clusterProfiler", "glue", "openxlsx",
  "dplyr", "org.Hs.eg.db", "DOSE", "ReactomePA", "TxDb.Hsapiens.UCSC.hg38.knownGene"
)
for (pkg in pkgs) suppressPackageStartupMessages(library(pkg, character.only = T))

setwd("~/projects/dux4/analysis/mgy/chipseq/")
##### create txdb using the gtf file #####
# gtf <- "~/doc/reference/gtf/gencode.v32.annotation.gtf"
# txdb <- makeTxDbFromGFF(gtf, format = "gtf")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
##### define input experiments #####
experiments <- c("DUX4_IGH", "DUX4_IGH_RAG1_RAG2", "RAG1_RAG2_H3")
##### read in bed files and filter them #####
beds_peak <- glue("chilin/{experiments}/{experiments}_peaks.narrowPeak") %>%
  setNames(experiments) %>%
  as.list()
peaks <- lapply(beds_peak, readPeakFile)
peaks_fil <- lapply(peaks, function(x) {return(x[abs(x$V7) >= 4 & x$V9 >= 3])}) ### fold change >= 4, p <= 0.001
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrixList <- lapply(peaks_fil, getTagMatrix, weightCol = "V5", windows = promoter)
##### peak heatmap for all samples #####
pdf("figs/peakheatmap.pdf", height = 7, width = 6)
peakHeatmap(peaks_fil, weightCol = "V5", TxDb = txdb, upstream = 3000, downstream = 3000, color = rainbow(length(peaks_fil)))
dev.off()
##### average profile for all samples #####
pdf("figs/avgprof.pdf", height = 7, width = 8)
plotAvgProf2(peaks_fil, weightCol = "V5", TxDb = txdb, upstream = 3000, downstream = 3000)
dev.off()

pdf("figs/avgprof_facet.pdf", height = 6, width = 7)
plotAvgProf(tagMatrixList, xlim = c(-3000, 3000), facet = "row", conf = 0.95, resample = 1000)
dev.off()
##### peak annotation for all samples #####
peaks_anno <- lapply(peaks_fil, annotatePeak, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
for (experiment in experiments) {
  write.table(
    as.data.frame(peaks_anno[[experiment]]),
    glue("chilin/{experiment}/{experiment}.annot.tab"),
    quote = F, row.names = F, sep = "\t"
  )
}
pdf("figs/annotation_peaks.pdf", height = 5, width = 9)
plotAnnoBar(peaks_anno)
dev.off()
pdf("figs/annotation_peaks_distance.pdf", height = 5, width = 9)
plotDistToTSS(peaks_anno, title="Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()
##### enrichment for all samples #####
genes <- lapply(peaks_anno, function(x) {as.data.frame(x)$geneId}) #  %>% gsub("\\.\\d*", "", .)
# genes_ensembl <- lapply(genes, function(x) {return(gsub("\\..*", "", x))})
# genes_entrez <- lapply(peaks_anno, function(x) {as.data.frame(x)$ENTREZID %>% .[!is.na(.)]})
# lapply(genes, function(x) {bitr(x, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") %>% .$ENTREZID})
compareCluster(
  geneCluster = genes,
  fun = "enrichGO",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  keyType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  ont = "BP"
) -> enrich_go
pdf("figs/enrichment_comparison.pdf", height = 10, width = 11)
dotplot(enrich_go, showCategory = 25, title = "GO BP Enrichment Analysis")
dev.off()

genes2 <- lapply(peaks, function(x) {
  seq2gene(x, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb = txdb)
})
genes2 <- lapply(genes2, function(x) {
  return(x[!grepl("ENST", x)])
})
compareCluster(
  geneCluster = genes2,
  fun = "enrichGO",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  keyType = "ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "BP"
) -> enrich_go2
pdf("figs/enrichment_comparison2.pdf", height = 7, width = 8)
dotplot(enrich_go2, showCategory = 10, title = "GO BP Enrichment Analysis")
dev.off()
##### venn plot for genes #####
pdf("figs/venn_genes.pdf", height = 5, width = 5)
vennplot(genes)
dev.off()
# txdb2 <- TxDb.Hsapiens.UCSC.hg38.knownGene
# enrichPeakOverlap(
#   queryPeak = beds_peak[[1]],
#   targetPeak = unlist(beds_peak[2:3]),
#   TxDb = txdb2,
#   pAdjustMethod = "BH",
#   nShuffle = 1,
#   mc.cores = 1
# )
##### coverage plot for each sample #####
for (experiment in experiments) {
  peak <- peaks_fil[[experiment]]
  pdf(as.character(glue("figs/{experiment}.covplot.pdf")), height = 8, width = 6)
  print(
    covplot(peak, weightCol = "V5", lower = 20, chrs = paste0("chr", 1:22))
  )
  dev.off()
}


genes <- lapply(peaks_anno[c("DUX4_IGH", "DUX4_IGH_RAG1_RAG2")], function(x) {as.data.frame(x)$geneId})
enrich_go <- compareCluster(
  geneCluster = genes,
  fun = "enrichGO",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  keyType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  ont = "BP"
)
pdf("figs/enrichment/go.pdf", width = 8, height = 10)
dotplot(enrich_go, title = "GO BP Enrichment Analysis", showCategory = 30)
dev.off()
write.xlsx(as.data.frame(enrich_go), "enrichment/go.xlsx")
enrich_kegg <- compareCluster(
  geneCluster = genes,
  fun = "enrichKEGG",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  # qvalueCutoff = 1,
  organism = "hsa",
  keyType = "kegg"
)
pdf("figs/enrichment/kegg.pdf", width = 8, height = 9)
dotplot(enrich_kegg, title = "KEGG Enrichment Analysis", showCategory = 30)
dev.off()
write.xlsx(as.data.frame(enrich_kegg), "enrichment/kegg.xlsx")
enrich_mkegg <- compareCluster(
  geneCluster = genes,
  fun = "enrichMKEGG",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  organism = "hsa",
  keyType = "kegg"
)
pdf("figs/enrichment/mkegg.pdf", width = 10, height = 10)
dotplot(enrich_mkegg, title = "MKEGG Enrichment Analysis", showCategory = 30)
dev.off()
write.xlsx(as.data.frame(enrich_mkegg), "enrichment/mkegg.xlsx")
enrich_wp <- compareCluster(
  geneCluster = genes,
  fun = "enrichWP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  organism = "Homo sapiens"
)
pdf("figs/enrichment/wp.pdf", width = 9, height = 10)
dotplot(enrich_wp, title = "WikiPathways Enrichment Analysis", showCategory = 30)
dev.off()
write.xlsx(as.data.frame(enrich_wp), "enrichment/wp.xlsx")
enrich_reactome <- compareCluster(
  geneCluster = genes,
  fun = "enrichPathway",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  organism = "human"
)
pdf("figs/enrichment/reactome.pdf", width = 9, height = 10)
dotplot(enrich_reactome, title = "Reactome Enrichment Analysis", showCategory = 30)
dev.off()
write.xlsx(as.data.frame(enrich_reactome), "enrichment/reactome.xlsx")


myEnrichPlot <- function(data, title) {
  data <- data[order(data$GeneRatio, decreasing = F), ]
  data$Description <- factor(data$Description, levels = data$Description)
  data <- data[order(data$GeneRatio, decreasing = T), ]
  ggplot(data, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = .5)) +
    scale_colour_gradient(limits = range(data$p.adjust), low = "red") +
    ylab(NULL) +
    ggtitle(title)
}

enrichment_selected <- read.xlsx("enrichment/interested.xlsx")

pdf("figs/enrichment/interested.pdf", width = 15, height = 15)
myEnrichPlot(enrichment_selected, "enrichment analysis for up binding genes")
dev.off()







