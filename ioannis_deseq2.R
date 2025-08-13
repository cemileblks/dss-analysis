library("DESeq2")
library(ggplot2)
library("pheatmap")
library("clusterProfiler")
library("DOSE")
library("org.Hs.eg.db")
library("cowplot")
library(EnhancedVolcano)
library("viridis") 
library(patchwork)
library(edgeR)
library("RColorBrewer")
library(readr)
library(tximport)
library(PCAtools)
library(ggrepel)
library(tidyverse)
library(biomaRt)
library(enrichplot)
library("ggnewscale")

setwd("/Users/pgautie2/Documents/Bioinfo/ioannis_dnmt3b")


salmon_table_ioannis<-read.table("salmon.merged.gene_counts_length_scaled.tsv",sep="\t",header=TRUE,row.names="gene_id")
cts_ioannis <- salmon_table_ioannis[,c(2:9)]
cts_ioannis<-as.matrix(round(cts_ioannis))
cts_ioannis
coldata_ioannis<-read.csv("metadata.csv",header=TRUE,row.names=1)
cts_ioannis <-cts_ioannis[,c(6:8,3:5,2,1)]
  
dds_ioannis<- DESeqDataSetFromMatrix(countData = cts_ioannis,
                                       colData = coldata_ioannis,
                                       design= ~ Genotype)
dds_ioannis <- dds_ioannis[ rowSums(counts(dds_ioannis)) > 1, ]
vsd_ioannis<-vst(dds_ioannis, blind=TRUE)
plotPCA(vsd_ioannis, intgroup=c("Genotype"))



#pcaData <- plotPCA(vsd_gast, intgroup=c("time","treatment"), returnData=TRUE)
#percentVar <- round(100 * attr(pcaData, "percentVar"))
#ggplot(pcaData, aes(PC1, PC2, color=time,shape=treatment)) +
#  geom_point(size=3) +
#  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#  coord_fixed() +  
#  geom_label_repel(aes(label = name), colour = I(alpha("black", 0.85)), size = 3,box.padding = 1.5,fill="white",max.overlaps = 30 );




dds_ioannis<-DESeq(dds_ioannis)
res_ioannis <- results(dds_ioannis)
summary(res_ioannis)

#### Relevel ESC #####

levels(dds_ioannis$Genotype)
dds_ioannis$Genotype<-relevel(dds_ioannis$Genotype,ref="WT")

dds_ioannis<-DESeq(dds_ioannis)

resultsNames(dds_ioannis)

res_Genotype_DNMT3B_vs_WT_noshrink<- results(dds_ioannis, contrast=c("Genotype","DNMT3B","WT"), independentFiltering=TRUE)
res_Genotype_DNMT3B_vs_WT <- lfcShrink(dds_ioannis, coef="Genotype_DNMT3B_vs_WT", type="apeglm",res=res_Genotype_DNMT3B_vs_WT_noshrink)
res_Genotype_DNMT3B_vs_WT_ordered  <- res_Genotype_DNMT3B_vs_WT[order(res_Genotype_DNMT3B_vs_WT$padj),]
#write.table(res_Genotype_DNMT3B_vs_WT_ordered,sep="\t",file="Deseq2_Gast_Genotype_DNMT3B_vs_WT_190623.txt")

#### OR annotation via package #####

res_Genotype_DNMT3B_vs_WT_ordered_sval <- lfcShrink(dds_ioannis, coef="Genotype_DNMT3B_vs_WT", type="apeglm",res=res_Genotype_DNMT3B_vs_WT_noshrink,svalue=TRUE)
res_Genotype_DNMT3B_vs_WT_ordered_merged<-merge(as.data.frame(res_Genotype_DNMT3B_vs_WT_noshrink),as.data.frame(res_Genotype_DNMT3B_vs_WT_ordered_sval),by=0)

rownames(res_Genotype_DNMT3B_vs_WT_ordered_merged)<-res_Genotype_DNMT3B_vs_WT_ordered_merged$Row.names
res_Genotype_DNMT3B_vs_WT_ordered_merged<-res_Genotype_DNMT3B_vs_WT_ordered_merged[,-1]

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_annotation <-getBM(attributes = c("ensembl_gene_id","gene_biotype","external_gene_name","chromosome_name","start_position","end_position"),filters = "ensembl_gene_id",values = rownames(res_Genotype_DNMT3B_vs_WT_ordered_merged),mart = ensembl)
rownames(gene_annotation)<-gene_annotation$ensembl_gene_id
res_Genotype_DNMT3B_vs_WT_ordered_merged_annot<-merge(as.data.frame(res_Genotype_DNMT3B_vs_WT_ordered_merged),gene_annotation,by=0)
rownames(res_Genotype_DNMT3B_vs_WT_ordered_merged_annot)<-res_Genotype_DNMT3B_vs_WT_ordered_merged_annot$Row.names
res_Genotype_DNMT3B_vs_WT_ordered_merged_annot<-res_Genotype_DNMT3B_vs_WT_ordered_merged_annot[order(res_Genotype_DNMT3B_vs_WT_ordered_merged_annot$svalue),]
#write.table(as_tibble(res_Genotype_DNMT3B_vs_WT_ordered_merged_annot[,-c(1,4,8,10,12)],rownames = "Ensembl_id"),sep="\t",row.names = F,file="Deseq2_Genotype_DNMT3B_vs_WT_sval_annot_190623.txt")



ncounts_ioannis<-counts(dds_ioannis, normalized=TRUE)
ntd_ioannis <- normTransform(dds_ioannis)

#write.table(ncounts,sep="\t",file="normalized_counts.txt")

#annotmarkers<-read.table("topgenes.txt",sep="\t",header = FALSE)
#annotmarkers<-as.vector(annotmarkers[,1])
#annotmarkers2<-read.table("top1000_EB14_Nodox_vs_Ch9_noDox_DOWN.txt",sep="\t",header = FALSE)
#annotmarkers2<-as.vector(annotmarkers2[,1])

#selectpluri<-assay(ntd[rownames(ntd) %in% annotmarkers,])
#selectpluri2<-cts[rownames(cts) %in% annotmarkers2,]

#pheatmap(selectpluri, cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, col= viridis(20),scale="row")



plotCounts(dds_ioannis, gene="ENSG00000169855", intgroup=c("Genotype"),main = "ROBO1")


vst <- assay(vst(dds))
p <- pca(vst, metadata = coldata, removeVar = 0.1)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p, showLoadings = TRUE,labSize = 3, pointSize = 3, sizeLoadingsNames = 3)
plotloadings(p, labSize = 2)



############## Go enrichment #############
########### All samples #######

topGenes2<-subset(res_Genotype_DNMT3B_vs_WT, padj <0.05)
topGenes<-row.names(subset(topGenes2, abs(log2FoldChange)> 2))



egoBP <- enrichGO(gene         = topGenes,
                  OrgDb         = org.Hs.eg.db,
                  universe = row.names(res_Genotype_DNMT3B_vs_WT),
                  keyType = "ENSEMBL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 10, maxGSSize = 5000,
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable = FALSE)
head(egoBP)
dotegoup<-dotplot(egoBP)+ ggtitle("topgenes_padj<0.05_abs(log2F)>2 - BP") +
  guides(color = guide_colorbar(order=2),size = guide_legend(order=1)) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.001,
                                   decimal.mark = '.'))

dotegoup 

egoMF <- enrichGO(gene         = topGenes,
                  OrgDb         = org.Hs.eg.db,
                  universe = row.names(res_Genotype_DNMT3B_vs_WT),
                  keyType = "ENSEMBL",
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 10, maxGSSize = 5000,
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable = FALSE)
head(egoMF)
dotegoup2<-dotplot(egoMF)+ ggtitle("topgenes_padj<0.05_abs(log2F)>2 - MF") +
  guides(color = guide_colorbar(order=2),size = guide_legend(order=1)) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.001,
                                   decimal.mark = '.'))

dotegoup2

egoCC <- enrichGO(gene         = topGenes,
                  OrgDb         = org.Hs.eg.db,
                  universe = row.names(res_Genotype_DNMT3B_vs_WT),
                  keyType = "ENSEMBL",
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 10, maxGSSize = 5000,
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable = FALSE)

head(egoCC)
dotegoup3<-dotplot(egoCC)+ ggtitle("topgenes_padj<0.05_abs(log2F)>2 - CC") +
  guides(color = guide_colorbar(order=2),size = guide_legend(order=1)) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.001,
                                   decimal.mark = '.'))

dotegoup3
topGenes3<-row.names(subset(topGenes2, abs(log2FoldChange)> 2))

ids = bitr(topGenes3,fromType="ENSEMBL",toType="ENTREZID", OrgDb="org.Hs.eg.db")

kk <- enrichKEGG(gene     = ids$ENTREZID,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05)
head(kk)
dotegoupkk<-dotplot(kk)+ ggtitle("cond_B_vs_A abs(log2)>4 - KEGG") +
  guides(color = guide_colorbar(order=2),size = guide_legend(order=1)) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.001,
                                   decimal.mark = '.'))

dotegoupkk

kegg_mmu04151_entrez<-c("18053","12444","16773","20393","30955","320207","11600","16774","16404","263803","12841","54635","18595","14254","12826","12827","16772","81877","14164","12985","18707","19699","18212","12829","12443","12828")
kegg_mmu04151_ensembl = bitr(kegg_mmu04151_entrez,fromType="ENTREZID",toType=c("SYMBOL"), OrgDb="org.Mm.eg.db")

############## Heatmaps #############

selectmarkersUP<-read.table("top50UP.txt",sep="\t",header = TRUE)
selectpluriUP<-ncounts[rownames(ncounts) %in% selectmarkersUP$ens,]
pheatmap(selectpluriUP, cluster_rows=FALSE, labels_row=selectmarkersUP$name,show_rownames=TRUE,cluster_cols=FALSE, col= viridis(20),scale="row",main = "Top 50 Upregulated Genes (normalised counts)")
pheatmap(assay(ntd)[selectmarkersUP$ens,], labels_row=selectmarkersUP$name,cluster_rows=FALSE, cluster_cols=FALSE,show_rownames=TRUE,col= viridis(20),main = "Top 50 Upregulated Genes (normTransform)")

selectmarkersDOWN<-read.table("top50DOWN.txt",sep="\t",header = TRUE)
selectpluriDOWN<-ncounts[rownames(ncounts) %in% selectmarkersDOWN$ens,]
pheatmap(selectpluriUP, cluster_rows=FALSE, labels_row=selectmarkersUP$name,show_rownames=TRUE,cluster_cols=FALSE, col= viridis(20),scale="row",main = "Top 50 Downregulated Genes (normalised counts)")
pheatmap(assay(ntd)[selectmarkersDOWN$ens,], labels_row=selectmarkersDOWN$name,cluster_rows=FALSE, cluster_cols=FALSE,show_rownames=TRUE,col= viridis(20),main = "Top 50 Downregulated Genes (normTransform)")

############## Volcano Plot #############

annot_table<-read.table("Deseq2_Genotype_DNMT3B_vs_WT_sval_annot_190623.txt",sep="\t",header = TRUE)
EnhancedVolcano(annot_table,
                lab = annot_table$external_gene_name,
                x = 'log2FoldChange.y',
                y = 'padj',
                pCutoff = 0.05,
                title = 'DNMT3B vs WT',
                FCcutoff = 2,
                xlim = c(-10, 10),
                labSize = 4.0,
                drawConnectors = TRUE)

