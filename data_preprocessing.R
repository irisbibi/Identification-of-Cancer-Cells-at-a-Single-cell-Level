setwd('~/python/02718_computational_medicine')

library(data.table)
library(SingleCellExperiment)
library(scater)
library(edgeR)

library(dplyr)
library(Seurat)
library(patchwork)

options(stringsAsFactors = FALSE)

############ GSE144735 ############
expr_matrix <- fread("GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt")
anno <- read.table(gzfile("GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz"), sep = '\t')

num_genes = dim(expr_matrix)[1]
gene_names = c(expr_matrix[1:num_genes, 1]$Index)
expr_matrix = expr_matrix[,-1]
row.names(expr_matrix) = gene_names

colnames(anno) = anno[1,]
anno = anno[-1,]
anno = anno[which(anno$Class != "Border"),]
index_rmBorder = c(anno$Index)
expr_matrix = expr_matrix[ , ..index_rmBorder]
row.names(expr_matrix) = gene_names

num_normal = dim(anno[which(anno$Class == "Normal"),])[1]
num_tumor = dim(anno[which(anno$Class == "Tumor"),])[1]

labels = anno$Class
labels = as.factor(labels)

seurat_obj <- CreateSeuratObject(counts = expr_matrix, project = "GSE144735", min.cells = 3, min.features = 200)
seurat_obj@meta.data$orig.ident <- labels
seurat_obj@meta.data$active.ident <- labels
seurat_obj@active.ident <- labels

##### QC
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "active.ident")

plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# filter cells that have unique feature counts over 6000 or less than 200
# filter cells that have >20% mitochondrial counts
seurat_obj_test <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)


##### Normalize
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)


##### Feature selection
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
# plot1 + plot2


##### Scale
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)


##### PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
# Examine and visualize PCA results a few different ways
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
# visualize
DimPlot(seurat_obj, reduction = "pca")

DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)

# determine the dimensionality
seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)

JackStrawPlot(seurat_obj, dims = 1:15)

ElbowPlot(seurat_obj)


##### Non-linear dimensional reduction (UMAP/tSNE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

DimPlot(seurat_obj, reduction = "umap")


##### Cluster biomarkers
Idents(object = seurat_obj) <- "orig.ident"
tumor.markers <- FindMarkers(seurat_obj, ident.1 = "Tumor", ident.2 = "Normal", min.pct = 0.25, method = "DESeq2")
head(tumor.markers, n = 5)

normal.markers <- FindMarkers(seurat_obj, ident.1 = "Normal", ident.2 = "Tumor", min.pct = 0.25, method = "DESeq2")
head(normal.markers, n = 5)

# DEG
DEGs = rownames(tumor.markers)

fileConn <- file("DEGs.txt")
writeLines(DEGs, fileConn)
close(fileConn)

#plot volcano
library(ggplot2)
DEG=as.data.frame(tumor.markers)
plot(DEG$avg_log2FC,-log10(DEG$p_val_adj))
logFC_cutoff=1
DEG$change = as.factor(ifelse(DEG$p_val_adj < 0.05 & abs(DEG$avg_log2FC) > logFC_cutoff,
                              ifelse(DEG$avg_log2FC > logFC_cutoff ,'UP','DOWN'),'NOT'))
table(DEG$change)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',]))
g = ggplot(data=DEG, aes(,x=avg_log2FC, y=-log10(p_val_adj),color=change)) + geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+ xlab("log2 fold change") + ylab("-log10 p-value") + 
  ggtitle( this_tile ) +   theme(plot.title = element_text(size=15,hjust = 0.5)) + 
  scale_colour_manual(values = c('blue','grey','red')) 	## corresponding to the levels(res$change)
print(g)

# find markers for every cluster compared to all remaining cells, report only the positive ones
all.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

DEG=as.data.frame(tumor.markers)
diff <- DEG
logFC <-diff$avg_log2FC
adj <- diff$p_val_adj
data <- data.frame(logFC=logFC,padj=adj)
data$sig[(data$padj > 0.05|data$padj=="NA")|(data$logFC < 0.5)& data$logFC > -0.5] <- "no"
data$sig[data$padj <= 0.05 & data$logFC >= 0.5] <- "up"
data$sig[data$padj <= 0.05 & data$logFC <= -0.5] <- "down"


x_lim <- max(logFC,-logFC)
library(ggplot2)
library(RColorBrewer)
pdf(file = "DEG_volcano.pdf",width=8,height=8)
theme_set(theme_bw())
p <- ggplot(data,aes(logFC,-1*log10(padj),
                     color = sig))+geom_point()+
  xlim(-5,5) +  labs(x="log2(FoldChange)",y="-log10(FDR)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-0.5,0.5),linetype=4)
# p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))+ylim(0,15)
p <- p  +guides(colour = FALSE)
# p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
p
dev.off()
print(p)


VlnPlot(seurat_obj, features = c("JCHAIN", "CFD", "IGHA1", "IGHA2", "IGKC"))
VlnPlot(seurat_obj, features = c("COL1A1", "COL3A1", "IGHG3", "BGN", "KRT19"))
VlnPlot(seurat_obj, features = DEGs[1:5])


FeaturePlot(seurat_obj, features = c("JCHAIN", "COL1A1", "CFD", "COL3A1", "IGHA1", "IGHG3"))
FeaturePlot(seurat_obj, features = DEGs[1:6])


##### write sparse matrix to csv
norm_data = seurat_obj[["RNA"]]@data
norm_data = norm_data[VariableFeatures(seurat_obj),]
write_sparse_csv <- function(x, file, ..., chunk = 100){
  passes <- nrow(x) %/% chunk
  remaining <- nrow(x) %% chunk
  if(passes > 0){
    inx <- seq_len(chunk)
    y <- x[inx, , drop = FALSE]
    y <- as.matrix(y)
    write.table(y, file, append = FALSE, sep = ",", col.names = !is.null(colnames(x)), ...)
    passes <- passes - 1L
    for(i in seq_len(passes)){
      inx <- inx + chunk
      y <- x[inx, , drop = FALSE]
      y <- as.matrix(y)
      write.table(y, file, append = TRUE, sep = ",", col.names = FALSE,  ...)
    }
    if(remaining > 0){
      inx <- inx + remaining
      y <- x[inx, , drop = FALSE]
      y <- as.matrix(y)
      write.table(y, file, append = TRUE, sep = ",", col.names = FALSE, ...)
    }
  } else if(remaining > 0){
    inx <- seq_len(remaining)
    y <- x[inx, , drop = FALSE]
    y <- as.matrix(y)
    write.table(y, file, append = FALSE, sep = ",", col.names = FALSE, ...)
  }
}
write_sparse_csv(norm_data, "norm_data_rmBorder_allFeat.csv")

write.csv(anno, "annotation_rmBorder.csv")




############ GSE200997 ############
expr_matrix <- fread("GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv")
anno <- read.table(gzfile("GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz"), sep = ',')

num_genes = dim(expr_matrix)[1]
gene_names = c(expr_matrix[1:num_genes, 1]$V1)
expr_matrix = expr_matrix[,-1]
row.names(expr_matrix) = gene_names

colnames(anno) = anno[1,]
anno = anno[-1,]

num_normal = dim(anno[which(anno$Condition == "Normal"),])[1]
num_tumor = dim(anno[which(anno$Condition == "Tumor"),])[1]
index_normal = c(anno[which(anno$Condition == "Normal"),][1:10000,1])
index_tumor = c(anno[which(anno$Condition == "Tumor"),][1:10000,1])
anno = rbind(anno[which(anno$Condition == "Normal"),][1:10000,], anno[which(anno$Condition == "Tumor"),][1:10000,])
# inex_keep = c(index_normal, index_tumor)
expr_matrix = data.table(expr_matrix[ , ..index_normal], expr_matrix[ , ..index_tumor])
row.names(expr_matrix) = gene_names

labels = c(rep("Normal", 10000), rep("Tumor", 10000))
labels = as.factor(labels)


library(dplyr)
library(Seurat)
library(patchwork)

seurat_obj <- CreateSeuratObject(counts = expr_matrix, project = "GSE200997", min.cells = 3, min.features = 200)
seurat_obj@meta.data$orig.ident <- labels
seurat_obj@meta.data$active.ident <- labels
seurat_obj@active.ident <- labels


##### QC
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "active.ident")

plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# filter cells that have unique feature counts over 6000 or less than 200
# filter cells that have >20% mitochondrial counts
# seurat_obj_test <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)


##### Normalize
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)


##### Feature selection
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
# plot1 + plot2


##### Scale
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

##### PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
# Examine and visualize PCA results a few different ways
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
# visualize
DimPlot(seurat_obj, reduction = "pca")

DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)

# determine the dimensionality
seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)

JackStrawPlot(seurat_obj, dims = 1:15)

ElbowPlot(seurat_obj)


##### Non-linear dimensional reduction (UMAP/tSNE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap")

# save normalized data
norm_data = seurat_obj[["RNA"]]@data
write_sparse_csv(norm_data, "GSE200997_norm_data.csv")

write.csv(anno, "GSE200997_annotation.csv")






############ GSE132465 ############

expr_matrix <- fread("GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt")
anno <- read.table(gzfile("GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"), sep = '\t')

num_genes = dim(expr_matrix)[1]
gene_names = c(expr_matrix[1:num_genes, 1]$Index)
expr_matrix = expr_matrix[,-1]
row.names(expr_matrix) = gene_names

colnames(anno) = anno[1,]
anno = anno[-1,]
# anno = anno[which(anno$Class != "Border"),]
# index_rmBorder = c(anno$Index)
# expr_matrix = expr_matrix[ , ..index_rmBorder]
# row.names(expr_matrix) = gene_names

num_normal = dim(anno[which(anno$Class == "Normal"),])[1]
num_tumor = dim(anno[which(anno$Class == "Tumor"),])[1]
index_normal = c(anno[which(anno$Class == "Normal"),][1:10000,1])
index_tumor = c(anno[which(anno$Class == "Tumor"),][1:10000,1])
anno = rbind(anno[which(anno$Class == "Normal"),][1:10000,], anno[which(anno$Class == "Tumor"),][1:10000,])
# inex_keep = c(index_normal, index_tumor)
expr_matrix = data.table(expr_matrix[ , ..index_normal], expr_matrix[ , ..index_tumor])
row.names(expr_matrix) = gene_names

labels = anno$Class
labels = as.factor(labels)

# umi <- SingleCellExperiment(
#   assays = list(counts = as.matrix(molecules)), 
#   colData = anno
# )

# keep_feature <- rowSums(counts(umi)) > 0
# umi <- umi[keep_feature, ]
# 
# sizeFactors(assay(umi)) <- edgeR::calcNormFactors(counts(assay(umi)), method = "TMM")


library(dplyr)
library(Seurat)
library(patchwork)

seurat_obj <- CreateSeuratObject(counts = expr_matrix, project = "GSE132465", min.cells = 3, min.features = 200)
seurat_obj@meta.data$orig.ident <- labels
seurat_obj@meta.data$active.ident <- labels
seurat_obj@active.ident <- labels


##### QC
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "active.ident")

plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# filter cells that have unique feature counts over 6000 or less than 200
# filter cells that have >20% mitochondrial counts
# seurat_obj_test <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)


##### Normalize
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

##### Feature selection
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
# plot1 + plot2


##### Scale
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

##### PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
# Examine and visualize PCA results a few different ways
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
# visualize
DimPlot(seurat_obj, reduction = "pca")

DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)

# determine the dimensionality
seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)

JackStrawPlot(seurat_obj, dims = 1:15)

ElbowPlot(seurat_obj)


##### Non-linear dimensional reduction (UMAP/tSNE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

DimPlot(seurat_obj, reduction = "umap")


##### Cluster biomarkers
Idents(object = seurat_obj) <- "orig.ident"
tumor.markers <- FindMarkers(seurat_obj, ident.1 = "Tumor", ident.2 = "Normal", min.pct = 0.25, method = "DESeq2")
head(tumor.markers, n = 5)

normal.markers <- FindMarkers(seurat_obj, ident.1 = "Normal", ident.2 = "Tumor", min.pct = 0.25, method = "DESeq2")
head(normal.markers, n = 5)

# DEG
DEG = as.data.frame(tumor.markers)
DEGs =c(rownames(DEG[which(DEG$avg_log2FC > 1),]), rownames(DEG[which(DEG$avg_log2FC < -1),]))
# DEGs = rownames(tumor.markers)

fileConn <- file("DEGs.txt")
writeLines(DEGs, fileConn)
close(fileConn)

#plot volcano
library(ggplot2)
DEG=as.data.frame(tumor.markers)
plot(DEG$avg_log2FC,-log10(DEG$p_val_adj))
logFC_cutoff=1
DEG$change = as.factor(ifelse(DEG$p_val_adj < 0.05 & abs(DEG$avg_log2FC) > logFC_cutoff,
                              ifelse(DEG$avg_log2FC > logFC_cutoff ,'UP','DOWN'),'NOT'))
table(DEG$change)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',]))
g = ggplot(data=DEG, aes(,x=avg_log2FC, y=-log10(p_val_adj),color=change)) + geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+ xlab("log2 fold change") + ylab("-log10 p-value") + 
  ggtitle( this_tile ) +   theme(plot.title = element_text(size=15,hjust = 0.5)) + 
  scale_colour_manual(values = c('blue','black','red')) 	## corresponding to the levels(res$change)
print(g)

# find markers for every cluster compared to all remaining cells, report only the positive ones
all.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

VlnPlot(seurat_obj, features = c("SRGN", "VIM", "CCL5", "IGFBP7", "CFD"))
VlnPlot(seurat_obj, features = c("LCN2", "KRT18", "S100P", "TFF3", "AGR2"))
VlnPlot(seurat_obj, features = DEGs[1:5])


FeaturePlot(seurat_obj, features = c("SRGN", "LCN2", "VIM", "KRT18", "CCL5", "S100P"))
FeaturePlot(seurat_obj, features = DEGs[1:6])


##### write sparse matrix to csv
norm_data = seurat_obj[["RNA"]]@data
# norm_data = norm_data[VariableFeatures(seurat_obj),]
write_sparse_csv <- function(x, file, ..., chunk = 100){
  passes <- nrow(x) %/% chunk
  remaining <- nrow(x) %% chunk
  if(passes > 0){
    inx <- seq_len(chunk)
    y <- x[inx, , drop = FALSE]
    y <- as.matrix(y)
    write.table(y, file, append = FALSE, sep = ",", col.names = !is.null(colnames(x)), ...)
    passes <- passes - 1L
    for(i in seq_len(passes)){
      inx <- inx + chunk
      y <- x[inx, , drop = FALSE]
      y <- as.matrix(y)
      write.table(y, file, append = TRUE, sep = ",", col.names = FALSE,  ...)
    }
    if(remaining > 0){
      inx <- inx + remaining
      y <- x[inx, , drop = FALSE]
      y <- as.matrix(y)
      write.table(y, file, append = TRUE, sep = ",", col.names = FALSE, ...)
    }
  } else if(remaining > 0){
    inx <- seq_len(remaining)
    y <- x[inx, , drop = FALSE]
    y <- as.matrix(y)
    write.table(y, file, append = FALSE, sep = ",", col.names = FALSE, ...)
  }
}
write_sparse_csv(norm_data, "GSE132465_norm_data_allFeat.csv")

write.csv(anno, "GSE132465_annotation.csv")
