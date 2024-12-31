library(Seurat)
library(ggplot2)
library(reshape2)

#Merged SRAT ----
CFG <- list()
CFG$QUAL <- 30 #Higher than 30
CFG$ndims <- 25
CFG$ndims.p0 <- 15
CFG$ndims.p1 <- 15
CFG$ndims.tumor <- 15

srat.lr <- readRDS(file = '~/nanopore_souporcell/srat-aneuploidy.rds')
srat.lr@meta.data <- srat.lr@meta.data %>%
  mutate(reference = case_when(
    patient == '0' & is_aneuploid == 'yes' ~ 'p0.tumor',
    patient == '0' & is_aneuploid == 'no' ~ 'p0.healthy',
    patient == '1' & is_aneuploid == 'yes' ~ 'p1.tumor',
    patient == '1' & is_aneuploid == 'no' ~ 'p1.healthy',
    TRUE ~ 'NA'
  ))

srat.meta <- data.frame(srat.lr@meta.data)

#vafs.merged<-read.table(file = '~/snv_calling/processing.vafs/vafs.merged.tsv', header = TRUE)
alt <- read.table(file = '~/sop_input/alt.mtx', header = TRUE)
alt$ID <- NULL
alt<- alt %>% unite(col='SNV',CHROM:ALT, sep=':')
rownames(alt) <- alt$SNV
rownames(alt) <- paste0('alt:', rownames(alt))
alt$SNV <- NULL
alt$QUAL <- alt[alt$QUAL > CFG$QUAL]
alt$QUAL <- NULL
alt <- as.matrix(alt)

ref <- read.table(file = '~/sop_input/ref.mtx', header = TRUE)
ref$ID <- NULL
ref<- ref %>% unite(col='SNV',CHROM:ALT, sep=':')
rownames(ref) <- ref$SNV
rownames(ref) <- paste0('ref:', rownames(ref))
ref$SNV <- NULL
ref <- ref[ref$QUAL > CFG$QUAL]
ref$QUAL <- NULL
ref <- as.matrix(ref)

counts.per.cell <- (rbind(alt, ref))
counts.per.cell <- as.matrix(counts.per.cell)

mitochondrial.snvs <- grepl("^alt:chrM|^ref:chrM", rownames(counts.per.cell))

mito.matrix <- counts.per.cell[mitochondrial.snvs, ]
nuc.matrix <- counts.per.cell[!mitochondrial.snvs, ]

nuc.matrix <- nuc.matrix[, colnames(nuc.matrix) %in% gsub('LX279_','',colnames(srat.lr))]
dim(nuc.matrix)


srat.snvs <- Seurat::CreateSeuratObject(counts = nuc.matrix, project = "ONT_snvs", meta = srat.meta)
colnames(srat.snvs) <- paste0('LX279_', colnames(srat.snvs))
srat.snvs <- NormalizeData(srat.snvs, normalization.method = "LogNormalize")

srat.snvs <- ScaleData(srat.snvs, features = rownames(srat.snvs), verbose = FALSE)

srat.snvs <- FindVariableFeatures(srat.snvs)


srat.snvs <- suppressMessages(SCTransform(srat.snvs, 
                                     vars.to.regres=NULL,
                                     vst.flavor='v2', ## version that solved some problems
                                     verbose=FALSE,
                                     variable.features.n=3000))

top10 <- head(VariableFeatures(srat.snvs), 10)

DefaultAssay(srat.snvs)
srat.snvs <- RunPCA(srat.snvs, npcs = 50)

elbow <- ElbowPlot(srat.snvs, reduction = "pca", ndims = 50)
elbow #geom_vline(xintercept = c(seq(from = 0, to = 50, by = 1)))

DimHeatmap(srat.snvs, dims = (1:12), cells = 100)

srat.snvs <- RunUMAP(srat.snvs, dims = 1:CFG$ndims, assay = 'RNA')

DimPlot(srat.snvs, group.by = 'reference')

srat.snvs <- FindNeighbors(srat.snvs, dims = 1:CFG$ndims)
## Computing nearest neighbor graph
## Computing SNN

srat.snvs <- FindClusters(srat.snvs, resolution = seq(0.4, 2, 0.2), algorithm = 1)

srat.snvs <- FindClusters(srat.snvs, resolution = 0.8, algorithm = 1)

columns <- grep("SCT_snn_res", names(srat.snvs@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the
## cluster ids)
for (col in columns) {
  cat("\n====\n", col, ":")
  print(table(srat.snvs@meta.data[[col]]))
}

clus.column <- paste0("SCT_snn_res.", 1.8)
Idents(srat) <- clus.column

clusters <- sort(unique(srat.snvs@meta.data[[clus.column]]))
## switching color scheme to make it look different:
cluster.colors <- Polychrome::alphabet.colors(length(clusters))
names(cluster.colors) <- as.character(clusters)

## and show result with the other umaps to get your
## bearings:
p_clus <- DimPlot(srat.snvs, reduction = "umap", pt.size = 0.5, group.by = clus.column,
                  cols = cluster.colors, label=TRUE)

p_clus
rm(p_clusbyreso)
gc()

#Error in check.length(gparname) : 
#  'gpar' element 'lwd' must not be length 0
#In addition: Warning message:
#  The `add` argument of `group_by()` is deprecated as of dplyr 1.0.0.
#ℹ Please use the `.add` argument instead.
#ℹ The deprecated feature was likely used in the dplyr package.
#Please report the issue at <https://github.com/tidyverse/dplyr/issues>.
#This warning is displayed once every 8 hours.
#Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated. 

## based on this I think it's still good to choose 0.8
clus.column <- paste0("SCT_snn_res.", 1.4)
Idents(srat.snvs) <- clus.column

typenames <- unique(srat.snvs@meta.data$singler)
singler.colors <- Polychrome::palette36.colors(length(typenames))
names(singler.colors) <- typenames

DimPlot(srat.snvs, group.by = 'reference') | FeaturePlot(srat.snvs, features = 'aneuploidy_score')

de_snvs <- FindAllMarkers(srat.snvs, assay = "RNA",
                           min.pct = 0.1, logfc.threshold = log(1.5))
## The resulting data.frame is too unwieldy to work with,
## the top 10 genes per cluster are enough:

de_snvs %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
head(top10)

DoHeatmap(subset(srat.snvs, downsample = 50), features = top10$gene,
          group.colors = cluster.colors, assay = "RNA", slot = "scale.data") +
  theme(axis.text.y = element_text(size = 6))

#saveRDS(srat.snvs, file = '~/srat.snvs-clusters.rds')
#srat.snvs <- readRDS(file = '~/srat.snvs-clusters.rds')

srat.tumor.p0 <- subset(srat.snvs, subset = reference == 'p0.tumor')
srat.tumor.p1 <- subset(srat.snvs, subset = reference == 'p1.tumor')

srat.ge.tumor.p0 <- subset(srat.lr, subset = reference == 'p0.tumor')
srat.ge.tumor.p1 <- subset(srat.lr, subset = reference == 'p1.tumor')

srat.lr <- AddMetaData(srat.lr, metadata = srat.snvs$SCT_snn_res.1.8, col.name = 'snvs.clustering.1.8')
DimPlot(srat.lr, group.by = 'snvs.clustering.1.8', label = TRUE, cols = rainbow(length(unique(srat.lr$snvs.clustering.1.8))))

#Tumor 0 ----
srat.tumor.p0 <- NormalizeData(srat.tumor.p0, normalization.method = "LogNormalize")

srat.tumor.p0 <- ScaleData(srat.tumor.p0, features = rownames(srat.tumor.p0), verbose = FALSE)

srat.tumor.p0 <- FindVariableFeatures(srat.tumor.p0)


srat.tumor.p0 <- suppressMessages(SCTransform(srat.tumor.p0, 
                                          vars.to.regres=NULL,
                                          vst.flavor='v2', ## version that solved some problems
                                          verbose=FALSE,
                                          variable.features.n=3000))

top10 <- head(VariableFeatures(srat.tumor.p0), 10)

DefaultAssay(srat.tumor.p0)
srat.tumor.p0 <- RunPCA(srat.tumor.p0, npcs = 50)

elbow <- ElbowPlot(srat.tumor.p0, reduction = "pca", ndims = 50)
elbow #geom_vline(xintercept = c(seq(from = 0, to = 50, by = 1)))

DimHeatmap(srat.tumor.p0, dims = (6:18), cells = 100)

srat.tumor.p0 <- RunUMAP(srat.tumor.p0, dims = 1:CFG$ndims.p0, assay = 'RNA')
srat.tumor.p0 <- FindNeighbors(srat.tumor.p0, dims = 1:CFG$ndims.p0)
## Computing nearest neighbor graph
## Computing SNN

## Using clustree ----
srat.tumor.p0 <- FindClusters(srat.tumor.p0, resolution = seq(0.4, 2, 0.2), algorithm = 1)

columns <- grep("SCT_snn_res", names(srat.snvs@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the
## cluster ids)
for (col in columns) {
  cat("\n====\n", col, ":")
  print(table(srat.tumor.p0@meta.data[[col]]))
}

clus.column <- paste0("SCT_snn_res.", 1)
Idents(srat.tumor.p0) <- clus.column

clusters <- sort(unique(srat.tumor.p0@meta.data[[clus.column]]))
## switching color scheme to make it look different:
cluster.colors <- Polychrome::alphabet.colors(length(clusters))
names(cluster.colors) <- as.character(clusters)

## and show result with the other umaps to get your
## bearings:
p_clus <- DimPlot(srat.tumor.p0, reduction = "umap", pt.size = 0.5, group.by = clus.column,
                  cols = cluster.colors, label=TRUE)
p_clus | DimPlot(srat.tumor.p0, group.by = 'singler')

p_clus
gc()

#Error in check.length(gparname) : 
#  'gpar' element 'lwd' must not be length 0
#In addition: Warning message:
#  The `add` argument of `group_by()` is deprecated as of dplyr 1.0.0.
#ℹ Please use the `.add` argument instead.
#ℹ The deprecated feature was likely used in the dplyr package.
#Please report the issue at <https://github.com/tidyverse/dplyr/issues>.
#This warning is displayed once every 8 hours.
#Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated. 

## based on this I think it's still good to choose 0.8
clus.column <- paste0("SCT_snn_res.", 1)
Idents(srat.tumor.p0) <- clus.column

typenames <- unique(srat.tumor.p0@meta.data$singler)
singler.colors <- Polychrome::palette36.colors(length(typenames))
names(singler.colors) <- typenames

p_clus | FeaturePlot(srat.tumor.p0, features = 'aneuploidy_score')

de_snvs <- FindAllMarkers(srat.tumor.p0, assay = "RNA",
                           min.pct = 0.1, logfc.threshold = log(1.5))
## The resulting data.frame is too unwieldy to work with,
## the top 10 genes per cluster are enough:

de_snvs %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
head(top10)

DoHeatmap(subset(srat.tumor.p0, downsample = 50), features = top10$gene,
          group.colors = cluster.colors, assay = "RNA", slot = "scale.data") +
  theme(axis.text.y = element_text(size = 6))
VlnPlot(srat.tumor.p0, features = c('nCount_RNA', 'aneuploidy_score'))

p_clus | FeaturePlot(srat.tumor.p0, features = 'nCount_RNA')
p_clus | FeaturePlot(srat.tumor.p0, features = 'aneuploidy_score')
p_clus | DimPlot(srat.tumor.p0, group.by = 'singler')

Idents(srat.ge.tumor.p0) <- srat.tumor.p0$SCT_snn_res.1
de_genes_snvs_p0 <- FindAllMarkers(srat.ge.tumor.p0, assay = "RNA",
               min.pct = 0.1, logfc.threshold = log(1.5))

de_genes_snvs_p0 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) -> top10.p0.ge
head(top10.p0.ge)

de_genes_snvs_p0 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10.p0.ge.top_n
head(top10.p0.ge.top_n)

DoHeatmap(subset(srat.ge.tumor.p0, downsample = 50), features = top10.p0.ge.top_n$gene,
          group.colors = cluster.colors, assay = "RNA", slot = "scale.data") +
  theme(axis.text.y = element_text(size = 6))

VlnPlot(srat.ge.tumor.p0, features = c('nCount_RNA', 'aneuploidy_score'))

#Tumor 1 ----
srat.tumor.p1 <- NormalizeData(srat.tumor.p1, normalization.method = "LogNormalize")

srat.tumor.p1 <- ScaleData(srat.tumor.p1, features = rownames(srat.tumor.p1), verbose = FALSE)

srat.tumor.p1 <- FindVariableFeatures(srat.tumor.p1)


srat.tumor.p1 <- suppressMessages(SCTransform(srat.tumor.p1, 
                                              vars.to.regres=NULL,
                                              vst.flavor='v2', ## version that solved some problems
                                              verbose=FALSE,
                                              variable.features.n=3000))

top10 <- head(VariableFeatures(srat.tumor.p1), 10)

DefaultAssay(srat.tumor.p1)
srat.tumor.p1 <- RunPCA(srat.tumor.p1, npcs = 50)

elbow <- ElbowPlot(srat.tumor.p1, reduction = "pca", ndims = 50)
elbow #geom_vline(xintercept = c(seq(from = 0, to = 50, by = 1)))

DimHeatmap(srat.tumor.p1, dims = (6:18), cells = 100)

srat.tumor.p1 <- RunUMAP(srat.tumor.p1, dims = 1:CFG$ndims.p1, assay = 'RNA')

srat.tumor.p1 <- FindNeighbors(srat.tumor.p1, dims = 1:CFG$ndims.p1)
## Computing nearest neighbor graph
## Computing SNN

srat.tumor.p1 <- FindClusters(srat.tumor.p1, resolution = seq(0.4, 2, 0.2), algorithm = 1)

columns <- grep("SCT_snn_res", names(srat.snvs@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the
## cluster ids)
for (col in columns) {
  cat("\n====\n", col, ":")
  print(table(srat.tumor.p1@meta.data[[col]]))
}

clus.column <- paste0("SCT_snn_res.", 1)
Idents(srat.tumor.p1) <- clus.column

clusters <- sort(unique(srat.tumor.p1@meta.data[[clus.column]]))
## switching color scheme to make it look different:
cluster.colors.p1.snvs <- Polychrome::alphabet.colors(length(clusters))
names(cluster.colors.p1.snvs) <- as.character(clusters)

## and show result with the other umaps to get your
## bearings:
p_clus_snvs_p1 <- DimPlot(srat.tumor.p1, reduction = "umap", pt.size = 0.5, group.by = clus.column,
                  cols = cluster.colors.p1.snvs, label=TRUE)
##p_clus_snvs_p1 ----
p_clus_snvs_p1
gc()

## based on this I think it's still good to choose 1
clus.column <- paste0("SCT_snn_res.", 1)
Idents(srat.tumor.p1) <- clus.column

typenames <- unique(srat.tumor.p1@meta.data$singler)
singler.colors <- Polychrome::palette36.colors(length(typenames))
names(singler.colors) <- typenames

p_clus_snvs_p1 | FeaturePlot(srat.tumor.p1, features = 'aneuploidy_score')

##DE snvs ----
de_snvs_p1 <- FindAllMarkers(srat.tumor.p1, assay = "RNA",
                           min.pct = 0.1, logfc.threshold = log(1.5))
## The resulting data.frame is too unwieldy to work with,
## the top 10 genes per cluster are enough:

de_snvs_p1 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10.p1.snvs
head(top10.p1.snvs)

#Looking at chr17
probably_BIRC5 <- de_snvs_p1[grepl("^alt:chr17:78", de_snvs_p1$gene), ]
head(probably_BIRC5)

DoHeatmap(subset(srat.tumor.p1, downsample = 50), features = top10.p1.snvs$gene,
          group.colors = cluster.colors, assay = "RNA", slot = "scale.data") +
  theme(axis.text.y = element_text(size = 6))

VlnPlot(srat.tumor.p1, features = c('nCount_RNA', 'aneuploidy_score'))

p_clus_snvs_p1 | FeaturePlot(srat.tumor.p1, features = 'nCount_RNA')
p_clus_snvs_p1 | FeaturePlot(srat.tumor.p1, features = 'aneuploidy_score')
p_clus_snvs_p1 | DimPlot(srat.tumor.p1, group.by = 'singler')

Idents(srat.ge.tumor.p1) <- srat.tumor.p1$SCT_snn_res.1
de_genes_snvs_p1 <- FindAllMarkers(srat.ge.tumor.p1, assay = "RNA",
                                   min.pct = 0.1, logfc.threshold = log(1.5))

de_genes_snvs_p1 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) -> top10.p1.ge
head(top10.p1.ge)

de_genes_snvs_p1 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10.p1.ge.top_n
head(top10.p1.ge.top_n)

Idents(srat.ge.tumor.p1) <- srat.ge.tumor.p1$SCT_snn_res.1
de_genes_ge <- FindAllMarkers(srat.ge.tumor.p1, assay = "RNA",
                                   min.pct = 0.1, logfc.threshold = log(1.5))
de_genes_ge %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) -> top10.p1.ge.only
head(top10.p1.ge.only)

DoHeatmap(subset(srat.ge.tumor.p1, downsample = 50), features = top10.p1.ge$gene,
          group.colors = cluster.colors, assay = "RNA", slot = "scale.data") +
  theme(axis.text.y = element_text(size = 6))

VlnPlot(srat.ge.tumor.p1, features = c('nCount_RNA', 'aneuploidy_score'))

srat.ge.tumor.p1 <- AddMetaData(srat.ge.tumor.p1, metadata = srat.tumor.p1$SCT_snn_res.1, col.name = 'clust_snvs')

srat.ge.tumor.p1 <- RunPCA(srat.ge.tumor.p1, npcs = 50)

elbow <- ElbowPlot(srat.ge.tumor.p1, reduction = "pca", ndims = 50)
elbow #geom_vline(xintercept = c(seq(from = 0, to = 50, by = 1)))

DimHeatmap(srat.ge.tumor.p1, dims = (6:18), cells = 100)

srat.ge.tumor.p1 <- RunUMAP(srat.ge.tumor.p1, dims = 1:CFG$ndims.p1, assay = 'RNA')

srat.ge.tumor.p1 <- FindNeighbors(srat.ge.tumor.p1, dims = 1:CFG$ndims.p1)
## Computing nearest neighbor graph
## Computing SNN

srat.ge.tumor.p1 <- FindClusters(srat.ge.tumor.p1, resolution = seq(0.4, 2, 0.2), algorithm = 1)

columns <- grep("SCT_snn_res", names(srat.snvs@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the
## cluster ids)
for (col in columns) {
  cat("\n====\n", col, ":")
  print(table(srat.ge.tumor.p1@meta.data[[col]]))
}

clus.column <- paste0("SCT_snn_res.", 1)
Idents(srat.tumor.p1) <- clus.column

clusters <- sort(unique(srat.ge.tumor.p1@meta.data[[clus.column]]))
## switching color scheme to make it look different:
cluster.colors <- Polychrome::alphabet.colors(length(clusters))
names(cluster.colors) <- as.character(clusters)

## and show result with the other umaps to get your
## bearings:
p_clus_ge_p1 <- DimPlot(srat.ge.tumor.p1, reduction = "umap", pt.size = 0.5, group.by = clus.column,
                  cols = cluster.colors, label=TRUE)
p_clus_ge_p1 | DimPlot(srat.ge.tumor.p1, group.by = 'clust_snvs', cols = cluster.colors, pt.size = 0.5, label = TRUE)
p_clus_ge_p1 | FeaturePlot(srat.ge.tumor.p1, features = 'percent_mito')

##Pathway analysis ----
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

if (websiteLive) dbs <- listEnrichrDbs()

if (websiteLive) head(dbs)

#Select databases
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
if (websiteLive) {
  enriched <- enrichr(head(de_genes_snvs_p1[de_genes_snvs_p1$cluster == 5, 'gene'], n = 30), dbs)
} #Choose genes to see enrichment

if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}

DimPlot(srat.ge.tumor.p1, group.by = 'Phase')
VlnPlot(srat.ge.tumor.p1, features = 'BIRC5')

#Tumors together ----
srat.tumor <- subset(srat.snvs, subset = is_aneuploid == 'yes')

srat.tumor <- NormalizeData(srat.tumor, normalization.method = "LogNormalize")

srat.tumor <- ScaleData(srat.tumor, features = rownames(srat.tumor), verbose = FALSE)

srat.tumor <- FindVariableFeatures(srat.tumor)


srat.tumor <- suppressMessages(SCTransform(srat.tumor, 
                                              vars.to.regres=NULL,
                                              vst.flavor='v2', ## version that solved some problems
                                              verbose=FALSE,
                                              variable.features.n=3000))

top10 <- head(VariableFeatures(srat.tumor), 10)

DefaultAssay(srat.tumor)
srat.tumor <- RunPCA(srat.tumor, npcs = 50)

elbow <- ElbowPlot(srat.tumor, reduction = "pca", ndims = 50)
elbow #geom_vline(xintercept = c(seq(from = 0, to = 50, by = 1)))

DimHeatmap(srat.tumor, dims = (6:18), cells = 100)

srat.tumor <- RunUMAP(srat.tumor, dims = 1:CFG$ndims.tumor, assay = 'RNA')

srat.tumor <- FindNeighbors(srat.tumor, dims = 1:CFG$ndims.tumor)
## Computing nearest neighbor graph
## Computing SNN

srat.tumor <- FindClusters(srat.tumor, resolution = seq(0.4, 2, 0.2), algorithm = 1)

columns <- grep("SCT_snn_res", names(srat.snvs@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the
## cluster ids)
for (col in columns) {
  cat("\n====\n", col, ":")
  print(table(srat.tumor@meta.data[[col]]))
}

clus.column <- paste0("SCT_snn_res.", 1)
Idents(srat.tumor) <- clus.column

clusters <- sort(unique(srat.tumor@meta.data[[clus.column]]))
## switching color scheme to make it look different:
cluster.colors <- Polychrome::alphabet.colors(length(clusters))
names(cluster.colors) <- as.character(clusters)

## and show result with the other umaps to get your
## bearings:
p_clus <- DimPlot(srat.tumor, reduction = "umap", pt.size = 0.5, group.by = clus.column,
                  cols = cluster.colors, label=TRUE)

p_clus
gc()

#Error in check.length(gparname) : 
#  'gpar' element 'lwd' must not be length 0
#In addition: Warning message:
#  The `add` argument of `group_by()` is deprecated as of dplyr 1.0.0.
#ℹ Please use the `.add` argument instead.
#ℹ The deprecated feature was likely used in the dplyr package.
#Please report the issue at <https://github.com/tidyverse/dplyr/issues>.
#This warning is displayed once every 8 hours.
#Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated. 

## based on this I think it's still good to choose 0.8
clus.column <- paste0("SCT_snn_res.", 1)
Idents(srat.tumor) <- clus.column

typenames <- unique(srat.tumor@meta.data$singler)
singler.colors <- Polychrome::palette36.colors(length(typenames))
names(singler.colors) <- typenames

DimPlot(srat.tumor, group.by = 'SCT_snn_res.1') | FeaturePlot(srat.tumor, features = 'aneuploidy_score')

de_snvs <- FindAllMarkers(srat.tumor, assay = "RNA",
                           min.pct = 0.1, logfc.threshold = log(1.5))
## The resulting data.frame is too unwieldy to work with,
## the top 10 genes per cluster are enough:

de_snvs %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
head(top10)

DoHeatmap(subset(srat.tumor, downsample = 50), features = top10$gene,
          group.colors = cluster.colors, assay = "RNA", slot = "scale.data") +
  theme(axis.text.y = element_text(size = 6))
VlnPlot(srat.tumor, features = c('nCount_RNA', 'aneuploidy_score'))

p_clus | DimPlot(srat.tumor, group.by = 'patient')
p_clus | FeaturePlot(srat.tumor, features = 'nCount_RNA')
p_clus | FeaturePlot(srat.tumor, features = 'aneuploidy_score')
p_clus | DimPlot(srat.tumor, group.by = 'singler')



## Frequency plot by hand
# Supongamos que tienes tus datos en srat.snvs$SCT_snn_res.1.8 y srat.lr$SCT_snn_res.1.4

# Calcula la tabla de contingencia
contingency_table <- table(srat.ge.tumor.p1$clust_snvs, srat.ge.tumor.p1$SCT_snn_res.1)

# Convierte la tabla de contingencia a porcentajes
percent_table <- prop.table(contingency_table, margin = 1) * 100

# Convierte la tabla de porcentajes a un data frame largo (long format)
percent_df <- as.data.frame(melt(percent_table))

# Renombra las columnas para que ggplot2 las entienda
colnames(percent_df) <- c("SNV_clusts", "GE_clusts", "Percentage")
percent_df$SNV_clusts <- factor(percent_df$SNV_clusts, levels = 0:8)
percent_df$GE_clusts <- factor(percent_df$GE_clusts, levels = 0:8)


clustnames <- unique(srat.ge.tumor.p1$clust_snvs)
barplot.cols <- Polychrome::palette36.colors(length(clustnames))
names(barplot.cols) <- clustnames

# Crea el barplot
ggplot(percent_df, aes(x = SNV_clusts, y = Percentage, fill = GE_clusts)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cell cycle phase per cluster",
       x = "SNV-based clusters",
       y = "Relative abundance") +
  scale_fill_manual(values = barplot.cols) +
  theme_minimal()
