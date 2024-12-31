#install.packages("devtools")
#devtools::install_bitbucket("princessmaximacenter/scutils")

##Sometimes there are problems with the Matrix package. This should solve it:
#remove.packages('Matrix')
#devtools::install_version("Matrix", version = "1.5.0", repos = "http://cran.us.r-project.org")

library(SCutils)
library(Matrix)  # needed for working with sparse matrices
library(Seurat)
library(ggplot2)
library(dplyr)  # general data manipulation
library(patchwork)  # composing multi-panel plots
library(org.Hs.eg.db)  # genome information
library(ggpubr)
library(magrittr)
library(stringr)
library(Polychrome)
library(GO.db)
library(ggraph)
library(clustree)
library(SingleR)
library(pheatmap)
library(infercnv)
library(celldex)
library(scDblFinder)

# Data loading----
CFG <- list()  # global

CFG$data_dir <- "~/LX279_gene_expression/ONT/analysis_nanopore_LX279/transcripts"  # where to read data from
CFG$output_dir <- "~/LX279_gene_expression/ONT/analysis_nanopore_LX279/transcripts"  # where to write data to

# explained later:
CFG$ndims <- 14
CFG$random_seed <- 2033
# Cutoffs used during quality control for filtering genes
# (explained later)
CFG$min_txpts <- 200
CFG$max_txpts <- 12000
CFG$max_pctmito <- 23.5

# Graphical parameters

CFG$gradient.colors <- viridis::viridis(101, direction = -1)

## NOTE: many colors change if the number of clusters
## change, so they are not fixed here.
lib <- "LX279"  # name of the library, also needed later on.
countsdir <- paste0(CFG$data_dir, "/LX279_sup_complete.transcript_expression.counts.tsv")  # what is actually on disk
counts <- read.delim(countsdir, row.names=1) #we do processing ourselves
counts <- as.matrix(counts)
counts <- as(counts, "sparseMatrix") 
#is of type 'dgCMatrix', a sparse matrix
cnts_per_cell <- Matrix::colSums(counts)  # note use of Matrix::colSums

summary(cnts_per_cell)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#253    1693    2500    3295    4010   24913 
#lower than gene_expression somehow, doesnt count intronic 
# and is more more restrictive than gene assignment

df <- data.frame(counts = cnts_per_cell)

ggplot(df, aes(x = counts)) + geom_density() + geom_rug() +
  scale_x_continuous(trans = "log2")

## Changing cell labels ----
cellnames <- colnames(counts)
colnames(counts) <- paste0(lib, "_", cellnames)


## Metadata ----
# with transcripts we have ensembl id, we can change it to HGNC (retrieved from biomart):
ensembl_2_geneID <- read.delim("~/LX279_gene_expression/ONT/analysis_nanopore_LX279/transcripts/ensembl_2_geneID")
dim(ensembl_2_geneID)
#277105      2

#Let's check how many unique isoforms do we have
ensID <- as.vector(rownames(counts))
length(ensID)
# 65126

#We want to keep only ensembl IDs that we have in our dataset, so:
ensembl_2_geneID_filtered <- ensembl_2_geneID[ensembl_2_geneID$Transcript.stable.ID %in% ensID, ]

#But some are duplicated for some reason:
without_duplicates<-unique(ensembl_2_geneID_filtered$Transcript.stable.ID)
duplicates<-which(duplicated(ensembl_2_geneID_filtered$Transcript.stable.ID))
ensembl_2_geneID_filtered <- ensembl_2_geneID_filtered[-duplicates,]
dim(ensembl_2_geneID_filtered)
#[1] 64044     2

#Some transcripts didnt have geneID
unIDed_transcripts <- setdiff(ensID, ensembl_2_geneID_filtered$Transcript.stable.ID)
length(unIDed_transcripts)
#1082

#There are also some transcripts without GeneID. We have to deal with them
temp <- data.frame(unIDed_transcripts)
temp[,'HGNC.symbol'] <- "" 
colnames(temp) <- c("Transcript.stable.ID", "HGNC.symbol")
ensembl_2_geneID_filtered <- data.frame(rbind(ensembl_2_geneID_filtered, temp))

change<-function(x) {
  if (x[2] == "") {
    x[2] = "NA"
    return(x)
  }
  return(x)
}
ensembl_2_geneID_filtered <- data.frame(t(apply(ensembl_2_geneID_filtered, 1, change)))
ensembl_2_geneID_filtered <- ensembl_2_geneID_filtered %>%
  arrange(Transcript.stable.ID)

geneID_ensembl <- paste0(ensembl_2_geneID_filtered$HGNC.symbol, "_",
                         ensembl_2_geneID_filtered$Transcript.stable.ID)

rownames(counts) <- geneID_ensembl

mitos <- grep("^mt-", rownames(counts), value = TRUE, ignore.case = TRUE)

## show them:
mitos
#[1] "MT-ND3_ENST00000361227"  "MT-ND4L_ENST00000361335" "MT-ND4_ENST00000361381"  "MT-ND1_ENST00000361390" 
#[5] "MT-ND2_ENST00000361453"  "MT-ND5_ENST00000361567"  "MT-CO1_ENST00000361624"  "MT-ND6_ENST00000361681" 
#[9] "MT-CO2_ENST00000361739"  "MT-CYB_ENST00000361789"  "MT-ATP8_ENST00000361851" "MT-ATP6_ENST00000361899"
#[13] "MT-CO3_ENST00000362079" 

percent_mito <- 100 * Matrix::colSums(counts[mitos, ])/cnts_per_cell
# note the use of Matrix::colSums again

## get the non-mito's for easy selection of 'the rest of
## the genes'
nuclear <- setdiff(rownames(counts), mitos)

log2_transcript_counts <- log2(1 + Matrix::colSums(counts[nuclear,]))

log2_feature_counts <- log2(1 + Matrix::colSums(counts[nuclear,] > 0))

## set up a data.frame:
meta <- data.frame(percent_mito = percent_mito,
                   log2_counts = log2_transcript_counts,
                   log2_features = log2_feature_counts,
                   lib = rep(lib, ncol(counts)))

## let's show how features ( = genes ) relate to number of
## transcripts
ggplot(meta, aes(x = log2_counts, y = log2_features)) + geom_point(color = "blue")

## Creating Seurat object ----
counts <- counts[nuclear, ]

srat <- Seurat::CreateSeuratObject(counts = counts, project = "nanopore",
                                   meta = meta)
## Warning: Feature names cannot have underscores ('_'), 
#replacing with dashes ('-'). This doesn't happen with nano data

## Later this we won't use the counts object anymore, so we
## should get rid of it and cleanup
rm(counts)
gc()
##used   (Mb) gc trigger   (Mb)   max used    (Mb)
#Ncells   9666999  516.3   17059325  911.1   17059325   911.1
#Vcells 185968360 1418.9  639368077 4878.0 2373940874 18111.8
srat
#An object of class Seurat 
#65113 features across 3821 samples within 1 assay
#Active assay: RNA (65113 features, 0 variable features)

#Much less features than illumina (36588 against 24740, roughly 14k less)
head(srat@meta.data)

##Removing doublets ----
#Barcode name is without the library name, so we have to add it
sop_nano <- read.delim("~/sop_input/nuc/clusters.tsv",
                       header = TRUE,
                       row.names = 1)
rownames(sop_nano) <- paste0(lib, "_", rownames(sop_nano))
rownames(sop_nano) <- gsub("-1", "", rownames(sop_nano), fixed = TRUE)

goodBc <- rownames(sop_nano)[rownames(sop_nano) %in% colnames(srat)]

srat@meta.data[goodBc, 'patient_ont'] <- sop_nano[goodBc, 'assignment']
srat@meta.data[goodBc, 'status'] <- sop_nano[goodBc, 'status']

srat <- subset(srat, subset = status == "singlet")

#Filtering cells ----
##Transcript counts ----
v <- VlnPlot(srat, "nCount_RNA", group.by = "lib")
lines <- seq(from = 0, to = 36000, by = 250)
hlines <- geom_hline(yintercept = lines, col = "grey", linetype = 1,
                     lwd = 0.1)

v + hlines
v + scale_y_continuous(trans = "log2") + hlines

## create a scatter graph
f_lin <-
  FeatureScatter(srat,
                 feature1 = "nCount_RNA",
                 feature2 = "nFeature_RNA",
                 pt.size = 0.5) + geom_hline(yintercept = CFG$max_pctmito,
                                             linetype = 2) +
  geom_vline(xintercept = c(CFG$min_txpts, CFG$max_txpts),
             linetype = 2)
f_lin
##nCount_RNA simply means number of RNA transcripts.
##In contrast, nFeature_RNA (used below) means: number of genes expressed
##(regardless of the number of transcripts per gene)

## also show the transcripts in logarithmic scale:

f_log <- f_lin + scale_x_continuous(trans = "log2") + scale_y_continuous(trans = "log2")

f_lin | f_log 

##Mitochondrial content ----
v <- VlnPlot(srat, "percent_mito", group.by = "lib")
mitolines <- seq(from = 0, to = 100, by = 5)
mitohlines <-
  geom_hline(
    yintercept = mitolines,
    col = "grey",
    linetype = 1,
    lwd = 0.1
  )
v + mitohlines

##Final selection ----
## actually set the cutoffs (best done at the beginning of
## your scripts to keep things tidy)
## plot this:

f <-
  FeatureScatter(srat,
                 feature1 = "nCount_RNA",
                 feature2 = "percent_mito",
                 pt.size = 1) + geom_hline(yintercept = CFG$max_pctmito, linetype = 2) +
  geom_vline(xintercept = c(CFG$min_txpts, CFG$max_txpts),
             linetype = 2)

f
f + scale_x_continuous(trans = "log2")

## show how many cells will we discard by applying the
## selection:
dim(srat)
#[1] 65113  3599
dim(subset(srat, nCount_RNA < CFG$min_txpts))
## [1] 65113    48
dim(subset(srat, nCount_RNA > CFG$max_txpts))
## [1] 65113    34
dim(subset(srat, percent_mito > CFG$max_pctmito))
## [1] 65113   404

## lastly, subset the seurat object and save if needed
srat <- subset(srat, subset = nCount_RNA >= CFG$min_txpts & nCount_RNA <=
                 CFG$max_txpts & percent_mito <= CFG$max_pctmito)

dim(srat)
##[1] 24740  3557

#Hemoglobin removal (should it have been done before?) ----

data(refdata_cellranger_GRCh38_3.0.0)

hb_genes <- genelists$hemo
hb_genes <- subset(ensembl_2_geneID_filtered, subset = HGNC.symbol %in% hb_genes)
hb_genes <- paste0(hb_genes$HGNC.symbol, "-", hb_genes$Transcript.stable.ID)
## alternatively: hb_genes <- lookup.hemogenes()

hb_genes <- intersect(hb_genes, rownames(srat))

hb_counts <- Matrix::colSums(srat@assays$RNA$counts[hb_genes,])

srat <- AddMetaData(srat, col.name = "log2hb_genes",
                    metadata = log2(1 + hb_counts))

srat <- AddMetaData(srat, col.name = "pct_hemo", metadata = 100 *
                      hb_counts/srat@meta.data$nCount_RNA)

## now show the previous plot along side a plot of the hemo
## content:
df <- subset(srat@meta.data)

p_ngenes <- ggplot(df, aes(x = log2_counts, y = log2_features,
                           color = pct_hemo)) + geom_point(size = 1, alpha = 1/2) +
  scale_color_gradient(
    low = "blue",
    high = "red",
    limits = c(0,
               5),
    oob = scales::squish
  )

p_pcthemo <- ggplot(df, aes(x = log2_counts, y = pct_hemo,
                            color = log2_features)) +
geom_point(size = 1, alpha = 1/2) +
  scale_color_gradient(
    low = "blue",
    high = "red",
    limits = c(2,
               16),
    oob = scales::squish
  )

p_ngenes | p_pcthemo

p <- ggplot(df, aes(x = nCount_RNA, y = nFeature_RNA)) +
  xlab("RNA Counts")+
  ylab("Features")+
  geom_point(size=0.5, color="red")
p

## and/or logarithmic:

p + scale_x_continuous(trans = "log2") + scale_y_continuous(trans = "log2")
#add anotation
ggplot(srat@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  xlab("RNA Counts")+
  ylab("Features")+
  geom_point(size=0.5, color="red")

CFG$max_pct_hemo <- 5

dim(srat)
## [1] 65113  3161
dim(subset(srat, pct_hemo > CFG$max_pct_hemo))
## [1] 24740    28
srat <- subset(srat, pct_hemo <= CFG$max_pct_hemo)

p_withouthemo <- ggplot(srat@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  xlab("RNA Counts")+
  ylab("Features")+
  geom_point(size=0.5, color="red")

p | p_withouthemo

## Saved srat ----
file <- paste0(CFG$output_dir, "/srat-cellsfiltered.rds")
saveRDS(file = file, srat)
#srat <-readRDS(file)

#Normalization and Dimensional Reduction ----
srat <- NormalizeData(srat, normalization.method = "LogNormalize")

srat <- ScaleData(srat, features = rownames(srat), verbose = FALSE)

srat <- FindVariableFeatures(srat)


srat <- suppressMessages(SCTransform(srat, 
                                     vars.to.regres=NULL,
                                     vst.flavor='v2', ## version that solved some problems
                                     verbose=FALSE,
                                     variable.features.n=3000))


##Saved srat ----

file<-paste0(CFG$output_dir, "/srat-normalized.rds")
saveRDS(file=file,srat)
#srat <- readRDS(file = file)
##To read the object back in, do
#srat <- readRDS("~/analysis_illumina_LX279/srat-normalized.rds")

top10 <- head(VariableFeatures(srat), 10)

p.varfeat <- VariableFeaturePlot(srat, selection.method = 'sct')
p.varfeat <- LabelPoints(plot = p.varfeat, points = top10, repel = TRUE)

## show plot without labels:
p.varfeat

DefaultAssay(srat)  # just check that it has 'SCT'
## [1] "SCT"
srat <- RunPCA(srat, npcs = 50)

DimPlot(srat, reduction = "pca", pt.size = 1)

elbow <- ElbowPlot(srat, reduction = "pca", ndims = 50)
elbow #geom_vline(xintercept = c(seq(from = 0, to = 50, by = 1)))

DimHeatmap(srat, dims = (1:12), cells = 100)

DimHeatmap(srat, dims = 12 + (1:12), cells = 100)

DimHeatmap(srat, dims = 24 + (1:12), cells = 100)

CFG$ndims <- 25

data(cc.genes)  # load the cell cycle genes (also available from SCutils' genelists)
data(refdata_cellranger_GRCh38_3.0.0)

s.genes <- subset(ensembl_2_geneID_filtered, subset = HGNC.symbol %in% cc.genes$s.genes)
s.genes <- paste0(s.genes$HGNC.symbol, "-", s.genes$Transcript.stable.ID)

g2m.genes <- subset(ensembl_2_geneID_filtered, subset = HGNC.symbol %in% cc.genes$g2m.genes)
g2m.genes <- paste0(g2m.genes$HGNC.symbol, "-", g2m.genes$Transcript.stable.ID)

srat <- CellCycleScoring(object = srat, s.features = s.genes,
                         g2m.features = g2m.genes, assay = "RNA")
## Warning: The following features are not present in the object:
#MLF1IP, not searching for symbol synonyms

stress_genes <- genelists$stress  # or lookup.stressgenes()
stress_genes <- subset(ensembl_2_geneID_filtered, subset = HGNC.symbol %in% stress_genes)
stress_genes <- paste0(stress_genes$HGNC.symbol, "-", stress_genes$Transcript.stable.ID)
srat <- AddModuleScore(srat, features = list(stress = stress_genes),
                       name = "stress", assay = "RNA")
## Warning: The following features are not present in the object:
#COMP, PRKN, AGR2, ERN2, FGF21, HSPB8, PACRG, BHLHA15, MIR199A1,
#MMP24-AS1-EDEM2, not searching for symbol synonyms

srat <- RunUMAP(srat, dims = 1:CFG$ndims)

p_phase <- DimPlot(srat, reduction = "umap", pt.size = 0.5, group.by = "Phase")

p_phase
#2 weird clusters of S and another of G2M

CFG$gradient.colors <- viridis::viridis(101, direction = -1)

norm1<-FeaturePlot(srat, pt.size = 0.5,
                   feature = c("nCount_RNA", "percent_mito",
                               "nFeature_RNA", "stress1"), order = TRUE)
norm1

# Deconfounding ----

## In this part we remove features, not cells
Scor <- metadataCorrelations(srat, "S.Score")

G2Mcor <- metadataCorrelations(srat, "G2M.Score")

additional <- derivedCellcycleGenes(Scor = Scor, Sgenes = s.genes,
                                    G2Mcor = G2Mcor, G2Mgenes = g2m.genes)
additional
additional$plot[[1]] | additional$plot[[2]] | additional$plot[[4]]

## also get the other deconfounders:
#stress_genes <- genelists$stress #OR: SCutils::lookup.stressgenes()
#hb_genes <- genelists$hemo  # as before
fem_genes <- genelists$female
fem_genes <- subset(ensembl_2_geneID_filtered, subset = HGNC.symbol %in% fem_genes)
fem_genes <- paste0(fem_genes$HGNC.symbol, "-", fem_genes$Transcript.stable.ID)

male_genes <- genelists$male
male_genes <- subset(ensembl_2_geneID_filtered, subset = HGNC.symbol %in% male_genes)
male_genes <- paste0(male_genes$HGNC.symbol, "-", male_genes$Transcript.stable.ID)

## combine these sets (also including stress_genes and
## hemoglobin genes)
remove <- unique(c(s.genes, g2m.genes, additional$S.derived,
                   additional$G2M.derived, stress_genes, hb_genes,
                   fem_genes, male_genes))

## check how many we will loose now:
length(intersect(VariableFeatures(srat), remove))
## [1] 202
## do the removal, for both the 'RNA' (i.e. LogNorm) and
## 'SCT' assays:
VariableFeatures(srat, assay = "RNA") <-
  setdiff(VariableFeatures(srat,assay = "RNA"), remove)

VariableFeatures(srat, assay = "SCT") <-
  setdiff(VariableFeatures(srat,assay = "SCT"), remove)
srat <- RunPCA(srat, verbose = TRUE, npcs = 50)
srat <- RunUMAP(srat, dims = 1:CFG$ndims)

p_type <- DimPlot(srat, reduction = "umap", pt.size = 0.5)

p_phase <- DimPlot(srat, reduction = "umap", pt.size = 0.5, group.by = "Phase")

## use patchwork again to get a 2 x 2 plot
p_phase

## and also the continuous variables:
FeaturePlot(
  srat,
  pt.size = 0.5,
  feature = c("nCount_RNA", "percent_mito",
              "nFeature_RNA", "stress1"),
  order = TRUE
)

#Clustering ----
srat <- FindNeighbors(srat, dims = 1:CFG$ndims)
## Computing nearest neighbor graph
## Computing SNN

srat <- FindClusters(srat, resolution = 0.4, algorithm = 1)
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 1296
## Number of edges: 37806
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9242
## Number of communities: 8
## Elapsed time: 0 seconds

srat <- FindClusters(srat, resolution = 0.8, algorithm = 1)
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 1296
## Number of edges: 37806
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8690
## Number of communities: 12
## Elapsed time: 0 seconds

srat <- FindClusters(srat, resolution = 1.2, algorithm = 1)
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 1296
## Number of edges: 37806
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8343
## Number of communities: 15
## Elapsed time: 0 seconds
srat <- FindClusters(srat, resolution = 1.6, algorithm = 1)
srat <- FindClusters(srat, resolution = 2.0, algorithm = 1)
columns <- grep("SCT_snn_res", names(srat@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the
## cluster ids)
for (col in columns) {
  cat("\n====\n", col, ":")
  print(table(srat@meta.data[[col]]))
}
##===
#SCT_snn_res.0.4 :
#  0   1   2   3   4   5   6   7   8   9  10  11 
#790 634 518 281 278 270 246 233  94  90  36  33 

#===
#  SCT_snn_res.0.8 :
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
#512 467 418 368 282 281 270 246 153 110  90  80  63  56  38  36  33 

#===
#  SCT_snn_res.1.2 :
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19 
#511 372 281 280 279 270 261 207 154 144 137 110 101  90  80  63  56  38  36  33 

## let's plot the different clusterings side by side to see
## the differences:

p_clusbyreso <- list()  # set up the list to hold the intermediate plots

for (reso in c("0.8", "1.2", "1.6", "2")) {
  
  clus.column <- paste0("SCT_snn_res.", reso)
  clusters <- sort(unique(srat@meta.data[[clus.column]]))
  
  ## switching color scheme to make it look different:
  clus.colors <- Polychrome::alphabet.colors(length(clusters))
  ## clus.colors <-
  ## Polychrome::kelly.colors(length(clusters))
  names(clus.colors) <- as.character(clusters)
  
  p <- DimPlot(srat, pt.size = 0.5, group.by = clus.column, cols = clus.colors) +
    labs(title = reso)
  p_clusbyreso[[reso]] <- p
}

## use patchwork to plot them side-by-side:
p_phase | p_clusbyreso[[2]] | p_clusbyreso[[3]]

## set the identities to the clusters coming from the 0.8
## resolution:

clus.column <- paste0("SCT_snn_res.", 2)
Idents(srat) <- clus.column

clusters <- sort(unique(srat@meta.data[[clus.column]]))
## switching color scheme to make it look different:
cluster.colors <- carto.pal(pal1 = "multi.pal", n1 =(length(clusters)))
names(cluster.colors) <- as.character(clusters)

## and show result with the other umaps to get your
## bearings:
p_clus <- DimPlot(srat, reduction = "umap", pt.size = 0.5, group.by = clus.column,
                  label=TRUE)

p_clus
rm(p_clusbyreso)
gc()

## Using clustree ----
#srat <- FindClusters(srat, resolution = seq(0.4, 2, 0.2), algorithm = 1)
#clustree::clustree(srat, prefix = "SCT_snn_res.", node_size = 2, edge_width = 0.5)
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
clus.column <- paste0("SCT_snn_res.", 2)
Idents(srat) <- clus.column

## scDblFinder ----
sce <- scDblFinder(GetAssayData(srat, assay = "RNA", slot = "counts"), 
                   iter = 10, 
                   includePCs = 14,
                   dims = 20)

srat$scDblFinder.score <- sce$scDblFinder.score
srat$scDblFinder.class <- sce$scDblFinder.class

p_sc.class <- DimPlot(srat, group.by = 'scDblFinder.class') +
  ggtitle('scDblrFinder class for LR')
p_sc.score <- FeaturePlot(srat, features = 'scDblFinder.score') +
  ggtitle('scDblFinder score for LR')
p_sc.class | p_sc.score

#Enrichment analysis ----
##Differently expressed genes ----
## first step: find the differentially expressed genes this
## takes ~ 7min
de_transcripts <- FindAllMarkers(srat, assay = "RNA", only.pos = TRUE,
                           min.pct = 0.1, logfc.threshold = log(1.5))
## The resulting data.frame is too unwieldy to work with,
## the top 10 genes per cluster are enough:

de_transcripts %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
head(top10)

# A tibble: 6 × 7
# Groups:   cluster [1]
#p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene  
#<dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr> 
#  1 0               2.51 0.928 0.183 0         0       CAMK4 
#2 9.81e-269       2.57 0.973 0.399 2.43e-264 0       NDFIP1
#3 3.21e-248       2.53 0.996 0.59  7.94e-244 0       ZEB1  
#4 9.44e-222       1.94 0.951 0.333 2.34e-217 0       RPS4Y1
#5 3.32e-215       2.37 0.756 0.157 8.21e-211 0       NELL2 
#6 1.62e-192       1.82 0.83  0.261 4.01e-188 0       ITPKB 

DoHeatmap(subset(srat, downsample = 50), features = top10$gene,
          group.colors = cluster.colors, assay = "RNA", slot = "scale.data") +
  theme(axis.text.y = element_text(size = 6))
DimPlot(srat,group.by = "patient")
## Enrichment analysis finally ----
## compareClusters needs a list() of clusters, not the
## grouped_df returned by dplyr::top_n. Convert the top10
## to a list:

df <- as.data.frame(top10)
clusters <- split(df$gene, df$cluster)

## show contents:
show(clusters)
## compareCluster crashes on gene symbols for some reason,
## so let's use EnsEMBL ids. The SCutils data() statement
## also attaches translation table 'gene_name2ensembl'
## which can be used to do the conversion

## create a little function that does the conversion for
## us.  We omit the original names for readability, and
## also ignore the NA's
lookup_ensg_id <- function(id) {
  na.omit(unname(gene_name2ensembl[id]))
}
## call it:
clusters <- lapply(clusters, lookup_ensg_id)

## show result:
show(clusters)
## Lastly, run it. This takes ~ 2 min
#enrichment <- compareCluster(clusters, fun = "enrichGO", OrgDb = "org.Hs.eg.db",
#                             keyType = "ENSEMBL", ont = "BP",
#                            minGSSize = 3, pAdjustMethod = "BH",
#                          pvalueCutoff = 0.01, qvalueCutoff = 0.05,
#                           readable = TRUE)

## so save it:
#file <- paste0(CFG$output_dir, "/enrichment.rds")
#saveRDS(file = file, enrichment)

## show results:
#clusterProfiler::dotplot(enrichment, showCategory = 4, 
#title = "GO biological process", font = 5)

p_patient <- DimPlot(srat, reduction = "umap", pt.size = 0.5,
                     group.by = "patient")
p_status <- DimPlot(srat, reduction = "umap", pt.size = 0.5,
                    group.by = "status")

p_patient | p_clus

##Saved srat ----
file <- paste0(CFG$output_dir, "/enrichment.rds")
saveRDS(file = file, srat)
#srat.lr <- readRDS(paste0(CFG$output_dir, "/enrichment.rds"))

#Cell identification ----
# (down)load the reference data:
#we have no reference for transcripts, so we just project cellID from genes onto
#our transcript clustering to see differences
hpca <- HumanPrimaryCellAtlasData()
hpca_names <- subset(ensembl_2_geneID_filtered, subset = hpca@NAMES %in% ensembl_2_geneID_filtered$HGNC.symbol)
hpca_names <- arrange(hpca_names, HGNC.symbol)
hpca_names <- paste0(hpca_names$HGNC.symbol, "-", hpca_names$Transcript.stable.ID)
hpca@NAMES <- hpca_names

srat.gene <-
  readRDS("~/analysis_nanopore_LX279/srat-cellsidentified.rds")

srat <- AddMetaData(srat, colnames(srat), col.name = 'bcs')
srat.gene <-
  AddMetaData(srat.gene, colnames(srat.gene), col.name = 'bcs')
merged.metadata <- merge(srat@meta.data,
                         srat.gene@meta.data[, c('singler', 'bcs')], by =
                           'bcs')
srat <- srat[, colnames(srat) %in% merged.metadata$bcs]
srat@meta.data$singler <- merged.metadata$singler.x

##Add metadata
typenames <- unique(srat@meta.data$singler)
singler.colors <- Polychrome::palette36.colors(length(typenames))
names(singler.colors) <- typenames

p_singler <- DimPlot(srat, reduction = "umap", pt.size = 0.5,
                     group.by = "singler", cols = singler.colors) + ggtitle(
                       "singleR LR")
p_pat <- DimPlot(srat, reduction = "umap", pt.size = 0.5,
                 group.by = "patient")

p_singler | p_pat

typenames <- unique(srat$singler)
type2.colors <- Polychrome::alphabet.colors(length(typenames))
names(type2.colors) <- typenames

p_type2 <- DimPlot(srat, reduction = "umap", pt.size = 0.5, group.by = "singler",
                   cols = singler.colors)

## show the final result, this time together with the other
## characteristics just because we can :-)
p_clus | p_type2

file <- paste0(CFG$output_dir, "/srat-cellsidentified.rds")
saveRDS(file = file, srat)
#srat <- readRDS(file)
#Tumor cells identification ----
## overview of numbers per celltype:
table(srat@meta.data$SCT_snn_res.0.8)

## which cells do we think we can trust as being healthy:
ref_types <- c("0","1","7","8","9","11","14","15")
maybe_tumor <- !(srat@meta.data$SCT_snn_res.0.8 %in% ref_types)
table(maybe_tumor)

## infercnv needs them as an external data file in a
## specific format, create that:
df <- data.frame(cell = colnames(srat), type = srat@meta.data$SCT_snn_res.0.8)

celltypes_file <- paste0(CFG$output_dir, "/celltypes.txt")
write.table(df, file = celltypes_file, sep = "\t", quote = FALSE,
            na = "", row.names = FALSE, col.names = FALSE)
## infercnv also needs to gene coordinates, in the
## tab-delimited format genename chrNN START END.  We can
## use the `genecoords` object from SCutils for this.

## this data needs some cleaning. Use only the canonical
## chromosomes, sort them, get rid of duplicates and put
## 'chr' in front of the chromosome:
wanted.chr <- c(as.character(1:22), "X", "Y")  # skip MT and KI2707 etc.
geneorder <- genecoords[genecoords$chr %in% wanted.chr, ]
geneorder$numeric_chr <- as.numeric(geneorder$chr)
## Warning: NAs introduced by coercion
geneorder[geneorder$chr %in% "X", "numeric_chr"] <- 30
geneorder[geneorder$chr %in% "Y", "numeric_chr"] <- 40
geneorder <- geneorder[order(geneorder$numeric_chr, geneorder$start),
]
geneorder <- geneorder[!duplicated(geneorder$gene_name), ]
geneorder <- with(geneorder, data.frame(gene = gene_name,
                                        chromosome = paste0("chr",chr),
                                        start = start,
                                        end = end))
## @@@@ WRONG DIR, FIX: geneorder_file <-
## paste0(CFG$data_dir, '/geneorder.txt')
geneorder_file <- paste0(CFG$output_dir, "/geneorder.txt")
write.table(geneorder, file = geneorder_file, sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE, na = "")

raw_counts <- GetAssayData(srat, assay = "RNA", slot = "counts")

common <- intersect(rownames(raw_counts), geneorder$gene)

## set up the infercnv object with all necessary
## information:
infercnv_obj <- infercnv::CreateInfercnvObject(
  delim = "\t",
  raw_counts_matrix = raw_counts[common,],
  annotations_file = celltypes_file,
  gene_order_file = geneorder_file,
  ref_group_name = ref_types
)
## infercnv writes its output to disk. It is a fair amount,
## so let's create a separate directory for that and write
## it there:
outdir <- paste0(CFG$output_dir, "/infercnv")
dir.create(outdir)
## running takes can take up to 10 min, best start this
## before the lecture:
infercnv_obj <- infercnv::run(infercnv_obj, out_dir = outdir,
                              cluster_by_groups = TRUE, denoise = TRUE, analysis_mode = "samples",
                              num_threads = 3, output_format = "png", window_length = 101,
                              save_final_rds = TRUE, plot_steps = FALSE, cutoff = 0.1,
                              sd_amplifier = 1.5, HMM = FALSE)
### @@@ this has changed in version 1.14, so the names of
### files are different again! FIX THIS
refexp <- as.matrix(read.delim(file = paste0(outdir, "/infercnv.references.txt"),
                               sep = " "))

obsexp <- as.matrix(read.delim(file = paste0(outdir, "/infercnv.observations.txt"),
                               sep = " "))

## first scale them so that SD is the 'unit' of Modif.
## expression:
obs.scaled <- scale(obsexp)
ref.scaled <- scale(refexp)

## trick by Jurrian de Kanter
.count_genes_in_runs <- function(x, length = 50, SDs = 2) {
  ## count total number of genes in stretches longer than
  ## LENGTH exceeding SDs std.devs from mean.
  runs <- rle(x > SDs)
  up <- sum(runs$lengths[runs$lengths > length & runs$values])
  runs <- rle(x < -SDs)
  down <- sum(runs$lengths[runs$lengths > length & runs$values])
  ## (don't use abs(x) > because pos and neg stretches
  ## might coalesce)
  sum(up, down)
}  # .count_genes_in_runs

ScoreCNV <- function(adj_expr, length = 70, SDs = 1.5) {
  apply(adj_expr, 2, .count_genes_in_runs, length = length,
        SDs = SDs)
}  ## ScoreCNV
len <- 50  # also try e.g. 30 or 70
SDs <- 1.3  # also try e.g. 1.5 or 2
obsscore <- ScoreCNV(obs.scaled, length = len)
refscore <- ScoreCNV(ref.scaled, SDs = SDs)

ttest <- t.test(obsscore, refscore)

plot(
  density(obsscore),
  main = sprintf("Aneuploidy score length=%d SDs=%.1f",
                 len, SDs),
  xlab = sprintf("t=%.2f", ttest$statistic)
)
rug(obsscore, ticksize = 0.01)
rug(refscore, ticksize = -0.01, col = "red")
lines(density(refscore), col = "red")
legend(x = "topright", legend = c("tumor", "reference"),
       col = c("black","red"), pch = NULL, lty = 1, bty = "n")
allscore <- c(obsscore, refscore)
allscore <- allscore[colnames(srat)]  # put in order of srat object.
srat <- AddMetaData(srat, allscore, col.name = "aneuploidy_score")

cutoff <- 600  # based on the density plot

aneuploid <- ifelse(srat@meta.data$aneuploidy_score < cutoff,
                    "no", "yes")
srat <- AddMetaData(srat, col.name = "is_aneuploid", aneuploid)

## Let's show all we have. If it's too much information,
## leave some of it out.

p_pat <- DimPlot(srat, reduction = "umap", pt.size = 0.5,
                 group.by = "patient")

p_celltype <- DimPlot(srat, reduction = "umap", pt.size = 0.5,
                      group.by = "singler", cols=singler.colors)

p_aneup <- FeaturePlot(srat, reduction = "umap", pt.size = 0.5,
                       features = "aneuploidy_score", cols = CFG$gradient.colors,
                       order = TRUE)

p_isaneup <- DimPlot(srat, reduction = "umap", pt.size = 0.5,
                     group.by = "is_aneuploid",cols = c("red","green"))

p_count <- FeaturePlot(srat, reduction = "umap", pt.size = 0.5,
                       features = "nCount_SCT", cols = CFG$gradient.colors,
                       order = TRUE, label = TRUE) + ggtitle("Counts per cell")

( p_pat | p_celltype)/(p_aneup | p_isaneup)

p_count | p_aneup

file<-paste0(CFG$output_dir, "/srat-aneuploidy.rds")
saveRDS(file=file,srat)

#Trash ----
## Per patient ----
srat0 <- subset(srat, subset = patient == "0")
srat1 <- subset(srat, subset = patient == "1")
singlersrat0 <- SingleR(
  test = GetAssayData(srat0, assay = "RNA",
                      slot = "data"),
  ref = hpca,
  labels = hpca$label.main
)
singlersrat1 <- SingleR(
  test = GetAssayData(srat1, assay = "RNA",
                      slot = "data"),
  ref = hpca,
  labels = hpca$label.main
)
#different patients
scor0 <- plotScoreHeatmap(singlersrat0, show.labels = TRUE, max.labels = 100,
                          show.pruned = FALSE, order.by = "clusters",
                          clusters = srat0@meta.data[[clus.column]],
                          annotation_colors = list(Clusters = cluster.colors))
scor1 <- plotScoreHeatmap(singlersrat1, show.labels = TRUE, max.labels = 100,
                          show.pruned = FALSE, order.by = "clusters",
                          clusters = srat1@meta.data[[clus.column]],
                          annotation_colors = list(Clusters = cluster.colors))
scor0
scor1

# With fine col
finesrat<-srat

fine_hpca <- SingleR(
  test = GetAssayData(finesrat, assay = "RNA",
                      slot = "data"),
  ref = hpca,
  labels = hpca$label.fine
)

finesrat@meta.data$singler <- fine_hpca$first.labels

cbind(table(fine_hpca$first.labels))

finetypenames <- unique(finesrat@meta.data$singler)
singler.colors <- Polychrome::palette36.colors(length(finetypenames))
names(singler.colors) <- finetypenames

p_fine <- DimPlot(finesrat, reduction = "umap", pt.size = 0.5,
                  group.by = "singler", cols = singler.colors) +
  ggtitle("fine singler")

## show it:

p_fine | p_singler | p_clus