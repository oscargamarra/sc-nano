#install.packages("devtools")
#devtools::install_bitbucket("princessmaximacenter/scutils")

##Sometimes there are problems with Matrix package, this should solve it:
#remove.packages('Matrix')
#devtools::install_version("Matrix", version = "1.5.0", repos = "http://cran.us.r-project.org")

#Packages used in this script
library(SCutils)
library(Seurat)
library(ggplot2)
library(Matrix)  # needed for working with sparse matrices
library(dplyr)  # general data manipulation
library(patchwork)  # composing multi-panel plots
library(org.Hs.eg.db)  # genome information
library(ggpubr)
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

#Overall list that to include shortcuts to data used
CFG <- list()  # global

CFG$data_dir <- "~/LX279_gene_expression/illumina"  # where to read data from
CFG$output_dir <- "~/LX279_gene_expression/illumina/analysis_illumina_LX279"  # where to write data to

# explained later:
CFG$ndims <- 15
CFG$random_seed <- 2033
# Cutoffs used during quality control for filtering genes
# (explained later)
CFG$min_txpts <- 1800
CFG$max_txpts <- 60000
CFG$max_pctmito <- 12.5

# Graphical parameters

CFG$gradient.colors <- viridis::viridis(101, direction = -1)
## NOTE: many colors change if the number of clusters
## change, so they are not fixed here.

lib <- "LX279"  # name of the library, also needed later on.
countsdir <- paste0(CFG$data_dir, "/filtered_matrix_LX279_short_reads")  # what is actually on disk
counts <- Read10X(data.dir =  countsdir)  # this
#is of type 'dgCMatrix', a sparse matrix
cnts_per_cell <- colSums(counts)  # note use of Matrix::colSums

summary(cnts_per_cell)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     501   8909   13254   17149   21828   103466

df <- data.frame(counts = cnts_per_cell)


ggplot(df, aes(x = counts)) + geom_density() + geom_rug() + scale_x_continuous(trans = "log2")

## Changing cell labels ----
cellnames <- colnames(counts)
colnames(counts) <- paste0(lib, "_", cellnames)

## Metadata ----
# get the mitochondrial genes based on their name:
mitos <- grep("^mt-", rownames(counts), value = TRUE, ignore.case = TRUE)

## show them:
mitos
##  [1] "MT-ND1"  "MT-ND2"  "MT-CO1"  "MT-CO2"  "MT-ATP8"
#"MT-ATP6" "MT-CO3"  "MT-ND3"  "MT-ND4L" "MT-ND4"  "MT-ND5"  "MT-ND6"  "MT-CYB"

percent_mito <- 100 * colSums(counts[mitos, ])/cnts_per_cell
# note the use of Matrix::colSums again

## get the non-mito's for easy selection of 'the rest of
## the genes'
nuclear <- setdiff(rownames(counts), mitos)

log2_transcript_counts <- log2(1 + Matrix::colSums(counts[nuclear,
]))

log2_feature_counts <- log2(1 + Matrix::colSums(counts[nuclear,
] > 0))

## set up a data.frame:
meta <- data.frame(percent_mito = percent_mito, log2_counts = log2_transcript_counts,
                   log2_features = log2_feature_counts, lib = rep(lib, ncol(counts)))

## let's show how features ( = genes ) relate to number of
## transcripts
ggplot(meta, aes(x = log2_counts, y = log2_features)) + geom_point(color = "blue")

## Creating Seurat object ----
counts <- counts[nuclear, ]

srat <- Seurat::CreateSeuratObject(counts = counts, project = "ALL",
                                   meta = meta)
## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

## Later this we won't use the counts object anymore, so we
## should get rid of it and cleanup
rm(counts)
gc()
##used  (Mb) gc trigger  (Mb)  max used  (Mb)
##Ncells  6453718 344.7   12043002 643.2   9696168 517.9
##Vcells 60713517 463.3  113996091 869.8 101583596 775.1
srat
##An object of class Seurat 
##36588 features across 4197 samples within 1 assay 
##Active assay: RNA (36588 features, 0 variable features)
##2 layers present: counts, data
str(srat)
head(srat@meta.data)

##Removing doublets ----


#We load souporcell output
sop <- read.delim(paste0(CFG$data_dir, "/clusters_k2_ilmn.tsv"), header=TRUE, row.names = 1)
rownames(sop) <- paste0(lib,"_",rownames(sop))

#Souporcell demultiplexes ALL cells (good and bad quality). We select only our filtered cells
goodBc <- rownames(sop)[rownames(sop) %in% colnames(srat)]

#Add to metadata assignment and status of each cell
srat@meta.data[goodBc,'patient'] <- sop[goodBc,'assignment']
srat@meta.data[goodBc,'status'] <- sop[goodBc,'status']

#Discard heterogenotypic doublets
srat <- subset(srat, subset = status == "singlet")

#We check whether everything is correct
unique(srat$status)
##[1] "singlet"

#Filtering cells----
##Transcript counts ----
v <- VlnPlot(srat, "nCount_RNA", group.by = "lib")
lines <- seq(from = 0, to = 1e+05, by = 2000)
hlines <- geom_hline(yintercept = lines, col = "grey", linetype = 1,
                     lwd = 0.1)

v + hlines  # plot now has those lines added
v + scale_y_continuous(trans = "log2") + hlines 

## Create a scatter graph. Play with thresholds for correct assessment
f_lin <-
  FeatureScatter(srat,
                 feature1 = "nCount_RNA",
                 feature2 = "nFeature_RNA",
                 pt.size = 0.5) + geom_hline(yintercept = CFG$max_pctmito,
                                             linetype = 2) +
  geom_vline(xintercept = c(CFG$min_txpts, CFG$max_txpts),
             linetype = 2)
f_lin

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
## Actually set the cutoffs (best done at the beginning of
## Your scripts to keep things tidy)
CFG$min_txpts
CFG$max_txpts
CFG$max_pctmito

## plot this:
f <-
  FeatureScatter(srat,
                 feature1 = "nCount_RNA",
                 feature2 = "percent_mito",
                 pt.size = 0.5) + geom_hline(yintercept = CFG$max_pctmito, linetype = 2) +
  geom_vline(xintercept = c(CFG$min_txpts, CFG$max_txpts),
             linetype = 2)

f
f + scale_x_continuous(trans = "log2")

## show how many cells will we discard by applying the selection:
dim(srat)
## [1] 36588  4017
dim(subset(srat, nCount_RNA < CFG$min_txpts))
## [1] 36588   124
dim(subset(srat, nCount_RNA > CFG$max_txpts))
## [1] 36588    52
dim(subset(srat, percent_mito > CFG$max_pctmito))
## [1] 36588    57
## lastly, subset the seurat object and save if needed

srat <- subset(srat, subset = nCount_RNA >= CFG$min_txpts & nCount_RNA <=
                 CFG$max_txpts & percent_mito <= CFG$max_pctmito)

dim(srat)
##[1] 36588  3751
##nCount_RNA simply means number of RNA transcripts.
##In contrast, nFeature_RNA (used below) means: number of genes expressed
##(regardless of the number of transcripts per gene)

#Hemoglobin removal (should it have been done before?) ----
#Load data for human genes
data(refdata_cellranger_GRCh38_3.0.0)

#Select hemoglobin genes
hb_genes <- genelists$hemo
## alternatively: hb_genes <- lookup.hemogenes()

hb_genes <- intersect(hb_genes, rownames(srat))

hb_counts <- Matrix::colSums(srat@assays$RNA$counts[hb_genes, ])

#Add it in metadata as log2 (+1 for cells with 0 counts)
srat <- AddMetaData(srat, col.name = "log2hb_genes", metadata = log2(1 +
                                                                       hb_counts))
#Also in relation to total number of counts, as a percentage
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
                            color = log2_features)) + geom_point(size = 1, alpha = 1/2) +
  scale_color_gradient(
    low = "blue",
    high = "red",
    limits = c(2,
               16),
    oob = scales::squish
  )

show(p_ngenes | p_pcthemo)


p <- ggplot(df, aes(x = nCount_RNA, y = nFeature_RNA)) +
  xlab("RNA Counts")+
  ylab("Features")+
  geom_point(size=0.5, color="red")

## ggplot2 plot objects won't plot automatically inside a
## loop.  We therefore have explicitly call it with
## show():

show(p)

## and/or logarithmic:

p + scale_x_continuous(trans = "log2") + scale_y_continuous(trans = "log2")

CFG$max_pct_hemo <- 5

#See how many cells are being removed with this threshold
dim(srat)
## [1] 36588  3680
dim(subset(srat, pct_hemo > CFG$max_pct_hemo))
## [1] 36588    24

#Aaaaand subset
srat <- subset(srat, pct_hemo <= CFG$max_pct_hemo)

p_withouthemo <- ggplot(srat@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  xlab("RNA Counts")+
  ylab("Features")+
  geom_point(size=0.5, color="red")

#Compare plots here. See if outliers have been removed after filtering
p | p_withouthemo

## maybe save
file <- paste0(CFG$output_dir, "/srat-cellsfiltered.rds")
saveRDS(file = file, srat)

#Normalization and Dimensional Reduction ----
#Two types of normalization. First, Lognormalize and scale
srat <- NormalizeData(srat, normalization.method = "LogNormalize")

srat <- ScaleData(srat, features = rownames(srat), verbose = FALSE)

srat <- FindVariableFeatures(srat)

#Seurat also has SCTransform, which is a more complex and creates a new assay.
#SCTransform uses gene abundance info to make a more robust normalization, 
srat <- suppressMessages(SCTransform(srat, 
                                     vars.to.regres=NULL,
                                     vst.flavor='v2', ## version that solved some problems
                                     verbose=FALSE,
                                     variable.features.n=3000))

##Last saved srat ----
file<-paste0(CFG$output_dir, "/srat-normalized.rds")
saveRDS(file=file,srat)

##To read the object back in, do
#srat <- readRDS(file)

#Using SCTransform turns the default assay from your srat onto 'SCT'
#That is why this VariableFeatures will be different from the previous one
srat <- FindVariableFeatures(srat)

top10 <- head(VariableFeatures(srat), 10)

p.varfeat <- VariableFeaturePlot(srat, selection.method = 'sct')

## show plot without labels:
p.varfeat

DefaultAssay(srat)  # just check that it has 'SCT'
## [1] "SCT"

#Calculate Principal Components
srat <- RunPCA(srat, npcs = 50)
DimPlot(srat, reduction = "pca", pt.size = 1)

#See SD of differents PCs. Select the number of PCs where the 'elbow' is
elbow <- ElbowPlot(srat, reduction = "pca", ndims = 50)
elbow #geom_vline(xintercept = c(seq(from = 0, to = 50, by = 1)))

#We can also select number of PCs by looking at how well our PCs allow us
#to distinguish between different sets of cells. Clean separation = include PC.
DimHeatmap(srat, dims = (1:24), cells = 100)

DimHeatmap(srat, dims = 24 + (1:12), cells = 100)

CFG$ndims <- 30

#Load cell cycle genes
data(cc.genes)  # load the cell cycle genes (also available from SCutils' genelists)

#Calculate a module score of cell cycle genes per cell
srat <- CellCycleScoring(object = srat, s.features = cc.genes$s.genes,
                         g2m.features = cc.genes$g2m.genes, assay = "RNA")

## Warning: The following features are not present in the object:
#MLF1IP, not searching for symbol synonyms

stress_genes <- genelists$stress  # or lookup.stressgenes()

srat <- AddModuleScore(srat, features = list(stress = stress_genes),
                       name = "stress", assay = "RNA")

## Warning: The following features are not present in the object:
#COMP, PRKN, AGR2, ERN2, FGF21, HSPB8, PACRG, BHLHA15, MIR199A1,
#MMP24-AS1-EDEM2, not searching for symbol synonyms

#We do UMAP for visualization using the selected dims from PCA
srat <- RunUMAP(srat, dims = 1:CFG$ndims)

p_type <- DimPlot(srat, reduction = "umap", pt.size = 0.5)

p_phase <- DimPlot(srat, reduction = "umap", pt.size = 0.5, group.by = "Phase")

(p_type | p_phase)

CFG$gradient.colors <- viridis::viridis(101, direction = -1)

norm1<-FeaturePlot(srat, pt.size = 0.5,
                   feature = c("nCount_RNA", "percent_mito",
                                "nFeature_RNA", "stress1"), order = TRUE)
norm1

# Removing deconfounders ----

## In this part we remove features, not cells.
#Calculate correlation of genes with cell metadata over all cells
Scor <- metadataCorrelations(srat, "S.Score")

G2Mcor <- metadataCorrelations(srat, "G2M.Score")

#the object 'additional' includes non-canonical genes that correlate with cc genes
additional <- derivedCellcycleGenes(Scor = Scor, Sgenes = genelists$s.genes,
                                    G2Mcor = G2Mcor, G2Mgenes = genelists$g2m.genes)
additional
additional$plot[[1]] | additional$plot[[2]] | additional$plot[[4]]

## also get the other deconfounders:
stress_genes <- genelists$stress #OR: SCutils::lookup.stressgenes()
hb_genes <- genelists$hemo  # as before
fem_genes <- genelists$female
male_genes <- genelists$male
srat <- AddModuleScore(srat, features = list(hb = hb_genes,
                                             fem = fem_genes,
                                             male = male_genes),
                       name = c('hb', 'fem', 'male'), assay = "RNA")

FeaturePlot(srat, features = 'fem2') | DimPlot(srat, group.by = 'patient')
## combine these sets (also including stress_genes and
## hemoglobin genes)
remove <- unique(c(cc.genes$s.genes, cc.genes$g2m.genes, additional$S.derived,
                   additional$G2M.derived, stress_genes, hb_genes))

## check how many we will loose now:
length(intersect(VariableFeatures(srat), remove))
## [1] 132
## do the removal, for both the 'RNA' (i.e. LogNorm) and
## 'SCT' assays:
VariableFeatures(srat, assay = "RNA") <- setdiff(VariableFeatures(srat,
                                                                  assay = "RNA"), remove)

VariableFeatures(srat, assay = "SCT") <- setdiff(VariableFeatures(srat,
                                                                  assay = "SCT"), remove)

#Rerun dimred
srat <- RunPCA(srat, verbose = TRUE, npcs = 50)
srat <- RunUMAP(srat, dims = 1:CFG$ndims)

p_type <- DimPlot(srat, reduction = "umap", pt.size = 0.5)

p_phase <- DimPlot(srat, reduction = "umap", pt.size = 0.5, group.by = "Phase")

## use patchwork again to get a 2 x 2 plot
(p_type | p_phase)

## and also the continuous variables:
FeaturePlot(
  srat,
  pt.size = 0.5,
  feature = c("nCount_RNA", "percent_mito",
              "pct_hemo", "stress1"), order = TRUE)

#Clustering ----
srat <- FindNeighbors(srat, dims = 1:CFG$ndims)
## Computing nearest neighbor graph
## Computing SNN

#Check different resos
srat <- FindClusters(srat, resolution = seq(0.4, 2, 0.2), algorithm = 1)

columns <- grep("SCT_snn_res", names(srat@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the cluster ids)
## See changes in number of clusters and the cells at each cluster
for (col in columns) {
  cat("\n====\n", col, ":")
  print(table(srat@meta.data[[col]]))
}
##===
##  SCT_snn_res.0.4 :
##  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15 
##685 514 487 324 305 296 268 165 138 135  90  71  62  41  41  34 

##===
##  SCT_snn_res.0.8 :
##  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19 
##683 348 325 305 296 291 197 161 139 137 124 116 111  90  80  70  62  44  41  36 

##===
##  SCT_snn_res.1.2 :
##  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23 
##342 337 325 305 296 296 214 171 141 137 130 119 111 109  90  77  76  71  71  64  53  44  41  36 

## let's plot the different clusterings side by side to see
## the differeces:

p_clusbyreso <- list()  # set up the list to hold the intermediate plots

for (reso in c("0.4", "0.8", "1.2")) {
  
  clus.column <- paste0("SCT_snn_res.", reso)
  clusters <- sort(unique(srat@meta.data[[clus.column]]))
  
  ## switching color scheme to make it look different:
  clus.colors <- Polychrome::palette36.colors(length(clusters))
  ## clus.colors <-
  ## Polychrome::kelly.colors(length(clusters))
  names(clus.colors) <- as.character(clusters)
  
  p <- DimPlot(srat, pt.size = 0.5, group.by = clus.column, cols = clus.colors) +
    labs(title = reso)
  p_clusbyreso[[reso]] <- p
}

## use patchwork to plot them side-by-side:
p_phase | p_clusbyreso[[2]] | p_clusbyreso[[3]]

## set the identies to the clusters coming from the 0.8
## resolution:

clus.column <- paste0("SCT_snn_res.", 1.6)
Idents(srat) <- clus.column

clusters <- sort(unique(srat@meta.data[[clus.column]]))
## switching color scheme to make it look different:
cluster.colors <- Polychrome::palette36.colors(length(clusters))
names(cluster.colors) <- as.character(clusters)

## and show result with the other umaps to get your
## bearings:
p_clus <- DimPlot(srat, reduction = "umap", pt.size = 0.5, group.by = clus.column, cols = cluster.colors, label = TRUE)

(p_type | p_clus)
rm(p_clusbyreso)
gc()

## Using clustree ----
#srat <- FindClusters(srat, resolution = seq(0.4, 2, 0.2), algorithm = 1)
clustree::clustree(srat, prefix = "SCT_snn_res.", node_size = 2, edge_width = 0.5)

## Re-check whether previous resolution looks good
clus.column <- paste0("SCT_snn_res.", 1.6)
Idents(srat) <- clus.column
p_clus <- DimPlot(srat, reduction = "umap", pt.size = 0.5, group.by = clus.column, label = TRUE)

## scDblFinder ----
#This tool allow us to detect doublets from same genotype
sce <- scDblFinder(GetAssayData(srat, assay = "RNA", slot = "counts"), 
            iter = 10, 
            includePCs = 14,
            dims = 20)

srat$scDblFinder.score <- sce$scDblFinder.score
srat$scDblFinder.class <- sce$scDblFinder.class

p_sc.class <- DimPlot(srat, group.by = 'scDblFinder.class') +
  ggtitle('scDblrFinder class for SR')
p_sc.score <- FeaturePlot(srat, features = 'scDblFinder.score') +
  ggtitle('scDblrFinder score for SR')
p_sc.class | p_sc.score

#As with souporcell's annotation, keep only singlets
srat <- subset(srat, subset = scDblFinder.class == 'singlet')

#Enrichment analysis ----
##Differently expressed genes ----
## first step: find the differentially expressed genes this
## takes ~ 7min
#Finding CDs
CDs<-grep('^CD[0-9]{1,3}$', rownames(srat), value = TRUE)

de_genes <- FindAllMarkers(srat, assay = "RNA", only.pos = TRUE,
                           min.pct = 0.1, logfc.threshold = log(1.5))
## The resulting data.frame is too unwieldy to work with,
## the top 10 genes per cluster are enough:

de_genes %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
head(top10)

DoHeatmap(subset(srat, downsample = 50), features = top10$gene,
          group.colors = cluster.colors, assay = "RNA", slot = "scale.data") +
  theme(axis.text.y = element_text(size = 6))

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
enrichment <- compareCluster(clusters, fun = "enrichGO", OrgDb = "org.Hs.eg.db",
                             keyType = "ENSEMBL", ont = "BP",
                             minGSSize = 3, pAdjustMethod = "BH",
                             pvalueCutoff = 0.01, qvalueCutoff = 0.05,
                             readable = TRUE)

## so save it:
file <- paste0(CFG$output_dir, "/enrichment.rds")
saveRDS(file = file, enrichment)
#srat <- readRDS("~/analysis_illumina_LX279/enrichment.rds")

## show results:
clusterProfiler::dotplot(enrichment, showCategory = 4, 
title = "GO biological process", font = 5)

#Cell identification ----
## (down)load the reference data:
hpca <- HumanPrimaryCellAtlasData()

singler_hpca <- SingleR(
  test = GetAssayData(srat, assay = "RNA",
                      slot = "data"),
  ref = hpca,
  labels = hpca$label.main
)
## plot the heatmap
plotScoreHeatmap(singler_hpca, show.labels = TRUE, max.labels = 100,
                 show.pruned = FALSE, order.by = "clusters",
                 clusters = srat@meta.data[[clus.column]])

##Add new annotation to metadata
srat@meta.data$singler <- singler_hpca$labels

#To explore cell composition:
cbind(table(singler_hpca$labels))

#plot it:
typenames <- unique(srat@meta.data$singler)
singler.colors <- Polychrome::palette36.colors(length(typenames))
names(singler.colors) <- typenames

p_singler.sr <- DimPlot(srat, reduction = "umap", pt.size = 0.5,
                     group.by = "singler", cols = singler.colors)

p_singler.sr | p_clus

file <- paste0(CFG$output_dir, "/srat-cellsidentified.rds")
saveRDS(file = file, srat)
#srat <- readRDS(file = file)



#InferCNV ----
## overview of numbers per celltype:
table(srat@meta.data$singler)

#Simplify cell types by reannotating underrepresented cell types
cell_percentage_ilmn <- sum(table(srat$singler))

srat$singler <- ifelse(
  table(srat$singler)[srat$singler] / cell_percentage_ont < 0.01,
  'Other',
  srat$singler
)

table(srat$singler)

#Simplification by patient
srat@meta.data <- srat@meta.data %>%
  mutate(
    cells_by_patient = case_when(
      patient == '0' & singler == 'B_cell' ~ 'P0 B cell',
      patient == '0' & singler == 'Monocyte' ~ 'P0 Monocyte',
      patient == '0' & singler == 'NK_cell' ~ 'P0 NK cell',
      patient == '0' & singler == 'Other' ~ 'P0 Other',
      patient == '0' & singler == 'Pro-B_cell_CD34+' ~ 'P0 Pro-B cell CD34+',
      patient == '0' & singler == 'T_cells' ~ 'P0 T cells',
      patient == '1' & singler == 'B_cell' ~ 'P1 B Cell',
      patient == '1' & singler == 'Monocyte' ~ 'P1 Monocyte',
      patient == '1' & singler == 'NK_cell' ~ 'P1 NK cell',
      patient == '1' & singler == 'Other' ~ 'P1 Other',
      patient == '1' & singler == 'Pro-B_cell_CD34+' ~ 'P1 Pro-B cell CD34+',
      patient == '1' & singler == 'T_cells' ~ 'P1 T cells',
      TRUE ~ 'NA'
    )
  )

## which cells do we think we can trust as being healthy:
ref_types <- c(
  'P0 B cell',
  'P1 B Cell',
  'P0 Monocyte',
  'P1 Monocyte',
  'P0 T cells',
  'P1 T cells',
  'P0 NK cell',
  'P1 NK cell'
)

##And which ones are suspected to be tumor
maybe_tumor <- !(srat@meta.data$cells_by_patient %in% ref_types)
table(maybe_tumor)

## infercnv needs them as an external data file in a
## specific format, create that:
df <- data.frame(cell = colnames(srat), type = srat@meta.data$cells_by_patient)

celltypes_file <- paste0(CFG$output_dir, "/celltypes.txt")
write.table(df, file = celltypes_file, sep = "\t", quote = FALSE,
            na = "", row.names = FALSE, col.names = FALSE)
## infercnv also needs to gene coordinates, in the
## tab-delimited format genename chrNN START END.  We can
## use the `genecoords` object from SCutils for this.

## this data needs some cleaning. Use only the canonical
## chromosomes, sort them, get rid of duplicates and put
## 'chr' in front of the chromosome:
data(refdata_gex_GRCh38_2020_A)
wanted.chr <- c(as.character(1:22), "X", "Y")  # skip MT and KI2707 etc.
genecoords$chr <- gsub("chr", "", genecoords$chr)
geneorder <- genecoords[genecoords$chr %in% wanted.chr, ]
geneorder$numeric_chr <- as.numeric(geneorder$chr)
## Warning: NAs introduced by coercion
geneorder[geneorder$chr %in% "X", "numeric_chr"] <- 30
geneorder[geneorder$chr %in% "Y", "numeric_chr"] <- 40
geneorder <- geneorder[order(geneorder$numeric_chr, geneorder$start),]
geneorder <- geneorder[!duplicated(geneorder$gene_name), ]
geneorder <- with(geneorder, data.frame(gene = gene_name,
                                        chromosome = chr,
                                        start = start,
                                        end = end))
## @@@@ WRONG DIR, FIX: geneorder_file <-
## paste0(CFG$data_dir, '/geneorder.txt')
geneorder_file <- paste0(CFG$output_dir, "/geneorder.txt")
write.table(geneorder, file = geneorder_file, sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE, na = "")

raw_counts <- GetAssayData(srat, assay = "RNA", slot = "counts")

raw_counts <- as.data.frame(raw_counts)
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
outdir <- paste0(CFG$output_dir, "/infercnv_by_patient")
dir.create(outdir)
## running takes can take up to 10 min, best start this
## before the lecture:
infercnv_obj <- infercnv::run(infercnv_obj, out_dir = outdir,
                              cluster_by_groups = TRUE, denoise = TRUE, analysis_mode = "samples",
                              num_threads = 3, output_format = "png", window_length = 101,
                              save_final_rds = TRUE, plot_steps = FALSE, cutoff = 0.1,
                              sd_amplifier = 1.5, HMM = FALSE, write_expr_matrix = TRUE)
### @@@ this has changed in version 1.14, so the names of
### files are different again! FIX THIS
refexp <- as.matrix(read.delim(file = paste0(outdir, "/infercnv.references.txt"),
                               sep = " "))
#idk why but the "-" is changed to "." so we change it again
colnames(refexp) <- gsub(".", "-", colnames(refexp), fixed = TRUE)
obsexp <- as.matrix(read.delim(file = paste0(outdir, "/infercnv.observations.txt"),
                               sep = " "))
colnames(obsexp) <- gsub(".", "-", colnames(obsexp), fixed = TRUE)


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

# Create data frames for observed and reference scores
df_obs <- data.frame(score = obsscore, group = "tumor")
df_ref <- data.frame(score = refscore, group = "reference")

# Create density estimates for both sets of scores
density_obs <- density(obsscore)
density_ref <- density(refscore)

df_density_obs <- data.frame(x = density_obs$x, y = density_obs$y, group = "tumor")
df_density_ref <- data.frame(x = density_ref$x, y = density_ref$y, group = "reference")

tumorlines <- seq(from = 0, to = 1000, by = 100)
tumorvlines <-
  geom_vline(
    xintercept = tumorlines,
    col = "grey",
    linetype = 1,
    lwd = 0.1
  )

# Create the base plot
ggplot() +
  geom_rug(
    data = df_obs,
    aes(x = score),
    sides = "b",
    color = "black",
    size = 0.1
  ) +
  geom_rug(
    data = df_ref,
    aes(x = score),
    sides = "t",
    color = "red",
    size = 0.1
  ) +
  geom_line(data = df_density_obs, aes(x = x, y = y, color = group), size = 1) +
  geom_line(data = df_density_ref, aes(x = x, y = y, color = group), size = 1) +
  labs(
    title = sprintf("Aneuploidy score length=%d SDs=%.1f", len, SDs),
    x = sprintf("t=%.2f", ttest$statistic),
    y = "Density"
  ) +
  scale_color_manual(values = c("tumor" = "black", "reference" = "red")) +
  theme(legend.position = "topright") +
  guides(color = guide_legend(title = NULL)) +
  tumorvlines

cutoff <- 800  # based on the density plot

allscore <- c(obsscore, refscore)
allscore <- allscore[colnames(srat)]  # put in order of srat object.

srat <- AddMetaData(srat, allscore, col.name = "aneuploidy_score")

aneuploid <- ifelse(
  srat@meta.data$singler == "Monocyte",
  "no",
  ifelse(srat@meta.data$aneuploidy_score < cutoff, "no", "yes")
)

srat <- AddMetaData(srat, col.name = "is_aneuploid", aneuploid)

srat@meta.data <- srat@meta.data %>%
  mutate(
    reference = case_when(
      patient == '0' & is_aneuploid == 'yes' ~ 'p0.tumor',
      patient == '0' & is_aneuploid == 'no' ~ 'p0.healthy',
      patient == '1' & is_aneuploid == 'yes' ~ 'p1.tumor',
      patient == '1' & is_aneuploid == 'no' ~ 'p1.healthy',
      TRUE ~ 'NA'
    )
  )

## Let's show all we have. If it's too much information,
## leave some of it out.

p_pat <- DimPlot(srat, reduction = "umap", pt.size = 0.5,
                     group.by = "patient")

p_celltype <- DimPlot(srat, reduction = "umap", pt.size = 0.5,
                      group.by = "singler",cols = singler.colors)

p_aneup <- FeaturePlot(srat, reduction = "umap", pt.size = 0.5,
                       features = "aneuploidy_score",
                       order = TRUE)

p_isaneup <- DimPlot(srat, reduction = "umap", pt.size = 0.5,
                     group.by = "is_aneuploid")

p_count <- FeaturePlot(srat, reduction = "umap", pt.size = 0.5,
                          features = "nCount_RNA", cols = CFG$gradient.colors,
                          order = TRUE, label = TRUE) + ggtitle("Counts per cell")



( p_pat | p_celltype)/(p_aneup | p_isaneup)

p_count | p_singler

##Last saved srat ----
file<-paste0(CFG$output_dir, "/srat-aneuploidy.rds")
saveRDS(file=file,srat)
srat <- readRDS(file=file)

#InferCNV patient 1 (0 in Ilmn) ----
srat.ilmn.0 <- subset(srat, subset = patient == '0')

ref_types <- c('B_cell','T_cells','NK_cell')
maybe_tumor <- !(srat.ilmn.0@meta.data$singler %in% ref_types)
table(maybe_tumor)

## infercnv needs them as an external data file in a
## specific format, create that:
df <- data.frame(cell = colnames(srat.ilmn.0), type = srat.ilmn.0@meta.data$singler)

celltypes_file <- paste0(CFG$output_dir, "/celltypes_p0.txt")
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
geneorder <- geneorder[order(geneorder$numeric_chr, geneorder$start),]
geneorder <- geneorder[!duplicated(geneorder$gene_name), ]
geneorder <- with(geneorder, data.frame(gene = gene_name,
                                        chromosome = paste0("chr",chr),
                                        start = start, end = end))
## @@@@ WRONG DIR, FIX: geneorder_file <-
## paste0(CFG$data_dir, '/geneorder.txt')
geneorder_file <- paste0(CFG$output_dir, "/geneorder_p0.txt")
write.table(geneorder, file = geneorder_file, sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE, na = "")

raw_counts <- GetAssayData(srat.ilmn.0, assay = "RNA", slot = "counts")

raw_counts <- as.data.frame(raw_counts)
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
outdir <- paste0(CFG$output_dir, "/infercnv_p0")
dir.create(outdir)
## running takes can take up to 10 min, best start this
## before the lecture:
infercnv_obj <- infercnv::run(infercnv_obj, out_dir = outdir,
                              cluster_by_groups = TRUE, denoise = TRUE, analysis_mode = "samples",
                              num_threads = 3, output_format = "png", window_length = 101,
                              save_final_rds = TRUE, plot_steps = FALSE, cutoff = 0.1,
                              sd_amplifier = 1.5, HMM = FALSE, write_expr_matrix = TRUE)
### @@@ this has changed in version 1.14, so the names of
### files are different again! FIX THIS
refexp <- as.matrix(read.delim(file = paste0(outdir, "/infercnv.references.txt"),
                               sep = " "))
#idk why but the "-" is changed to "." so we change it again
colnames(refexp) <- gsub(".", "-", colnames(refexp), fixed = TRUE)
obsexp <- as.matrix(read.delim(file = paste0(outdir, "/infercnv.observations.txt"),
                               sep = " "))
colnames(obsexp) <- gsub(".", "-", colnames(obsexp), fixed = TRUE)


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

# Create data frames for observed and reference scores
df_obs <- data.frame(score = obsscore, group = "tumor")
df_ref <- data.frame(score = refscore, group = "reference")

# Create density estimates for both sets of scores
density_obs <- density(obsscore)
density_ref <- density(refscore)

df_density_obs <- data.frame(x = density_obs$x, y = density_obs$y, group = "tumor")
df_density_ref <- data.frame(x = density_ref$x, y = density_ref$y, group = "reference")

tumorlines <- seq(from = 0, to = 1000, by = 100)
tumorvlines <-
  geom_vline(
    xintercept = tumorlines,
    col = "grey",
    linetype = 1,
    lwd = 0.1
  )

# Create the base plot
ggplot() +
  geom_rug(
    data = df_obs,
    aes(x = score),
    sides = "b",
    color = "black",
    size = 0.1
  ) +
  geom_rug(
    data = df_ref,
    aes(x = score),
    sides = "t",
    color = "red",
    size = 0.1
  ) +
  geom_line(data = df_density_obs, aes(x = x, y = y, color = group), size = 1) +
  geom_line(data = df_density_ref, aes(x = x, y = y, color = group), size = 1) +
  labs(
    title = sprintf("Aneuploidy score length=%d SDs=%.1f", len, SDs),
    x = sprintf("t=%.2f", ttest$statistic),
    y = "Density"
  ) +
  scale_color_manual(values = c("tumor" = "black", "reference" = "red")) +
  theme(legend.position = "topright") +
  guides(color = guide_legend(title = NULL)) +
  tumorvlines

cutoff <- 750  # based on the density plot

allscore <- c(obsscore, refscore)
allscore <- allscore[colnames(srat.ilmn.0)]  # put in order of srat.ilmn.0 object.

srat.ilmn.0 <- AddMetaData(srat.ilmn.0, allscore, col.name = "aneuploidy_score")

aneuploid <- ifelse(
  srat.ilmn.0@meta.data$singler == "Monocyte",
  "no",
  ifelse(srat.ilmn.0@meta.data$aneuploidy_score < cutoff, "no", "yes")
)

srat.ilmn.0 <- AddMetaData(srat.ilmn.0, col.name = "is_aneuploid", aneuploid)

## Let's show all we have. If it's too much information,
## leave some of it out.

p_pat <- DimPlot(srat.ilmn.0, reduction = "umap", pt.size = 0.5,
                 group.by = "patient")

p_celltype <- DimPlot(srat.ilmn.0, reduction = "umap", pt.size = 0.5,
                      group.by = "singler",cols = singler.colors)

p_aneup <- FeaturePlot(srat.ilmn.0, reduction = "umap", pt.size = 0.5,
                       features = "aneuploidy_score",
                       order = TRUE)

p_isaneup <- DimPlot(srat.ilmn.0, reduction = "umap", pt.size = 0.5,
                     group.by = "is_aneuploid")

p_count <- FeaturePlot(srat.ilmn.0, reduction = "umap", pt.size = 0.5,
                       features = "nCount_RNA", cols = CFG$gradient.colors,
                       order = TRUE, label = TRUE) + ggtitle("Counts per cell")



( p_pat | p_celltype)/(p_aneup | p_isaneup)

p_count | p_singler



#InferCNV patient 1 (0 in Ilmn) ----
srat.ilmn.1 <- subset(srat, subset = patient == '1')

ref_types <- c('B_cell','T_cells','NK_cell')
maybe_tumor <- !(srat.ilmn.1@meta.data$singler %in% ref_types)
table(maybe_tumor)

## infercnv needs them as an external data file in a
## specific format, create that:
df <- data.frame(cell = colnames(srat.ilmn.1), type = srat.ilmn.1@meta.data$singler)

celltypes_file <- paste0(CFG$output_dir, "/celltypes_p1.txt")
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
geneorder <- geneorder[order(geneorder$numeric_chr, geneorder$start),]
geneorder <- geneorder[!duplicated(geneorder$gene_name), ]
geneorder <- with(geneorder, data.frame(gene = gene_name,
                                        chromosome = paste0("chr",chr),
                                        start = start, end = end))
## @@@@ WRONG DIR, FIX: geneorder_file <-
## paste0(CFG$data_dir, '/geneorder.txt')
geneorder_file <- paste0(CFG$output_dir, "/geneorder_p1.txt")
write.table(geneorder, file = geneorder_file, sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE, na = "")

raw_counts <- GetAssayData(srat.ilmn.1, assay = "RNA", slot = "counts")

raw_counts <- as.data.frame(raw_counts)
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
outdir <- paste0(CFG$output_dir, "/infercnv_p1")
dir.create(outdir)
## running takes can take up to 10 min, best start this
## before the lecture:
infercnv_obj <- infercnv::run(infercnv_obj, out_dir = outdir,
                              cluster_by_groups = TRUE, denoise = TRUE, analysis_mode = "samples",
                              num_threads = 3, output_format = "png", window_length = 101,
                              save_final_rds = TRUE, plot_steps = FALSE, cutoff = 0.1,
                              sd_amplifier = 1.5, HMM = FALSE, write_expr_matrix = TRUE)
### @@@ this has changed in version 1.14, so the names of
### files are different again! FIX THIS
refexp <- as.matrix(read.delim(file = paste0(outdir, "/infercnv.references.txt"),
                               sep = " "))
#idk why but the "-" is changed to "." so we change it again
colnames(refexp) <- gsub(".", "-", colnames(refexp), fixed = TRUE)
obsexp <- as.matrix(read.delim(file = paste0(outdir, "/infercnv.observations.txt"),
                               sep = " "))
colnames(obsexp) <- gsub(".", "-", colnames(obsexp), fixed = TRUE)


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

# Create data frames for observed and reference scores
df_obs <- data.frame(score = obsscore, group = "tumor")
df_ref <- data.frame(score = refscore, group = "reference")

# Create density estimates for both sets of scores
density_obs <- density(obsscore)
density_ref <- density(refscore)

df_density_obs <- data.frame(x = density_obs$x, y = density_obs$y, group = "tumor")
df_density_ref <- data.frame(x = density_ref$x, y = density_ref$y, group = "reference")

tumorlines <- seq(from = 0, to = 1000, by = 100)
tumorvlines <-
  geom_vline(
    xintercept = tumorlines,
    col = "grey",
    linetype = 1,
    lwd = 0.1
  )

# Create the base plot
ggplot() +
  geom_rug(
    data = df_obs,
    aes(x = score),
    sides = "b",
    color = "black",
    size = 0.1
  ) +
  geom_rug(
    data = df_ref,
    aes(x = score),
    sides = "t",
    color = "red",
    size = 0.1
  ) +
  geom_line(data = df_density_obs, aes(x = x, y = y, color = group), size = 1) +
  geom_line(data = df_density_ref, aes(x = x, y = y, color = group), size = 1) +
  labs(
    title = sprintf("Aneuploidy score length=%d SDs=%.1f", len, SDs),
    x = sprintf("t=%.2f", ttest$statistic),
    y = "Density"
  ) +
  scale_color_manual(values = c("tumor" = "black", "reference" = "red")) +
  theme(legend.position = "topright") +
  guides(color = guide_legend(title = NULL)) +
  tumorvlines

cutoff <- 750  # based on the density plot

allscore <- c(obsscore, refscore)
allscore <- allscore[colnames(srat.ilmn.1)]  # put in order of srat.ilmn.1 object.

srat.ilmn.1 <- AddMetaData(srat.ilmn.1, allscore, col.name = "aneuploidy_score")

aneuploid <- ifelse(
  srat.ilmn.1@meta.data$singler == "Monocyte",
  "no",
  ifelse(srat.ilmn.1@meta.data$aneuploidy_score < cutoff, "no", "yes")
)

srat.ilmn.1 <- AddMetaData(srat.ilmn.1, col.name = "is_aneuploid", aneuploid)

srat.ilmn.1@meta.data <- srat.ilmn.1@meta.data %>%
  mutate(
    reference = case_when(
      patient == '0' & is_aneuploid == 'yes' ~ 'p0.tumor',
      patient == '0' & is_aneuploid == 'no' ~ 'p0.healthy',
      patient == '1' & is_aneuploid == 'yes' ~ 'p1.tumor',
      patient == '1' & is_aneuploid == 'no' ~ 'p1.healthy',
      TRUE ~ 'NA'
    )
  )

## Let's show all we have. If it's too much information,
## leave some of it out.

p_pat <- DimPlot(srat.ilmn.1, reduction = "umap", pt.size = 0.5,
                 group.by = "patient")

p_celltype <- DimPlot(srat.ilmn.1, reduction = "umap", pt.size = 0.5,
                      group.by = "singler",cols = singler.colors)

p_aneup <- FeaturePlot(srat.ilmn.1, reduction = "umap", pt.size = 0.5,
                       features = "aneuploidy_score",
                       order = TRUE)

p_isaneup <- DimPlot(srat.ilmn.1, reduction = "umap", pt.size = 0.5,
                     group.by = "is_aneuploid")

p_count <- FeaturePlot(srat.ilmn.1, reduction = "umap", pt.size = 0.5,
                       features = "nCount_RNA", cols = CFG$gradient.colors,
                       order = TRUE, label = TRUE) + ggtitle("Counts per cell")



( p_pat | p_celltype)/(p_aneup | p_isaneup)

p_count | p_singler



#Saving bcs ----
#For illumina we have to add 'CB:Z:' in front of everything just in case
XX_tumor_bc <- rownames(srat@meta.data %>% filter(patient == 0) %>% filter(is_aneuploid == 'yes'))
XX_tumor_bc<-gsub('LX279_','CB:Z:',XX_tumor_bc)
XX_tumor_bc<-gsub('-1','', XX_tumor_bc)
XX_healthy_bc <- rownames(srat@meta.data %>% filter(patient == 0) %>% filter(is_aneuploid == 'no'))
XX_healthy_bc<-gsub('LX279_','CB:Z:',XX_healthy_bc)
XX_healthy_bc<-gsub('-1','', XX_healthy_bc)

XY_tumor_bc <- rownames(srat@meta.data %>% filter(patient == 1)%>% filter(is_aneuploid == 'yes'))
XY_tumor_bc<-gsub('LX279_','CB:Z:',XY_tumor_bc)
XY_tumor_bc<-gsub('-1','', XY_tumor_bc)
XY_healthy_bc <- rownames(srat@meta.data %>% filter(patient == 1)%>% filter(is_aneuploid == 'no'))
XY_healthy_bc<-gsub('LX279_','CB:Z:',XY_healthy_bc)
XY_healthy_bc<-gsub('-1','', XY_healthy_bc)

#bcs <- list(XX_tumor_bc, XX_healthy_bc, XY_tumor_bc, XY_healthy_bc)
#for (x in list(XX_tumor_bc, XX_healthy_bc, XY_tumor_bc, XY_healthy_bc)){
#  file_con <- file(paste0(x,'.txt'), "w")
#  writeLines(x, file_con)
#  close(file_con)
#}

file_con <- file('LR_XX_tumor_bc.txt', "w")
writeLines(XX_tumor_bc, file_con)
close(file_con)

file_con <- file('LR_XX_healthy_bc.txt', "w")
writeLines(XX_healthy_bc, file_con)
close(file_con)


file_con <- file('LR_XY_tumor_bc.txt', "w")
writeLines(XY_tumor_bc, file_con)
close(file_con)

file_con <- file('LR_XY_healthy_bc.txt', "w")
writeLines(XY_healthy_bc, file_con)
close(file_con)

