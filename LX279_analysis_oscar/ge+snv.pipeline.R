.libPaths(
  c(
    "/hpc/pmc_holstege/rstudio/oscar_clonetracer/R/x86_64-pc-linux-gnu-library/4.4",
    "/usr/local/lib/R/site-library",
    "/usr/local/lib/R/library"
  )
)

library(SCutils)
library(Matrix)  # needed for working with sparse matrices
library(Seurat)
library(ggplot2)
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
library(clusterProfiler)
library(tidyverse)
library(enrichR)
library(reshape2)
library(SeuratObject)

CFG <- list()  # global

CFG$data_dir <- "~/nanopore_souporcell"  # where to read data from
CFG$output_dir <- "~/nanopore_souporcell"  # where to write data to
CFG$project <- "library name or brief annotation of lib"

# explained later:
CFG$ndims <- 20
CFG$random_seed <- 2033
CFG$reso <- 1.6

# Cutoffs used during quality control for filtering genes
# (explained later)
CFG$min_txpts <- 1000
CFG$max_txpts <- 23000
CFG$max_pctmito <- 10

# Graphical parameters
CFG$gradient.colors <- viridis::viridis(101, direction = -1)

#For SNVs
CFG$QUAL <- 30 #Higher than 30
CFG$ndims.snvs <- 25
CFG$reso.snv <- 1.8
CFG$reso.snv.p0 <- 1.4
CFG$reso.snv.p1 <- 1.8
CFG$ndims.p0 <- 14
CFG$ndims.p1 <- 25
CFG$ndims.tumor <- 15
CFG$reso.ge.p0 <-
CFG$reso.ge.p1 <- 1.4

#Preparation for souporcell ----
#First 'cycle'/step/level. Matrix including all cells to split by patient.
#Run souporcell with snvs coming from scNano
sop.input <- './sop_input/'

alt.mtx <- read.table(paste0(sop.input, 'alt.mtx'), header = TRUE)
ref.mtx <- read.table(paste0(sop.input, 'ref.mtx'), header = TRUE)

vaf.per.cell <- read.table(file = '~/vaf.per.cell.tsv', header = TRUE)
cutoff <- list()
cutoff$QUAL <- 30 #Higher than 30
cutoff$sds <- 0.3 #SNVs selected must have more than n% SD

#String manipulation
vaf.per.cell$ID <- NULL
vaf.per.cell <- vaf.per.cell %>% unite(col = 'SNV', CHROM:ALT, sep = ':')
rownames(vaf.per.cell) <- vaf.per.cell$SNV
vaf.per.cell$SNV <- NULL

vaf.per.cell <- vaf.per.cell[vaf.per.cell$QUAL > cutoff$QUAL, ]
vaf.per.cell$QUAL <- NULL
vaf.per.cell <- as.matrix(vaf.per.cell)

#Split by nuc and mito SNVs
mitochondrial.snvs <- grepl("^chrM", rownames(vaf.per.cell))
mito.matrix <- vaf.per.cell[mitochondrial.snvs, ]
nuc.matrix <- vaf.per.cell[!mitochondrial.snvs, ]

#Nuclear SNVs
nuc.unfiltered <- rowMeans(nuc.matrix, na.rm = TRUE)
nuc.unfiltered <- as.data.frame(nuc.unfiltered)
colnames(nuc.unfiltered) <- 'mean'
nuc.unfiltered$sd <- rowSds(nuc.matrix, na.rm = TRUE)
(sum(is.na(nuc.matrix)) / (ncol(nuc.matrix) * nrow(nuc.matrix))) * 100
#[1] 80.06532 percentage of NAs in our unprocessed data
nuc.selected.snvs <- nuc.matrix[nuc.unfiltered$sd > cutoff$sds, ]
table(nuc.selected.snvs)

#Alt
alt.mtx$ID <- NULL
alt.mtx <- alt.mtx %>% unite(col = 'SNV', CHROM:ALT, sep = ':')
rownames(alt.mtx) <- alt.mtx$SNV
alt.mtx$SNV <- NULL
alt.mtx$QUAL <- NULL
nuc.alt.mtx <- alt.mtx[rownames(alt.mtx) %in% rownames(nuc.selected.snvs), ]

#Need barcodes for souporcell
barcodes <- colnames(alt.mtx)

rownames(nuc.alt.mtx) <- NULL
colnames(nuc.alt.mtx) <- NULL

#Ref
ref.mtx$ID <- NULL
ref.mtx <- ref.mtx %>% unite(col = 'SNV', CHROM:ALT, sep = ':')
rownames(ref.mtx) <- ref.mtx$SNV
ref.mtx$SNV <- NULL
ref.mtx$QUAL <- NULL
nuc.ref.mtx <- ref.mtx[rownames(ref.mtx) %in% rownames(nuc.selected.snvs), ]
View(nuc.ref.mtx)

#Saving matrices
#Little trick is needed here, and then change manually the -1
nuc.ref.mtx <- as.matrix(nuc.ref.mtx)
nuc.alt.mtx <- as.matrix(nuc.alt.mtx)

index.nuc.alt <- nuc.alt.mtx != 0
index.nuc.ref <- nuc.ref.mtx != 0

nuc.ref.mtx[index.nuc.alt] <- -1
nuc.alt.mtx[index.nuc.ref] <- -1

nuc.sparse.alt.mtx <- as(as.matrix(nuc.alt.mtx), "sparseMatrix")
nuc.sparse.ref.mtx <- as(as.matrix(nuc.ref.mtx), "sparseMatrix")

writeMM(nuc.sparse.alt.mtx,
        file = paste0(sop.input, 'input_filtered/market.alt.mtx'))
writeMM(nuc.sparse.ref.mtx,
        file = paste0(sop.input, 'input_filtered/market.ref.mtx'))
writeLines(barcodes, paste0(sop.input, "input_filtered/barcodes.txt"))

#To change matrix, see '~/sop_input/changing.matrix.sh'

# Data loading----
lib <- "LX279"  # name of the library, also needed later on.
dir <- paste0(CFG$data_dir, "/", lib)  # what is actually on disk
counts <- read.delim(
  paste0(
    CFG$data_dir,
    "/LX279_sup_complete.gene_expression.counts.tsv"
  ),
  row.names = 1
) #we do processing ourselves
counts <- as.matrix(counts)
counts <- as(counts, "sparseMatrix")
#is of type 'dgCMatrix', a sparse matrix
cnts_per_cell <- Matrix::colSums(counts)  # note use of Matrix::colSums

summary(cnts_per_cell)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#393    3787    5384    7143    9143   36228
#Quite low as compared to Illumina,

df <- data.frame(counts = cnts_per_cell)

ggplot(df, aes(x = counts)) + geom_density() + geom_rug() +
  scale_x_continuous(trans = "log2")

## Changing cell labels ----
cellnames <- colnames(counts)
colnames(counts) <- paste0(lib, "_", cellnames)


## Metadata ----
# get the mitochondrial genes based on their name:
mitos <- grep("^mt-",
              rownames(counts),
              value = TRUE,
              ignore.case = TRUE)

## show them:
mitos
##  [1] "MT-ND1"  "MT-ND2"  "MT-CO1"  "MT-CO2"  "MT-ATP8"
#"MT-ATP6" "MT-CO3"  "MT-ND3"  "MT-ND4L" "MT-ND4"  "MT-ND5"  "MT-ND6"  "MT-CYB"

percent_mito <- 100 * Matrix::colSums(counts[mitos, ]) / cnts_per_cell
# note the use of Matrix::colSums again

## get the non-mito's for easy selection of 'the rest of
## the genes'
nuclear <- setdiff(rownames(counts), mitos)

log2_transcript_counts <- log2(1 + Matrix::colSums(counts[nuclear, ]))

log2_feature_counts <- log2(1 + Matrix::colSums(counts[nuclear, ] > 0))

## set up a data.frame:
meta <- data.frame(
  percent_mito = percent_mito,
  log2_counts = log2_transcript_counts,
  log2_features = log2_feature_counts,
  lib = rep(lib, ncol(counts))
)

## let's show how features ( = genes ) relate to number of
## transcripts
ggplot(meta, aes(x = log2_counts, y = log2_features)) + geom_point(color = "blue")

## Creating Seurat object ----
counts <- counts[nuclear, ]

srat <- Seurat::CreateSeuratObject(counts = counts,
                                   project = CFG$project,
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
#24740 features across 3821 samples within 1 assay
#Active assay: RNA (24740 features, 0 variable features)

head(srat@meta.data)

##Removing doublets ----
#Barcode name is without the library name, so we have to add it
#Souporcell with ONT, requires vcf with scNano
sop_nano <- read.delim("~/sop_input/nuc/clusters.tsv",
                       header = TRUE,
                       row.names = 1)
rownames(sop_nano) <- paste0(lib, "_", rownames(sop_nano))
rownames(sop_nano) <- gsub("-1", "", rownames(sop_nano), fixed = TRUE)

goodBc <- rownames(sop_nano)[rownames(sop_nano) %in% colnames(srat)]

srat@meta.data[goodBc, 'patient_ont'] <- sop_nano[goodBc, 'assignment']
srat@meta.data[goodBc, 'status'] <- sop_nano[goodBc, 'status']

srat <- subset(srat, subset = status == "singlet")

##Souporcell with Ilmn (if available)
sop_ilmn <- read.delim("~/clusters_k2_ilmn.tsv",
                       header = TRUE,
                       row.names = 1)
rownames(sop_ilmn) <- paste0("LX279_", rownames(sop_ilmn))
rownames(sop_ilmn) <- gsub('-1', '', rownames(sop_ilmn))

goodBc <- rownames(sop_ilmn)[rownames(sop_ilmn) %in% colnames(srat)]

srat@meta.data[goodBc, 'patient_ilmn'] <- sop_ilmn[goodBc, 'assignment']
srat$patient_ilmn <- recode(srat$patient_ilmn, `0` = "1", `1` = "0")

#if available, check overlap
table(srat$patient_ont, srat$patient_ilmn)

#Filtering cells ----
##Transcript counts ----
v <- VlnPlot(srat, "nCount_RNA", group.by = "lib")
lines <- seq(from = 0, to = 36000, by = 2000)
hlines <- geom_hline(
  yintercept = lines,
  col = "grey",
  linetype = 1,
  lwd = 0.1
)

v + hlines
v + scale_y_continuous(trans = "log2") + hlines

## create a scatter graph
f_lin <-
  FeatureScatter(srat,
                 feature1 = "nCount_RNA",
                 feature2 = "nFeature_RNA",
                 pt.size = 0.5) + geom_hline(yintercept = CFG$max_pctmito, linetype = 2) +
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
#[1] 24740  3599
dim(subset(srat, nCount_RNA < CFG$min_txpts))
## [1] 24740   99
dim(subset(srat, nCount_RNA > CFG$max_txpts))
## [1] 24740   50
dim(subset(srat, percent_mito > CFG$max_pctmito))
## [1] 24740   316

## lastly, subset the seurat object and save if needed
srat <- subset(
  srat,
  subset = nCount_RNA >= CFG$min_txpts & nCount_RNA <=
    CFG$max_txpts & percent_mito <= CFG$max_pctmito
)

dim(srat)
##[1] 24740  3557

#Hemoglobin removal

data(refdata_cellranger_GRCh38_3.0.0)

hb_genes <- genelists$hemo
## alternatively: hb_genes <- lookup.hemogenes()

hb_genes <- intersect(hb_genes, rownames(srat))

hb_counts <- Matrix::colSums(srat@assays$RNA$counts[hb_genes, ])

srat <- AddMetaData(srat, col.name = "log2hb_genes", metadata = log2(1 + hb_counts))

srat <- AddMetaData(
  srat,
  col.name = "pct_hemo",
  metadata = 100 *
    hb_counts / srat@meta.data$nCount_RNA
)

## now show the previous plot along side a plot of the hemo
## content:
df <- subset(srat@meta.data)

p_ngenes <- ggplot(df, aes(x = log2_counts, y = log2_features, color = pct_hemo)) + geom_point(size = 1, alpha = 1 /
                                                                                                 2) +
  scale_color_gradient(
    low = "blue",
    high = "red",
    limits = c(0, 5),
    oob = scales::squish
  )

p_pcthemo <- ggplot(df, aes(x = log2_counts, y = pct_hemo, color = log2_features)) +
  geom_point(size = 1, alpha = 1 / 2) +
  scale_color_gradient(
    low = "blue",
    high = "red",
    limits = c(2, 16),
    oob = scales::squish
  )

p_ngenes | p_pcthemo


p <- ggplot(df, aes(x = nCount_RNA, y = nFeature_RNA)) +
  xlab("RNA Counts") +
  ylab("Features") +
  geom_point(size = 0.5, color = "red")
p

## and/or logarithmic:

p + scale_x_continuous(trans = "log2") + scale_y_continuous(trans = "log2")
#add anotation
ggplot(srat@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  xlab("RNA Counts") +
  ylab("Features") +
  geom_point(size = 0.5, color = "red")

CFG$max_pct_hemo <- 5

dim(srat)
## [1] 24740  3557
dim(subset(srat, pct_hemo > CFG$max_pct_hemo))
## [1] 24740    28
srat <- subset(srat, pct_hemo <= CFG$max_pct_hemo)

p_withouthemo <- ggplot(srat@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  xlab("RNA Counts") +
  ylab("Features") +
  geom_point(size = 0.5, color = "red")

p | p_withouthemo

## Saved srat ----
file <- paste0(CFG$output_dir, "/srat-cellsfiltered.rds")
#saveRDS(file = file, srat)
#srat <-readRDS(file)

#Normalization and Dimensional Reduction ----
srat <- NormalizeData(srat, normalization.method = "LogNormalize")

srat <- ScaleData(srat, features = rownames(srat), verbose = FALSE)

srat <- FindVariableFeatures(srat)


srat <- suppressMessages(
  SCTransform(
    srat,
    vars.to.regres = NULL,
    vst.flavor = 'v2',
    ## version that solved some problems
    verbose = FALSE,
    variable.features.n = 3000
  )
)


##Saved srat ----

file <- paste0(CFG$output_dir, "/srat-normalized.rds")
#saveRDS(file=file,srat)
#srat <- readRDS(file = file)
##T o read the object back in, do
#srat <- readRDS("~/analysis_illumina_LX279/srat-normalized.rds")

## select the most highly variable genes `VariableFeatures`
## returns them sorted by decreasing variance.)
top10 <- head(VariableFeatures(srat), 10)

p.varfeat <- VariableFeaturePlot(srat, selection.method = 'sct')

p.varfeat <- LabelPoints(plot = p.varfeat,
                         points = top10,
                         repel = TRUE)

p.varfeat

DefaultAssay(srat)  # just check that it has 'SCT'
## [1] "SCT"
srat <- RunPCA(srat, npcs = 50)

DimPlot(srat,
        reduction = "pca",
        group.by = 'patient_ont',
        pt.size = 0.5)

elbow <- ElbowPlot(srat, reduction = "pca", ndims = 50)
elbow #geom_vline(xintercept = c(seq(from = 0, to = 50, by = 1)))

DimHeatmap(srat, dims = (1:18), cells = 100)

DimHeatmap(srat, dims = 18 + (1:12), cells = 100)

data(cc.genes)  # load the cell cycle genes (also available from SCutils' genelists)

srat <- CellCycleScoring(
  object = srat,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  assay = "RNA"
)

#See state of our cells before deconfounding
RidgePlot(
  srat,
  features = c("PCNA", "TOP2A", "MCM6", "MKI67"),
  ncol = 2,
  group.by = 'Phase'
)

## Warning: The following features are not present in the object:
#MLF1IP, not searching for symbol synonyms

stress_genes <- genelists$stress  # or lookup.stressgenes()
srat <- AddModuleScore(
  srat,
  features = list(stress = stress_genes),
  name = "stress",
  assay = "RNA"
)
## Warning: The following features are not present in the object:
#COMP, PRKN, AGR2, ERN2, FGF21, HSPB8, PACRG, BHLHA15, MIR199A1,
#MMP24-AS1-EDEM2, not searching for symbol synonyms

srat <- RunUMAP(srat, dims = 1:CFG$ndims)

p_phase <- DimPlot(srat,
                   reduction = "umap",
                   pt.size = 0.5,
                   group.by = "Phase")
#check if there are cells being clustered based on their cell cycle phase
p_phase

norm1 <- FeaturePlot(
  srat,
  pt.size = 0.5,
  feature = c("nCount_RNA", "percent_mito", "nFeature_RNA", "stress1"),
  order = TRUE
)
norm1

# Deconfounding ----
## In this part we remove features, not cells
Scor <- metadataCorrelations(srat, "S.Score")

G2Mcor <- metadataCorrelations(srat, "G2M.Score")

additional <- derivedCellcycleGenes2(
  Scor = Scor,
  Sgenes = genelists$s.genes,
  G2Mcor = G2Mcor,
  G2Mgenes = genelists$g2m.genes
)
additional
additional$plot[[1]] | additional$plot[[2]] | additional$plot[[4]]

## also get the other deconfounders:
stress_genes <- genelists$stress #OR: SCutils::lookup.stressgenes()
hb_genes <- genelists$hemo  # as before
fem_genes <- genelists$female
male_genes <- genelists$male

## combine these sets (also including stress_genes and
## hemoglobin genes)
remove <- unique(
  c(
    cc.genes$s.genes,
    cc.genes$g2m.genes,
    additional$S.derived,
    additional$G2M.derived,
    stress_genes,
    hb_genes,
    fem_genes,
    male_genes
  )
)

## check how many we will loose now:
length(intersect(VariableFeatures(srat), remove))
## [1] 126
## do the removal, for both the 'RNA' (i.e. LogNorm) and
## 'SCT' assays:
VariableFeatures(srat, assay = "RNA") <-
  setdiff(VariableFeatures(srat, assay = "RNA"), remove)

VariableFeatures(srat, assay = "SCT") <-
  setdiff(VariableFeatures(srat, assay = "SCT"), remove)
srat <- RunPCA(srat, verbose = TRUE, npcs = 50)
srat <- RunUMAP(srat, dims = 1:CFG$ndims)

p_phase <- DimPlot(srat,
                   reduction = "umap",
                   pt.size = 0.5,
                   group.by = "Phase")

p_phase

## and also the continuous variables:
FeaturePlot(
  srat,
  pt.size = 0.5,
  feature = c("nCount_RNA", "percent_mito", "nFeature_RNA", "stress1"),
  order = TRUE
)

#If lib has male and female, extra info that might be useful
FeaturePlot(srat, features = fem_genes) |
  DimPlot(srat, group.by = 'patient_ont')

#Clustering ----
srat <- FindNeighbors(srat, dims = 1:CFG$ndims)
## Computing nearest neighbor graph
## Computing SNN

## Using clustree ----
srat <- FindClusters(srat, resolution = seq(0.4, 2, 0.2), algorithm = 1)
clustree::clustree(srat,
                   prefix = "SCT_snn_res.",
                   node_size = 2,
                   edge_width = 0.5)
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
##
## Number of nodes: 1296
## Number of edges: 37806
##
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9242
## Number of communities: 8
## Elapsed time: 0 seconds

## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
##
## Number of nodes: 1296
## Number of edges: 37806
##
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8690
## Number of communities: 12
## Elapsed time: 0 seconds

columns <- grep("SCT_snn_res", names(srat@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the
## cluster ids)
for (col in columns) {
  cat("\n====\n", col, ":")
  print(table(srat@meta.data[[col]]))
}

#===
#  SCT_snn_res.0.4 :
#  0   1   2   3   4   5   6   7   8   9  10  11
#688 527 473 470 246 232 198 125  85  79  38  32

#===
#  SCT_snn_res.0.6 :
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
#577 502 304 282 246 213 198 166 135 128 109  85  81  61  38  35  33

#===
#  SCT_snn_res.0.8 :
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17
#576 502 304 278 246 198 188 166 135 132 110  85  81  61  38  35  33  25

#===
#  SCT_snn_res.1 :
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
#551 499 304 271 266 208 166 137 127 110 102  94  85  81  61  38  35  33  25

#===
#  SCT_snn_res.1.2 :
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
#498 336 303 262 208 194 167 160 137 131 123 119 102  94  85  81  61  38  35  33  26

#===
#  SCT_snn_res.1.4 :
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
#497 333 298 262 209 194 172 160 135 131 126 119 104  94  85  81  61  38  35  33  26

#===
#  SCT_snn_res.1.6 :
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23
#352 292 232 207 207 166 160 156 148 148 140 139 133 127 100  94  85  81  61  38  35  33  33  26

#===
#  SCT_snn_res.1.8 :
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23
#331 309 232 206 195 189 166 160 156 148 135 127 126 120 104  94  85  81  61  38  35  33  33  29

#===
#  SCT_snn_res.2 :
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
#313 234 213 197 192 166 161 156 148 142 137 134 126 124 101  94  88  85  81  76  61  38  35  33  33  25

## let's plot the different clusterings side by side to see
## the differences:

p_clusbyreso <- list()  # set up the list to hold the intermediate plots

for (reso in seq(0.4, 2, 0.2)) {
  clus.column <- paste0("SCT_snn_res.", reso)
  clusters <- sort(unique(srat@meta.data[[clus.column]]))
  
  ## switching color scheme to make it look different:
  clus.colors <- Polychrome::alphabet.colors(length(clusters))
  ## clus.colors <-
  ## Polychrome::kelly.colors(length(clusters))
  names(clus.colors) <- as.character(clusters)
  
  p <- DimPlot(srat,
               pt.size = 0.5,
               group.by = clus.column,
               cols = clus.colors) +
    labs(title = reso)
  p_clusbyreso[[as.character(reso)]] <- p
}

## use patchwork to plot them side-by-side:
p_phase | p_clusbyreso[[8]] | p_clusbyreso[[9]]

## set the identities to the clusters coming from the 0.8
## resolution:

#Choose reso based on num of clusters and
clus.column <- paste0("SCT_snn_res.", CFG$reso)
Idents(srat) <- clus.column

clusters <- sort(unique(srat@meta.data[[clus.column]]))
## switching color scheme to make it look different:
cluster.colors <- ggsci::pal_igv(palette = "default")(length(clusters))
names(cluster.colors) <- as.character(clusters)

## and show result with the other umaps to get your
## bearings:
p_clus <- DimPlot(
  srat,
  reduction = "umap",
  pt.size = 0.5,
  group.by = clus.column,
  cols = cluster.colors,
  label = TRUE
) + NoLegend()

p_clus
rm(p_clusbyreso)
gc()

## scDblFinder ----
sce.lr <- scDblFinder(
  GetAssayData(srat, assay = "RNA", slot = "counts"),
  iter = 10,
  includePCs = 14,
  dims = 20
)

srat$scDblFinder.score <- sce.lr$scDblFinder.score
srat$scDblFinder.class <- sce.lr$scDblFinder.class

p_sc.class <- DimPlot(srat, group.by = 'scDblFinder.class') +
  ggtitle('scDblrFinder class for LR')
p_sc.score <- FeaturePlot(srat, features = 'scDblFinder.score') +
  ggtitle('scDblFinder score for LR')

p_sc.class | p_sc.score

srat <- subset(srat, subset = scDblFinder.class == 'singlet')

#Enrichment analysis ----
##Differently expressed genes ----
## first step: find the differentially expressed genes this
## takes ~ 7min
de_genes <- FindAllMarkers(
  srat,
  assay = "RNA",
  only.pos = FALSE,
  min.pct = 0.1,
  logfc.threshold = log(1.5)
)
## The resulting data.frame is too unwieldy to work with,
## the top 10 genes per cluster are enough:

de_genes %>%
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

DoHeatmap(
  subset(srat, downsample = 50),
  features = top10$gene,
  group.colors = cluster.colors,
  assay = "RNA",
  slot = "scale.data"
) +
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
#enrichment <- compareCluster(
#  clusters,
#  fun = "enrichGO",
#  OrgDb = "org.Hs.eg.db",
#  keyType = "ENSEMBL",
#  ont = "BP",
#  minGSSize = 3,
#  pAdjustMethod = "BH",
#  pvalueCutoff = 0.01,
#  qvalueCutoff = 0.05,
#  readable = TRUE
#)

## so save it:
file <- paste0(CFG$output_dir, "/enrichment.rds")
#saveRDS(file = file, enrichment)

## show results:
clusterProfiler::dotplot(enrichment,
                         showCategory = 4,
                         title = "GO biological process",
                         font = 5)


##Saved srat ----
file <- paste0(CFG$output_dir, "/srat-enriched.rds")
#saveRDS(file = file, srat)
#srat <- readRDS(paste0(CFG$output_dir, "/enrichment.rds"))

#Cell identification ----
# (down)load the reference data:
hpca <- HumanPrimaryCellAtlasData()

singler_hpca <- SingleR(
  test = GetAssayData(srat, assay = "RNA", slot = "data"),
  ref = hpca,
  labels = hpca$label.main
)

plotScoreHeatmap(
  singler_hpca,
  show.labels = TRUE,
  max.labels = 100,
  show.pruned = FALSE,
  order.by = "clusters",
  clusters = srat@meta.data[[clus.column]],
  annotation_colors = list(Clusters = cluster.colors)
)

##Add metadata
srat@meta.data$singler <- singler_hpca$labels
cbind(table(singler_hpca$labels))


typenames <- unique(srat@meta.data$singler)
singler.colors <- Polychrome::palette36.colors(length(typenames))
names(singler.colors) <- typenames

p_type2 <- DimPlot(
  srat,
  reduction = "umap",
  pt.size = 0.5,
  group.by = "singler",
  cols = singler.colors
)

## show the final result, this time together with the other
## characteristics just because we can :-)
p_clus | p_type2


##Saved srat ----
file <- paste0(CFG$output_dir, "/srat-cellsidentified.rds")
#saveRDS(file = file, srat)
#srat <- readRDS(file)

#To investigate what cells are the ones with high NA in our transcriptomic analysis
#srat@meta.data[, 'NAcontent'] <-
#  metadata.cells[gsub('LX279_', '', colnames(srat)), 'NAs']
#FeaturePlot(srat, features = 'NAcontent') | p_singler
#FeaturePlot(srat, features = 'NAcontent') | FeaturePlot(srat,
#                                                        features = 'nCount_RNA')

#ggplot(srat@meta.data, aes(x = NAcontent, y = nCount_RNA)) +
#  geom_point(size = 1) +
#  labs(title = "Scatterplot of percentage of non-covered SNVs against number of counts",
#         x = "Percentage of NAs",
#         y = "nCount_RNA") + geom_smooth(method = "lm", color = "red")


#Tumor cells identification ----
## overview of numbers per celltype:
table(srat@meta.data$SCT_snn_res.1.4)


## which cells do we think we can trust as being healthy:
ref_types <- c('B_cell', 'T_cells', 'NK_cell')
maybe_tumor <- !(srat$singler %in% ref_types)
table(maybe_tumor)

## infercnv needs them as an external data file in a
## specific format, create that:
df <- data.frame(cell = colnames(srat), type = srat$singler)

celltypes_file <- paste0(CFG$output_dir, "/celltypes.txt")
write.table(
  df,
  file = celltypes_file,
  sep = "\t",
  quote = FALSE,
  na = "",
  row.names = FALSE,
  col.names = FALSE
)
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
geneorder <- geneorder[order(geneorder$numeric_chr, geneorder$start), ]
geneorder <- geneorder[!duplicated(geneorder$gene_name), ]
geneorder <- with(geneorder,
                  data.frame(
                    gene = gene_name,
                    chromosome = chr,
                    start = start,
                    end = end
                  ))
## @@@@ WRONG DIR, FIX: geneorder_file <-
## paste0(CFG$data_dir, '/geneorder.txt')
geneorder_file <- paste0(CFG$output_dir, "/geneorder.txt")
write.table(
  geneorder,
  file = geneorder_file,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  na = ""
)

raw_counts <- GetAssayData(srat, assay = "RNA", slot = "counts")

common <- intersect(rownames(raw_counts), geneorder$gene)

## set up the infercnv object with all necessary
## information:
infercnv_obj <- infercnv::CreateInfercnvObject(
  delim = "\t",
  raw_counts_matrix = raw_counts[common, ],
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
infercnv_obj <- infercnv::run(
  infercnv_obj,
  out_dir = outdir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  analysis_mode = "samples",
  num_threads = 2,
  output_format = "png",
  window_length = 101,
  save_final_rds = TRUE,
  plot_steps = FALSE,
  cutoff = 0.1,
  sd_amplifier = 1.5,
  HMM = FALSE,
  write_expr_matrix = TRUE
)
### @@@ this has changed in version 1.14, so the names of
### files are different again! FIX THIS
refexp <-
  as.matrix(read.delim(
    file = paste0(outdir, "/infercnv.references.txt"),
    sep = " "
  ))

obsexp <-
  as.matrix(read.delim(
    file = paste0(outdir, "/infercnv.observations.txt"),
    sep = " "
  ))

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

ScoreCNV <- function(adj_expr,
                     length = 70,
                     SDs = 1.5) {
  apply(adj_expr,
        2,
        .count_genes_in_runs,
        length = length,
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

#Based on this plot we can choose a cutoff
cutoff <- 700

allscore <- c(obsscore, refscore)
allscore <- allscore[colnames(srat)]  # put in order of srat object.

srat <- AddMetaData(srat, allscore, col.name = "aneuploidy_score")

aneuploid <- ifelse(
  srat@meta.data$singler == "Monocyte",
  "no",
  ifelse(srat@meta.data$aneuploidy_score < cutoff, "no", "yes")
)
srat <- AddMetaData(srat, col.name = "is_aneuploid", aneuploid)

## Let's show all we have. If it's too much information,
## leave some of it out.

p_pat <- DimPlot(srat,
                 reduction = "umap",
                 pt.size = 0.5,
                 group.by = "patient_ont")

p_celltype <- DimPlot(
  srat,
  reduction = "umap",
  pt.size = 0.5,
  group.by = "singler",
  cols = singler.colors
)

p_aneup <- FeaturePlot(
  srat,
  reduction = "umap",
  pt.size = 0.5,
  features = "aneuploidy_score",
  order = TRUE
)

p_isaneup <- DimPlot(
  srat,
  reduction = "umap",
  pt.size = 0.5,
  group.by = "is_aneuploid",
  cols = c("red", "green")
)

p_count <- FeaturePlot(
  srat,
  reduction = "umap",
  pt.size = 0.5,
  features = "nCount_SCT",
  order = TRUE,
  label = TRUE
) + ggtitle("Counts per cell")

(p_pat | p_celltype) / (p_aneup | p_isaneup)

p_count | p_aneup

srat@meta.data <- srat@meta.data %>%
  mutate(
    reference = case_when(
      patient_ont == '0' & is_aneuploid == 'yes' ~ 'p0.tumor',
      patient_ont == '0' & is_aneuploid == 'no' ~ 'p0.healthy',
      patient_ont == '1' & is_aneuploid == 'yes' ~ 'p1.tumor',
      patient_ont == '1' & is_aneuploid == 'no' ~ 'p1.healthy',
      TRUE ~ 'NA'
    )
  )

file <- paste0(CFG$output_dir, "/srat-aneuploidy.rds")
#saveRDS(file=file, srat)

#SNV ----
#Creating
srat.meta <- data.frame(srat@meta.data)

#String manipulation for alt & ref. We want a matrix without row/colnames
#and only counts as values
alt <- read.table(file = '~/sop_input/alt.mtx', header = TRUE)
alt$ID <- NULL
alt <- alt %>% unite(col = 'SNV', CHROM:ALT, sep = ':')
rownames(alt) <- alt$SNV
rownames(alt) <- paste0('alt:', rownames(alt))
alt$SNV <- NULL
alt$QUAL <- alt[alt$QUAL > CFG$QUAL]
alt$QUAL <- NULL
alt <- as.matrix(alt)

ref <- read.table(file = '~/sop_input/ref.mtx', header = TRUE)
ref$ID <- NULL
ref <- ref %>% unite(col = 'SNV', CHROM:ALT, sep = ':')
rownames(ref) <- ref$SNV
rownames(ref) <- paste0('ref:', rownames(ref))
ref$SNV <- NULL
ref <- ref[ref$QUAL > CFG$QUAL]
ref$QUAL <- NULL
ref <- as.matrix(ref)

counts.per.cell <- (rbind(alt, ref))
counts.per.cell <- as.matrix(counts.per.cell)

#Only gonna use nuclear reads
mitochondrial.snvs <- grepl("^alt:chrM|^ref:chrM", rownames(counts.per.cell))

mito.matrix <- counts.per.cell[mitochondrial.snvs, ]
nuc.matrix <- counts.per.cell[!mitochondrial.snvs, ]

nuc.matrix <- nuc.matrix[, colnames(nuc.matrix) %in% gsub('LX279_', '', colnames(srat))]
dim(nuc.matrix)

srat.snvs <- Seurat::CreateSeuratObject(counts = nuc.matrix,
                                        project = "ONT_snvs",
                                        meta = srat.meta)
colnames(srat.snvs) <- paste0('LX279_', colnames(srat.snvs))

srat.snvs <- NormalizeData(srat.snvs, normalization.method = "LogNormalize")

srat.snvs <- ScaleData(srat.snvs,
                       features = rownames(srat.snvs),
                       verbose = FALSE)

srat.snvs <- FindVariableFeatures(srat.snvs)

top10.snvs <- head(VariableFeatures(srat.snvs), 10)


srat.snvs <- suppressMessages(
  SCTransform(
    srat.snvs,
    vars.to.regres = NULL,
    vst.flavor = 'v2',
    ## version that solved some problems
    verbose = FALSE,
    variable.features.n = 3000
  )
)

p.varsnvs <- VariableFeaturePlot(srat.snvs, selection.method = 'sct')

p.varsnvs <- LabelPoints(plot = p.varsnvs,
                         points = top10.snvs,
                         repel = TRUE)

p.varsnvs

DefaultAssay(srat.snvs)
srat.snvs <- RunPCA(srat.snvs, npcs = 50)

elbow <- ElbowPlot(srat.snvs, reduction = "pca", ndims = 50)
elbow #geom_vline(xintercept = c(seq(from = 0, to = 50, by = 1)))

DimHeatmap(srat.snvs, dims = (15:27), cells = 100)

srat.snvs <- RunUMAP(srat.snvs, dims = 1:CFG$ndims.snvs, assay = 'RNA')

DimPlot(srat.snvs, group.by = 'reference')

srat.snvs <- FindNeighbors(srat.snvs, dims = 1:CFG$ndims.snvs)
## Computing nearest neighbor graph
## Computing SNN

srat.snvs <- FindClusters(srat.snvs,
                          resolution = seq(0.4, 2, 0.2),
                          algorithm = 1)

columns <- grep("SCT_snn_res", names(srat.snvs@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the
## cluster ids)
for (col in columns) {
  cat("\n====\n", col, ":")
  print(table(srat.snvs@meta.data[[col]]))
}

clustree::clustree(
  srat.snvs,
  prefix = "SCT_snn_res.",
  node_size = 2,
  edge_width = 0.5
)

#Once again, decision based on clustree
CFG$reso.snv
clus.column.snv <- paste0("SCT_snn_res.", CFG$reso.snv)
Idents(srat) <- clus.column.snv

clust.snvs <- sort(unique(srat.snvs@meta.data[[clus.column]]))
## switching color scheme to make it look different:
clust.cols.snvs <- Polychrome::alphabet.colors(length(clust.snvs))
names(clust.cols.snvs) <- as.character(clust.snvs)

## and show result with the other umaps to get your
## bearings:
p_clus_snv <- DimPlot(
  srat.snvs,
  reduction = "umap",
  pt.size = 0.5,
  group.by = clus.column.snv,
  cols = clust.cols.snvs,
  label = TRUE
)

p_clus_snv
gc()

de_snvs <- FindAllMarkers(
  srat.snvs,
  assay = "RNA",
  min.pct = 0.1,
  logfc.threshold = log(1.5)
)
## The resulting data.frame is too unwieldy to work with,
## the top 10 genes per cluster are enough:

de_snvs %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = pct.1) -> top3.de.snvs
head(top10.de.snvs)

##Obtaining SCT selected snvs matrix ----
SCT.scaled.snvs <- srat.snvs@assays$RNA$scale.data
genes_union <- union(top3.de.snvs.pct$gene, top10.pat.snvs$gene)
subset_SCT <- SCT.scaled.snvs[match(genes_union, rownames(SCT.scaled.snvs)), , drop = FALSE]
colnames(subset_SCT) <- gsub('LX279_', '', colnames(subset_SCT))
dim(subset_SCT)

##Ireen saving the day ----
Idents(srat.snvs) <- "SCT_snn_res.1.4"
srat.snvs.small <- subset(srat.snvs, downsample = 50)
# Which genes to plot
genes_union <- union(top3.de.snvs.pct$gene, top10.pat.snvs$gene)
#Get matrix for Heatmap
dataset <- as.matrix(GetAssayData(srat.snvs.small, assay = "RNA", slot = "scale.data")[genes_union, ])



DoHeatmap(
  subset(srat.snvs, downsample = 50),
  features = top3.de.snvs.pct$gene,
  #group.colors = clust.cols.snvs,
  assay = "RNA",
  slot = "scale.data",
  group.by = c('reference', 'ident')
) +
  theme(axis.text.y = element_text(size = 6))

#saveRDS(srat.snvs, file = '~/srat.snvs-clusters.rds')
srat.snvs <- readRDS(file = '~/srat.snvs-clusters.rds')

srat.tumor.p0 <- subset(srat.snvs, subset = reference == 'p0.tumor')
srat.tumor.p1 <- subset(srat.snvs, subset = reference == 'p1.tumor')

srat.ge.tumor.p0 <- subset(srat, subset = reference == 'p0.tumor')
srat.ge.tumor.p1 <- subset(srat, subset = reference == 'p1.tumor')

#Add clust based on snvs to GE. Check reso!!!
CFG$reso.snv
srat <- AddMetaData(srat,
                    metadata = srat.snvs$SCT_snn_res.1.8,
                    col.name = 'snvs.clustering')
DimPlot(srat,
        group.by = 'snvs.clustering',
        label = TRUE,
        cols = rainbow(length(unique(
          srat$snvs.clustering
        ))))

#Simplify annotation
cell_percentage_snv <- sum(table(srat.snvs$singler))

srat.snvs$singler <- ifelse(
  table(srat.snvs$singler)[srat.snvs$singler] / cell_percentage_snv < 0.01,
  'Other',
  srat.snvs$singler
) + 

table(srat.snvs$singler)

DimPlot(srat.snvs, group.by = 'singler', pt.size = 0.5)

#Tumor 0 ----
#SNV-based clustering
srat.tumor.p0 <- NormalizeData(srat.tumor.p0,
                               normalization.method = "LogNormalize",
                               assay = 'RNA')

srat.tumor.p0 <- ScaleData(
  srat.tumor.p0,
  features = rownames(srat.tumor.p0),
  verbose = FALSE,
  assay = 'RNA'
)

srat.tumor.p0 <- FindVariableFeatures(srat.tumor.p0, assay = 'RNA')


srat.tumor.p0 <- suppressMessages(
  SCTransform(
    srat.tumor.p0,
    vars.to.regres = NULL,
    vst.flavor = 'v2',
    ## version that solved some problems
    verbose = FALSE,
    variable.features.n = 3000
  )
)


top10.snvs.p0 <- head(VariableFeatures(srat.tumor.p0), 10)

p0.varsnvs <- VariableFeaturePlot(srat.tumor.p0, selection.method = 'sct')

p0.varsnvs <- LabelPoints(plot = p0.varsnvs,
                          points = top10.snvs.p0,
                          repel = TRUE)

p0.varsnvs

DefaultAssay(srat.tumor.p0)
srat.tumor.p0 <- RunPCA(srat.tumor.p0, npcs = 50)

elbow <- ElbowPlot(srat.tumor.p0, reduction = "pca", ndims = 50)
elbow #geom_vline(xintercept = c(seq(from = 0, to = 50, by = 1)))

#Play with ndims here to see which one is the optimal
DimHeatmap(srat.tumor.p0, dims = (14:28), cells = 100)

srat.tumor.p0 <- RunUMAP(srat.tumor.p0, dims = 1:CFG$ndims.p0)
srat.tumor.p0 <- FindNeighbors(srat.tumor.p0, dims = 1:CFG$ndims.p0)
## Computing nearest neighbor graph
## Computing SNN

## Using clustree ----
srat.tumor.p0 <- FindClusters(srat.tumor.p0,
                              resolution = seq(0.4, 2, 0.2),
                              algorithm = 1)

snv.cols.p0 <- grep("SCT_snn_res", names(srat.snvs@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the
## cluster ids)
for (col in snv.cols.p0) {
  cat("\n====\n", col, ":")
  print(table(srat.tumor.p0@meta.data[[col]]))
}

clustree::clustree(
  srat.tumor.p0,
  prefix = "SCT_snn_res.",
  node_size = 2,
  edge_width = 0.5
)

CFG$reso.snv.p0
clus.column.snv.p0 <- paste0("SCT_snn_res.", CFG$reso.snv.p0)
Idents(srat.tumor.p0) <- clus.column.snv.p0

clust.snvs.p0 <- sort(unique(srat.tumor.p0@meta.data[[clus.column.snv.p0]]))
## switching color scheme to make it look different:
cluster.colors.p0.snvs <- Polychrome::alphabet.colors(length(clust.snvs.p0))
names(cluster.colors.p0.snvs) <- as.character(clust.snvs.p0)

#Extra Plots
p_clus_snv_p0 <- DimPlot(
  srat.tumor.p0,
  reduction = "umap",
  pt.size = 0.5,
  group.by = clus.column.snv.p0,
  cols = cluster.colors,
  label = TRUE
)
p_clus_snv_p0 |
  FeaturePlot(srat.tumor.p0, features = 'aneuploidy_score', order = TRUE)

VlnPlot(srat.tumor.p0, features = c('nCount_RNA', 'aneuploidy_score'))

#I dont know how to make this variable, but check the reso, so it matches
#the Idents in the GE
CFG$reso.snv.p0
Idents(srat.ge.tumor.p0) <- srat.tumor.p0$SCT_snn_res.1.4

##DE SNVs based on SNV clust ----
de_snvs_snv_p0 <- FindAllMarkers(
  srat.tumor.p0,
  assay = "RNA",
  min.pct = 0.1,
  logfc.threshold = log(1.5)
)
## The resulting data.frame is too unwieldy to work with,
## the top 10 genes per cluster are enough:

de_snvs_snv_p0 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10.de.snv.snv
head(top10.de.snv.snv)

DoHeatmap(
  subset(srat.tumor.p0, downsample = 50),
  features = top10.de.snv.snv$gene,
  group.colors = clust.cols.snvs,
  assay = "RNA",
  slot = "scale.data"
) +
  theme(axis.text.y = element_text(size = 6))

##DEG based on SNV clusters----
de_genes_snvs_p0 <- FindAllMarkers(
  srat.ge.tumor.p0,
  assay = "RNA",
  min.pct = 0.1,
  logfc.threshold = log(1.5)
)
de_genes_snvs_p0 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) -> top10.p0.ge
head(top10.p0.ge)

de_genes_snvs_p0 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10.p0.ge.top_n
head(top10.p0.ge.top_n)

DoHeatmap(
  subset(srat.ge.tumor.p0, downsample = 50),
  features = top10.p0.ge$gene,
  group.colors = cluster.colors,
  assay = "RNA",
  slot = "scale.data"
) +
  theme(axis.text.y = element_text(size = 6))

#Check reso and change accordingly
CFG$reso.snv.p0
srat.ge.tumor.p0 <- AddMetaData(srat.ge.tumor.p0,
                                metadata = srat.tumor.p0$SCT_snn_res.1.4,
                                col.name = 'snv_clusters')

DoHeatmap(
  subset(srat.ge.tumor.p0, downsample = 50),
  features = top10.p0.ge$gene,
  group.colors = cluster.colors,
  assay = "RNA",
  slot = "scale.data"
) +
  theme(axis.text.y = element_text(size = 6))

VlnPlot(srat.ge.tumor.p0,
        features = c('nCount_RNA', 'aneuploidy_score'))


##GE based clusters (tumor only) ----
srat.ge.tumor.p0 <- NormalizeData(srat.ge.tumor.p0,
                                  normalization.method = "LogNormalize",
                                  assay = 'RNA')

srat.ge.tumor.p0 <- ScaleData(
  srat.ge.tumor.p0,
  features = rownames(srat.ge.tumor.p0),
  verbose = FALSE,
  assay = 'RNA'
)

srat.ge.tumor.p0 <- FindVariableFeatures(srat.ge.tumor.p0, assay = 'RNA')


srat.ge.tumor.p0 <- suppressMessages(
  SCTransform(
    srat.ge.tumor.p0,
    vars.to.regres = NULL,
    vst.flavor = 'v2',
    ## version that solved some problems
    verbose = FALSE,
    variable.features.n = 3000
  )
)

top10.ge.p0 <- head(VariableFeatures(srat.ge.tumor.p0), 10)

p0.vargenes <- VariableFeaturePlot(srat.ge.tumor.p0, selection.method = 'sct')

p0.vargenes <- LabelPoints(plot = p0.vargenes,
                           points = top10.ge.p0,
                           repel = TRUE)

p0.vargenes

DefaultAssay(srat.ge.tumor.p0)
srat.ge.tumor.p0 <- RunPCA(srat.ge.tumor.p0, npcs = 50)

elbow <- ElbowPlot(srat.ge.tumor.p0, reduction = "pca", ndims = 50)
elbow #geom_vline(xintercept = c(seq(from = 0, to = 50, by = 1)))

DimHeatmap(srat.ge.tumor.p0, dims = (12:26), cells = 100)
#Quickly check number of ndims
CFG$ndims.p0
srat.ge.tumor.p0 <- RunUMAP(srat.ge.tumor.p0, dims = 1:CFG$ndims.p0)
srat.ge.tumor.p0 <- FindNeighbors(srat.ge.tumor.p0, dims = 1:CFG$ndims.p0)
## Computing nearest neighbor graph
## Computing SNN

srat.ge.tumor.p0 <- FindClusters(srat.ge.tumor.p0,
                                 resolution = seq(0.4, 2, 0.2),
                                 algorithm = 1)

snv.cols.p0 <- grep("SCT_snn_res", names(srat.ge.tumor.p0@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the
## cluster ids)
for (col in snv.cols.p0) {
  cat("\n====\n", col, ":")
  print(table(srat.ge.tumor.p0@meta.data[[col]]))
}

clustree::clustree(
  srat.ge.tumor.p0,
  prefix = "SCT_snn_res.",
  node_size = 2,
  edge_width = 0.5
)

#Check reso and change to whatever works better
CFG$reso.ge.p0
clus.column.ge.p0 <- paste0("SCT_snn_res.", CFG$reso.ge.p0)
Idents(srat.ge.tumor.p0) <- clus.column.ge.p0

clust.ge.p0 <- sort(unique(srat.ge.tumor.p0@meta.data[[clus.column.ge.p0]]))
## switching color scheme to make it look different:
cluster.colors.p0.ge <- Polychrome::alphabet.colors(length(clust.ge.p0))
names(cluster.colors.p0.ge) <- as.character(clust.ge.p0)

#Extra Plots
p_clus_ge_p0 <- DimPlot(
  srat.ge.tumor.p0,
  reduction = "umap",
  pt.size = 0.5,
  group.by = clus.column.ge.p0,
  cols = cluster.colors.p0.ge,
  label = TRUE
)
p_clus_ge_p0 |
  FeaturePlot(srat.ge.tumor.p0, features = 'aneuploidy_score')

VlnPlot(srat.ge.tumor.p0,
        features = c('nCount_RNA', 'aneuploidy_score'))

#Check reso and change accordingly
CFG$reso.snv.p0
srat.ge.tumor.p0 <- AddMetaData(srat.ge.tumor.p0,
                                metadata = srat.tumor.p0$SCT_snn_res.1.4,
                                col.name = 'snv_clusters')
##Annotation of SNV based clusters on GE based UMAP
DimPlot(
  srat.ge.tumor.p0,
  group.by = 'snv_clusters',
  cols = cluster.colors.p0.snvs,
  pt.size = 0.5
)

CFG$reso.ge.p0
Idents(srat.ge.tumor.p0) <- srat.ge.tumor.p0$SCT_snn_res.1
de_genes_ge <- FindAllMarkers(
  srat.ge.tumor.p0,
  assay = "RNA",
  min.pct = 0.1,
  logfc.threshold = log(1.5)
)
de_genes_ge %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) -> top10.p0.ge.only
head(top10.p0.ge.only)

DoHeatmap(
  subset(srat.ge.tumor.p0, downsample = 50),
  features = top10.p0.ge.only$gene,
  group.colors = cluster.colors,
  assay = "RNA",
  slot = "scale.data"
) +
  theme(axis.text.y = element_text(size = 6))

VlnPlot(srat.ge.tumor.p0,
        features = c('nCount_RNA', 'aneuploidy_score'))


##Pathway analysis ----
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
}

if (websiteLive)
  dbs <- listEnrichrDbs()

if (websiteLive)
  head(dbs)

#Select dbs and which cluster you want to look at
dbs <- c(
  "GO_Molecular_Function_2023",
  "GO_Cellular_Component_2023",
  "GO_Biological_Process_2023"
)
if (websiteLive) {
  enriched <- enrichr(head(de_genes_snvs_p0[de_genes_snvs_p0$cluster == 4, 'gene'], n = 500), dbs)
} #Choose genes to see enrichment

if (websiteLive) {
  plotEnrich(
    enriched[[3]],
    showTerms = 20,
    numChar = 40,
    y = "Count",
    orderBy = "P.value"
  )
}

#Tumor 1 ----
#SNV-based clustering
srat.tumor.p1 <- NormalizeData(srat.tumor.p1,
                               normalization.method = "LogNormalize",
                               assay = 'RNA')

srat.tumor.p1 <- ScaleData(
  srat.tumor.p1,
  features = rownames(srat.tumor.p1),
  verbose = FALSE,
  assay = 'RNA'
)

srat.tumor.p1 <- FindVariableFeatures(srat.tumor.p1, assay = 'RNA')


srat.tumor.p1 <- suppressMessages(
  SCTransform(
    srat.tumor.p1,
    vars.to.regres = NULL,
    vst.flavor = 'v2',
    ## version that solved some problems
    verbose = FALSE,
    variable.features.n = 3000
  )
)


top10.snvs.p1 <- head(VariableFeatures(srat.tumor.p1), 10)

p1.varsnvs <- VariableFeaturePlot(srat.tumor.p1, selection.method = 'sct')

p1.varsnvs <- LabelPoints(plot = p1.varsnvs,
                          points = top10.snvs.p1,
                          repel = TRUE)

p1.varsnvs

DefaultAssay(srat.tumor.p1)
srat.tumor.p1 <- RunPCA(srat.tumor.p1, npcs = 50)

elbow <- ElbowPlot(srat.tumor.p1, reduction = "pca", ndims = 50)
elbow #geom_vline(xintercept = c(seq(from = 0, to = 50, by = 5)))

#Play with ndims here to see which one is the optimal
DimHeatmap(srat.tumor.p1, dims = (14:28), cells = 100)

CFG$ndims.p1
srat.tumor.p1 <- RunUMAP(srat.tumor.p1, dims = 1:CFG$ndims.p1)
srat.tumor.p1 <- FindNeighbors(srat.tumor.p1, dims = 1:CFG$ndims.p1)
## Computing nearest neighbor graph
## Computing SNN

## Using clustree ----
srat.tumor.p1 <- FindClusters(srat.tumor.p1,
                              resolution = seq(0.4, 2, 0.2),
                              algorithm = 1)

snv.cols.p1 <- grep("SCT_snn_res", names(srat.snvs@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the
## cluster ids)
for (col in snv.cols.p1) {
  cat("\n====\n", col, ":")
  print(table(srat.tumor.p1@meta.data[[col]]))
}

clustree::clustree(
  srat.tumor.p1,
  prefix = "SCT_snn_res.",
  node_size = 2,
  edge_width = 0.5
)

CFG$reso.snv.p1 <- 1
clus.column.snv.p1 <- paste0("SCT_snn_res.", CFG$reso.snv.p1)
Idents(srat.tumor.p1) <- clus.column.snv.p1

clust.snvs.p1 <- sort(unique(srat.tumor.p1@meta.data[[clus.column]]))
## switching color scheme to make it look different:
cluster.colors.p1.snvs <- Polychrome::alphabet.colors(length(clust.snvs.p1))
names(cluster.colors.p1.snvs) <- as.character(clust.snvs.p1)

#Extra Plots
p_clus_snv_p1 <- DimPlot(
  srat.tumor.p1,
  reduction = "umap",
  pt.size = 0.5,
  group.by = clus.column.snv.p1,
  cols = cluster.colors.p1.snvs,
  label = TRUE
)
p_clus_snv_p1 | DimPlot(srat.tumor.p1, group.by = 'Phase')

VlnPlot(srat.tumor.p1, features = c('nCount_RNA', 'aneuploidy_score'))

#I dont know how to make this variable, but check the reso, so it matches
#the Idents in the GE
CFG$reso.snv.p1
Idents(srat.ge.tumor.p1) <- srat.tumor.p1$SCT_snn_res.1

##DE SNVs based on SNV clusts ----
de_snvs_snv_p1 <- FindAllMarkers(
  srat.tumor.p1,
  assay = "RNA",
  min.pct = 0.1,
  logfc.threshold = log(1.5)
)
## The resulting data.frame is too unwieldy to work with,
## the top 10 genes per cluster are enough:

de_snvs_snv_p1 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10.de.snv.snv
head(top10.de.snv.snv)

DoHeatmap(
  subset(srat.tumor.p1, downsample = 50),
  features = top10.de.snv.snv$gene,
  group.colors = clust.cols.snvs,
  assay = "RNA",
  slot = "scale.data"
) +
  theme(axis.text.y = element_text(size = 6))


##DEG based on SNV clusters----
de_genes_snvs_p1 <- FindAllMarkers(
  srat.ge.tumor.p1,
  assay = "RNA",
  min.pct = 0.1,
  logfc.threshold = log(1.5)
)
de_genes_snvs_p1 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) -> top10.p1.ge
head(top10.p1.ge)

de_genes_snvs_p1 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10.p1.ge.top_n
head(top10.p1.ge.top_n)

#Check reso and change accordingly
CFG$reso.snv.p1
srat.ge.tumor.p1 <- AddMetaData(srat.ge.tumor.p1,
                                metadata = srat.tumor.p1$SCT_snn_res.1,
                                col.name = 'snv_clusters')

DoHeatmap(
  subset(srat.ge.tumor.p1, downsample = 50),
  features = top10.p1.ge$gene,
  group.colors = cluster.colors,
  assay = "RNA",
  slot = "scale.data"
) +
  theme(axis.text.y = element_text(size = 6))

VlnPlot(srat.ge.tumor.p1,
        features = c('nCount_RNA', 'aneuploidy_score'))


##GE based clusters (tumor only) ----
srat.ge.tumor.p1 <- NormalizeData(srat.ge.tumor.p1,
                                  normalization.method = "LogNormalize",
                                  assay = 'RNA')

srat.ge.tumor.p1 <- ScaleData(
  srat.ge.tumor.p1,
  features = rownames(srat.ge.tumor.p1),
  verbose = FALSE,
  assay = 'RNA'
)

srat.ge.tumor.p1 <- FindVariableFeatures(srat.ge.tumor.p1, assay = 'RNA')


srat.ge.tumor.p1 <- suppressMessages(
  SCTransform(
    srat.ge.tumor.p1,
    vars.to.regres = NULL,
    vst.flavor = 'v2',
    ## version that solved some problems
    verbose = FALSE,
    variable.features.n = 3000
  )
)


top10.ge.p1 <- head(VariableFeatures(srat.ge.tumor.p1), 10)

p1.vargenes <- VariableFeaturePlot(srat.ge.tumor.p1, selection.method = 'sct')

p1.vargenes <- LabelPoints(plot = p1.vargenes,
                           points = top10.ge.p1,
                           repel = TRUE)

p1.vargenes

DefaultAssay(srat.ge.tumor.p1)
srat.ge.tumor.p1 <- RunPCA(srat.ge.tumor.p1, npcs = 50)

elbow <- ElbowPlot(srat.ge.tumor.p1, reduction = "pca", ndims = 50)
elbow #geom_vline(xintercept = c(seq(from = 0, to = 50, by = 1)))

DimHeatmap(srat.ge.tumor.p1, dims = (12:26), cells = 100)
#Quickly check number of ndims
CFG$ndims.p1
srat.ge.tumor.p1 <- RunUMAP(srat.ge.tumor.p1, dims = 1:CFG$ndims.p1)
srat.ge.tumor.p1 <- FindNeighbors(srat.ge.tumor.p1, dims = 1:CFG$ndims.p1)
## Computing nearest neighbor graph
## Computing SNN

srat.ge.tumor.p1 <- FindClusters(srat.ge.tumor.p1,
                                 resolution = seq(0.4, 2, 0.2),
                                 algorithm = 1)

snv.cols.p1 <- grep("SCT_snn_res", names(srat.ge.tumor.p1@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the
## cluster ids)
for (col in snv.cols.p1) {
  cat("\n====\n", col, ":")
  print(table(srat.ge.tumor.p1@meta.data[[col]]))
}

clustree::clustree(
  srat.ge.tumor.p1,
  prefix = "SCT_snn_res.",
  node_size = 2,
  edge_width = 0.5
)

#Check reso and change to whatever works better
CFG$reso.ge.p1
clus.column.ge.p1 <- paste0("SCT_snn_res.", CFG$reso.ge.p1)
Idents(srat.ge.tumor.p1) <- clus.column.ge.p1

clust.ge.p1 <- sort(unique(srat.ge.tumor.p1@meta.data[[clus.column.ge.p1]]))
## switching color scheme to make it look different:
cluster.colors.p1.ge <- Polychrome::alphabet.colors(length(clust.ge.p1))
names(cluster.colors.p1.ge) <- as.character(clust.ge.p1)

#Extra Plots
p_clus_ge_p1 <- DimPlot(
  srat.ge.tumor.p1,
  reduction = "umap",
  pt.size = 0.5,
  group.by = clus.column.ge.p1,
  cols = cluster.colors.p1.ge,
  label = TRUE
)
p_clus_ge_p1 |
  DimPlot(srat.ge.tumor.p1, group.by = 'Phase', pt.size = 0.5)

VlnPlot(srat.ge.tumor.p1,
        features = c('nCount_RNA', 'aneuploidy_score'))

##DE SNVs based on GE clusters ----
CFG$reso.ge.p1
clus.column.ge.p1 <- paste0("SCT_snn_res.", CFG$reso.ge.p1)
Idents(srat.tumor.p1) <- clus.column.ge.p1

de_snvs_ge_p0 <- FindAllMarkers(
  srat.tumor.p1,
  assay = "RNA",
  min.pct = 0.1,
  logfc.threshold = log(1.5)
)
de_snvs_ge_p0 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) -> top10.p0.snv

##Annotation of SNV based clusters on GE based UMAP
DimPlot(
  srat.ge.tumor.p1,
  group.by = 'snv_clusters',
  cols = cluster.colors.p1.snvs,
  pt.size = 0.5
)

##DEG based on GE clusts ----
CFG$reso.ge.p1
Idents(srat.ge.tumor.p1) <- srat.ge.tumor.p1$SCT_snn_res.1.4
de_genes_ge_p1 <- FindAllMarkers(
  srat.ge.tumor.p1,
  assay = "RNA",
  min.pct = 0.1,
  logfc.threshold = log(1.5)
)
de_genes_ge_p1 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) -> top10.p1.ge.only
head(top10.p1.ge.only)

DoHeatmap(
  subset(srat.ge.tumor.p1, downsample = 50),
  features = top10.p1.ge.only$gene,
  group.colors = cluster.colors,
  assay = "RNA",
  slot = "scale.data"
) +
  theme(axis.text.y = element_text(size = 6))

VlnPlot(srat.ge.tumor.p1,
        features = c('nCount_RNA', 'aneuploidy_score'))


##Pathway analysis ----
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
}

if (websiteLive)
  dbs <- listEnrichrDbs()

if (websiteLive)
  head(dbs)


#Select dbs and which cluster you want to look at
dbs <- c(
  "GO_Molecular_Function_2023",
  "GO_Cellular_Component_2023",
  "GO_Biological_Process_2023"
)
if (websiteLive) {
  enriched <- enrichr(head(de_genes_snvs_p1[de_genes_snvs_p1$cluster == 5, 'gene'], n = 500), dbs)
} #Choose genes to see enrichment

if (websiteLive) {
  plotEnrich(
    enriched[[3]],
    showTerms = 20,
    numChar = 40,
    y = "Count",
    orderBy = "P.value"
  )
}

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
# Supongamos que tienes tus datos en srat.snvs$SCT_snn_res.1.8 y srat$SCT_snn_res.1.4

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
