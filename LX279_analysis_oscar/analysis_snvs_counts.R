#remove.packages('Matrix')
#devtools::install_version("Matrix", version = "1.5.0", repos = "http://cran.us.r-project.org")

#detach("package:Matrix", unload = TRUE)

local.lib <- '/hpc/pmc_holstege/rstudio/oscar_clonetracer/R/x86_64-pc-linux-gnu-library/4.4'

library(Seurat)
library(Matrix)
library(ComplexHeatmap)
library(magick)
library(reshape2)
library(ggplot2)
library(impute)
library(matrixStats)
library(stringr)
library(reshape2) #for melt function
library(cowplot) #for plot_grid
library(tidyr) #for unite func
library(RColorBrewer)
library(umap)
library(tidyverse)
library("ggsci")
library("gridExtra")

#Loading data ----
#List all the cutoffs in the beginning for easier reading
cutoff <- list()
cutoff$QUAL <- 30 #Higher than 30
cutoff$nonzero <- 10 #Maximum of n% of NAs, in other words, highly covered
cutoff$cells <- 80 #Cells can have a maximum of n% of NAs for the total of SNVs
cutoff$sds <- 0.5 #SNVs selected must have more than n% SD
cutoff$cumulative.var <- 80 #For selection of PCs, amount of cumulative Var needed 

#Will use srat from ONT to filter out cells based on RNA assay
srat.lr <- readRDS(file = '~/analysis_nanopore_LX279/genes/srat-aneuploidy.rds')

#Check whether underrepresented cell types have been collapsed into a single class 'Other'
table(srat.lr$singler) 

cell_percentage_ont <- sum(table(srat.lr$singler))

# Renaming just in case it wasn't done
srat.lr$singler <- ifelse(
  table(srat.lr$singler)[srat.lr$singler] / cell_percentage_ont < 0.01,
  'Other',
  srat.lr$singler
)

table(srat.lr$singler)

#Store metadata from RNA based srat
srat.meta <- data.frame(row.names = gsub('LX279_','',colnames(srat.lr)),
                        patient = srat.lr@meta.data$patient_ont,
                        is_aneuploid = srat.lr@meta.data$is_aneuploid,
                        cell_type = srat.lr@meta.data$singler)

#vafs.merged<-read.table(file = '~/snv_calling/processing.vafs/vafs.merged.tsv', header = TRUE)
alt <- read.table(file = '~/LX279_snvs/sop_input/alt.mtx', header = TRUE)

#Modify to keep only counts and SNV info
alt$ID <- NULL
alt<- alt %>% unite(col='SNV',CHROM:ALT, sep=':')
rownames(alt) <- alt$SNV
rownames(alt) <- paste0('alt:', rownames(alt))
alt$SNV <- NULL
alt$QUAL <- alt[alt$QUAL > cutoff$QUAL]
alt$QUAL <- NULL
alt <- as.matrix(alt)

#Same for ref
ref <- read.table(file = '~/LX279_snvs/sop_input/ref.mtx', header = TRUE)
ref$ID <- NULL
ref<- ref %>% unite(col='SNV',CHROM:ALT, sep=':')
rownames(ref) <- ref$SNV
rownames(ref) <- paste0('ref:', rownames(ref))
ref$SNV <- NULL
ref <- ref[ref$QUAL > cutoff$QUAL]
ref$QUAL <- NULL
ref <- as.matrix(ref)

#Create allele count matrix
counts.per.cell <- (rbind(alt, ref))
counts.per.cell <- as.matrix(counts.per.cell)

cnts_per_cell <- Matrix::colSums(counts.per.cell)# note use of Matrix::colSums
cnts_per_snv <- rowSums(counts.per.cell)
summary(cnts_per_cell)
summary(cnts_per_snv)
hist(cnts_per_cell, breaks = 50)
hist(cnts_per_snv, breaks = 1000)

#String manipulation ----
#Split by nuc and mito SNVs
mitochondrial.snvs <- grepl("^alt:chrM|^ref:chrM", rownames(counts.per.cell))

mito.matrix <- counts.per.cell[mitochondrial.snvs, ]
nuc.matrix <- counts.per.cell[!mitochondrial.snvs, ]
#top.matrix <- counts.per.cell[rownames(counts.per.cell) %in% union(top10.de.snvs$gene, top10.pat.snvs$gene), ]

#Genotype ----
bcs <- gsub('LX279_','',colnames(srat.lr))
selected.bcs <- nuc.matrix[, colnames(nuc.matrix) %in% bcs]

#Nuc patient 0
cnts_per_cell <- colSums(selected.bcs)
cnts_per_snv <- rowSums(selected.bcs)
summary(cnts_per_cell)
summary(cnts_per_snv)

#Save information about our matrix as metadata
nuc.unfiltered <- rowMeans(selected.bcs, na.rm = TRUE)
nuc.unfiltered <- as.data.frame(nuc.unfiltered)
colnames(nuc.unfiltered) <- 'mean'
#SD of each SNV
nuc.unfiltered$sd <- rowSds(selected.bcs, na.rm = TRUE)
#Number of cells with supporting reads per SNV
nuc.unfiltered$nonzero <- (apply(selected.bcs, 1, function(row) sum(row != 0, na.rm = TRUE)) / ncol(selected.bcs)) * 100

dim(selected.bcs)
#[1] 2633 3691
summary(colSums(selected.bcs))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#33     891    1454    2036    2660   13313 
summary(rowSums(selected.bcs))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0   135.2   376.0  1250.9   720.0 83313.0
normalized <- LogNormalize(selected.bcs)

#filtered.bcs stores information of the selected.bcs matrix, namely, info
#after deleting low quality cells.
nuc.normalized <- rowMeans(normalized, na.rm = TRUE)
nuc.normalized <- as.data.frame(nuc.normalized)
colnames(nuc.normalized) <- 'mean'
nuc.normalized$sd <- rowSds(normalized, na.rm = TRUE)
nuc.normalized$nonzero <- (apply(normalized, 1, function(row) sum(row != 0, na.rm = TRUE)) / ncol(normalized)) * 100

dim(nuc.normalized)
#[1] 2633 3691
summary(colSums(normalized))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#33     891    1454    2036    2660   13313 
summary(rowSums(normalized))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0   135.2   376.0  1250.9   720.0 83313.0


#selected.snvs is our matrix after taking only high coverage SNVs AFTER
#filtering low qual cells
filtered <- normalized[nuc.normalized$sd > cutoff$sds,]

##Imputation ----
nuc.srat.meta <- srat.meta[rownames(srat.meta) %in% colnames(filtered),]
nuc.srat.meta <- srat.meta[colnames(filtered),]
nuc.srat.annotation <- HeatmapAnnotation(df = nuc.srat.meta, col = list(
  is_aneuploid = c('no' = 'red','yes' = 'blue'),
  patient_ont = c('0' = 'green'),
  cell_type = c('B_cell' = '#FED439FF',
                'Monocyte' = '#709AE1FF',
                'NK_cell' = '#8A9197FF',
                'Other' = '#D2AF81FF',
                'Pro-B_cell_CD34+' = '#FD7446FF',
                'T_cells' = '#D5E4A2FF')
)
)

set.seed(123)
nuc.ht <- ComplexHeatmap::Heatmap(filtered,
                                  cluster_rows = TRUE,
                                  cluster_columns = TRUE,
                                  use_raster = TRUE,
                                  col = brewer.pal(7, "Greens"),
                                  show_row_names = FALSE,
                                  show_column_names = FALSE,
                                  row_title = 'SNVs',
                                  column_title = 'Barcodes',
                                  column_km = 2,
                                  column_km_repeats = 10,
                                  top_annotation = nuc.srat.annotation)
nuc.ht <- draw(nuc.ht)

patient0 <- colnames(filtered)[column_order(nuc.ht)[[2]]]
patient1 <- colnames(filtered)[column_order(nuc.ht)[[1]]]

## Adding metadata to srat ----
patient0 <- paste0('LX279_', patient0)
patient1 <- paste0('LX279_', patient1)

patient.counts <- rep("NA", nrow(srat.lr@meta.data))
names(patient.counts) <- rownames(srat.lr@meta.data)

#Annotation is different from Souporcell's, so maybe we need to swap patients
patient.counts[WhichCells(srat.lr, cells = patient0)] <- "patient0"
patient.counts[WhichCells(srat.lr, cells = patient1)] <- "patient1"

#See where each cell has been classified
srat.lr <- AddMetaData(srat.lr, col.name = "patient.counts", metadata = patient.counts)


#Patient 0 ----
srat0 <- subset(srat.lr, subset = patient_ont == "0")
srat1 <- subset(srat.lr, subset = patient_ont == "1")

bcs <- gsub('LX279_','',colnames(srat0))
pat0.selected.bcs <- nuc.matrix[,colnames(nuc.matrix) %in% bcs]

#Nuc patient 0
pat0.bcs <- gsub('LX279_','',colnames(srat0))
pat0.selected.bcs <- nuc.matrix[,colnames(nuc.matrix) %in% pat0.bcs]

#Nuc patient 0
#Metadata per patient now
pat0.cnts_per_cell <- colSums(pat0.selected.bcs)
pat0.cnts_per_snv <- rowSums(pat0.selected.bcs)
summary(pat0.cnts_per_cell)
summary(pat0.cnts_per_snv)
hist(pat0.cnts_per_cell, breaks = 50)
hist(pat0.cnts_per_snv, breaks = 1000)

#This one is unfiltered, we store its metadata
pat0.nuc.unfiltered <- rowMeans(pat0.selected.bcs, na.rm = TRUE)
pat0.nuc.unfiltered <- as.data.frame(pat0.nuc.unfiltered)
colnames(pat0.nuc.unfiltered) <- 'mean'
pat0.nuc.unfiltered$sd <- rowSds(pat0.selected.bcs, na.rm = TRUE)
pat0.nuc.unfiltered$nonzero <- (apply(pat0.selected.bcs, 1, function(row) sum(row != 0, na.rm = TRUE)) / ncol(pat0.selected.bcs)) * 100

dim(pat0.selected.bcs)
#[1] 2633 3691
summary(colSums(pat0.selected.bcs))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#33     891    1454    2036    2660   13313 
summary(rowSums(pat0.selected.bcs))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0   135.2   376.0  1250.9   720.0 83313.0
pat0.normalized <- LogNormalize(pat0.selected.bcs)

#nuc.normalized stores information of the normalized matrix, namely, info
#after deleting low quality cells.
pat0.nuc.normalized <- rowMeans(pat0.normalized, na.rm = TRUE)
pat0.nuc.normalized <- as.data.frame(pat0.nuc.normalized)
colnames(pat0.nuc.normalized) <- 'mean'
pat0.nuc.normalized$sd <- rowSds(pat0.normalized, na.rm = TRUE)
pat0.nuc.normalized$nonzero <- (apply(pat0.normalized, 1, function(row) sum(row != 0, na.rm = TRUE)) / ncol(pat0.normalized)) * 100

dim(pat0.nuc.normalized)
#[1] 2633 3691
summary(colSums(pat0.normalized))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#33     891    1454    2036    2660   13313 
summary(rowSums(pat0.normalized))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0   135.2   376.0  1250.9   720.0 83313.0

hist(pat0.nuc.normalized$mean, breaks = 100)
hist(pat0.nuc.normalized$sd, breaks = 100)
hist(pat0.nuc.normalized$nonzero, breaks = 100)


#filtered is our matrix after taking only SNVs highly informative, with high SD
pat0.filtered <- pat0.normalized[pat0.nuc.normalized$sd > cutoff$sds,]

##Clustering ----
p0.nuc.srat.meta <- srat.meta[rownames(srat.meta) %in% colnames(pat0.filtered),]
p0.nuc.srat.meta <- srat.meta[colnames(pat0.filtered),]
p0.nuc.srat.annotation <- HeatmapAnnotation(df = p0.nuc.srat.meta, col = list(
  is_aneuploid = c('no' = 'red','yes' = 'blue'),
  patient_ont = c('0' = 'green'),
  cell_type = c('B_cell' = '#FED439FF',
                'Monocyte' = '#709AE1FF',
                'NK_cell' = '#8A9197FF',
                'Other' = '#D2AF81FF',
                'Pro-B_cell_CD34+' = '#FD7446FF',
                'T_cells' = '#D5E4A2FF')
)
)

set.seed(123)
p0.nuc.ht <- ComplexHeatmap::Heatmap(pat0.filtered,
                                     cluster_rows = TRUE,
                                     cluster_columns = TRUE,
                                     use_raster = TRUE,
                                     col = brewer.pal(7, "Greens"),
                                     show_row_names = FALSE,
                                     show_column_names = FALSE,
                                     row_title = 'SNVs',
                                     column_title = 'Barcodes',
                                     column_km = 2,
                                     column_km_repeats = 10,
                                     top_annotation = p0.nuc.srat.annotation)
p0.nuc.ht <- draw(p0.nuc.ht)

p0.healthy <- colnames(pat0.filtered)[column_order(p0.nuc.ht)[[2]]]
p0.tumor <- colnames(pat0.filtered)[column_order(p0.nuc.ht)[[1]]]

## Adding metadata to srat ----
p0.healthy <- paste0('LX279_', p0.healthy)
p0.tumor <- paste0('LX279_', p0.tumor)

test.counts <- rep("NA", nrow(srat.lr@meta.data))
names(test.counts) <- rownames(srat.lr@meta.data)

#In the reference the patients are inversed
test.counts[WhichCells(srat.lr, cells = p0.healthy)] <- "p0.tumor"
test.counts[WhichCells(srat.lr, cells = p0.tumor)] <- "p0.healthy"

#Patient 1 ----
#Subseting the full matrix (with ALL SNVs and HIGH QUAL cells) for cells
#from patient 1 coming from our clustering of SNVs
pat1.bcs <- gsub('LX279_','',colnames(srat1))
pat1.selected.bcs <- nuc.matrix[,colnames(nuc.matrix) %in% pat1.bcs]

#Nuc patient 0
pat1.cnts_per_cell <- colSums(pat1.selected.bcs)
pat1.cnts_per_snv <- rowSums(pat1.selected.bcs)
summary(pat1.cnts_per_cell)
summary(pat1.cnts_per_snv)
hist(pat1.cnts_per_cell, breaks = 50)
hist(pat1.cnts_per_snv, breaks = 1000)

pat1.nuc.unfiltered <- rowMeans(pat1.selected.bcs, na.rm = TRUE)
pat1.nuc.unfiltered <- as.data.frame(pat1.nuc.unfiltered)
colnames(pat1.nuc.unfiltered) <- 'mean'
pat1.nuc.unfiltered$sd <- rowSds(pat1.selected.bcs, na.rm = TRUE)
pat1.nuc.unfiltered$nonzero <- (apply(pat1.selected.bcs, 1, function(row) sum(row != 0, na.rm = TRUE)) / ncol(pat1.selected.bcs)) * 100

dim(pat1.selected.bcs)
#[1] 2633 3691
summary(colSums(pat1.selected.bcs))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#33     891    1454    2036    2660   13313 
summary(rowSums(pat1.selected.bcs))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0   135.2   376.0  1250.9   720.0 83313.0
pat1.normalized <- LogNormalize(pat1.selected.bcs)

#filtered.bcs stores information of the selected.bcs matrix, namely, info
#after deleting low quality cells.
pat1.nuc.normalized <- rowMeans(pat1.normalized, na.rm = TRUE)
pat1.nuc.normalized <- as.data.frame(pat1.nuc.normalized)
colnames(pat1.nuc.normalized) <- 'mean'
pat1.nuc.normalized$sd <- rowSds(pat1.normalized, na.rm = TRUE)
pat1.nuc.normalized$nonzero <- (apply(pat1.normalized, 1, function(row) sum(row != 0, na.rm = TRUE)) / ncol(pat1.normalized)) * 100

dim(pat1.nuc.normalized)
#[1] 2633 3691
summary(colSums(pat1.normalized))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#33     891    1454    2036    2660   13313 
summary(rowSums(pat1.normalized))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0   135.2   376.0  1250.9   720.0 83313.0

hist(pat1.nuc.normalized$mean, breaks = 100)
hist(pat1.nuc.normalized$sd, breaks = 100)
hist(pat1.nuc.normalized$nonzero, breaks = 100)


#selected.snvs is our matrix after taking only high coverage SNVs AFTER
#filtering low qual cells
pat1.filtered <- pat1.normalized[pat1.nuc.normalized$sd > cutoff$sds,]

##Clustering ----
p1.nuc.srat.meta <- srat.meta[rownames(srat.meta) %in% colnames(pat1.filtered),]
p1.nuc.srat.meta <- srat.meta[colnames(pat1.filtered),]
p1.nuc.srat.annotation <- HeatmapAnnotation(df = p1.nuc.srat.meta, col = list(
  is_aneuploid = c('no' = 'green','yes' = 'red'),
    patient = c('1' = 'purple'),
    cell_type = c('B_cell' = '#FED439FF',
                  'Monocyte' = '#709AE1FF',
                  'NK_cell' = '#8A9197FF',
                  'Other' = '#D2AF81FF',
                  'Pro-B_cell_CD34+' = '#FD7446FF',
                  'T_cells' = '#D5E4A2FF')
  )
)

set.seed(123)
p1.nuc.ht <- ComplexHeatmap::Heatmap(pat1.filtered,
                                     cluster_rows = TRUE,
                                     cluster_columns = TRUE,
                                     use_raster = TRUE,
                                     col = brewer.pal(7, "Greens"),
                                     show_row_names = FALSE,
                                     show_column_names = FALSE,
                                     row_title = 'SNVs',
                                     column_title = 'Barcodes',
                                     column_km = 2,
                                     column_km_repeats = 10,
                                     top_annotation = p1.nuc.srat.annotation)
p1.nuc.ht <- draw(p1.nuc.ht)

p1.tumor <- colnames(pat1.filtered)[column_order(p1.nuc.ht)[[2]]]
p1.healthy <- colnames(pat1.filtered)[column_order(p1.nuc.ht)[[1]]]

## Adding metadata to srat ----
p1.healthy <- paste0('LX279_', p1.healthy)
p1.tumor <- paste0('LX279_', p1.tumor)

#In the reference the patients are inversed
test.counts[WhichCells(srat.lr, cells = p1.healthy)] <- "p1.tumor"
test.counts[WhichCells(srat.lr, cells = p1.tumor)] <- "p1.healthy"

srat.lr <- AddMetaData(srat.lr, col.name = "test.counts", metadata = test.counts)

#AUROC ----
table(srat.lr$test.counts, srat.lr$reference)

groundtruth_binary <- ifelse(srat.lr$is_aneuploid == "yes", 1, 0)

srat.lr$test.counts_binary <- ifelse(srat.lr$test.counts %in% c("p0.healthy", "p1.healthy"), 0, 
                              ifelse(srat.lr$test.counts %in% c("p0.tumor", "p1.tumor"), 1, NA))


valid_indices <- !is.na(srat.lr$test.counts_binary)  # Get indices where prediction is not NA
predicted <- srat.lr$test.counts_binary[valid_indices]  # Keep only valid predictions
groundtruth_valid <- groundtruth_binary[valid_indices]  # Keep corresponding ground truth values

# Step 7: Generate ROC Curve and calculate AUC
roc_curve <- roc(groundtruth_valid, predicted)
roc_data <- data.frame(
  specificity = rev(roc_curve$specificities),  # Reverse for correct plotting
  sensitivity = rev(roc_curve$sensitivities)   # Reverse for correct plotting
)

# Step 9: Calculate AUC and display
auc_value <- auc(roc_curve)
# Add AUC to the plot
ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "blue", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = paste("ROC Curve (AUC =", round(auc_value, 2), ")"), 
       x = "False Positive Rate (1 - Specificity)", 
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal()



DimPlot(srat.lr,
        reduction = "umap",
        pt.size = 0.5,
        group.by = "test.counts") + labs(title = 'Genotype and aneuploidy inference ONT SNV counts') + ggsci::scale_color_jco()

DimPlot(
  srat.lr,
  group.by = 'test.counts',
  cols = c('#AC92EB', '#4FC1E8', '#A0D568', '#FFCE54', '#ED5564')
) + labs(title = 'Patient and cell condition annotation using counts - knn')|
  DimPlot(
    srat.lr,
    group.by = 'reference',
    cols = c('#4FC1E8', '#A0D568', '#FFCE54', '#ED5564')
  )+ labs(title = 'Reference (souporcell-Ilmn + InferCNV - ONT')

table(srat.lr$test.counts, srat.lr$reference)



#Creating metacells ----
nuc.srat.annotation <- HeatmapAnnotation(
  Reference = nuc.srat.meta$reference,
  Cluster = nuc.srat.meta$cluster,
  col = list(
    Reference = c("p0.healthy" = "red", "p0.tumor" = "blue", 
                  "p1.healthy" = "green", "p1.tumor" = "yellow"),
    Cluster = structure(
      circlize::rand_color(17),  # Generar 17 colores aleatorios para los clústeres (0-16)
      names = as.character(0:16)  # Asegurarse que los clústeres tengan nombres '0' a '16'
    )
  )
)

# Ordenar las columnas de la matriz según la referencia y el clúster
ord <- order(nuc.srat.meta$reference)
mat_ordered <- pat0.selected.bcs[, ord]

# Crear el heatmap con las columnas ordenadas y anotaciones
Heatmap(mat_ordered, 
        top_annotation = nuc.srat.annotation, 
        show_column_names = FALSE,
        show_row_names = TRUE,
        cluster_columns = FALSE,
        use_raster = TRUE)


# Heatmap without cols clustering ----
#Order cells following reference and clust
srat.meta$bcs <- rownames(srat.meta)
srat.meta.sample <- srat.meta[rownames(srat.meta) %in% gsub('LX279_', '',colnames(dataset)),]
srat.meta.sample <- srat.meta.sample[order(srat.meta.sample$reference, srat.meta.sample$cluster),]


ordered_bcs <- srat.meta.sample$bcs
subset_SCT <- subset_SCT[,colnames(subset_SCT) %in% ordered_bcs]

dim(subset_SCT)
# Genera la anotación para las columnas

nuc.srat.annotation <- HeatmapAnnotation(
  Reference = srat.meta.sample$reference,
  Cluster = srat.meta.sample$cluster,
  col = list(Reference = c('p0.healthy' = '#66c2a5', 
                           'p1.healthy' = '#fc8d62',
                           'p0.tumor' = '#8da0cb',
                           'p1.tumor' = '#e78ac3')),
    Cluster = structure(
      circlize::rand_color(17),  # 
      names = as.character(0:16)  # 
  )
)

limited_data <- pmax(pmin(subset_SCT, 2.5), -2.5)
#
Heatmap(dataset,
        top_annotation = nuc.srat.annotation, 
        cluster_columns = FALSE, # Keep order of matrix
        show_column_names = FALSE,
        cluster_rows = FALSE,
        col = c('black', 'yellow'))



Idents(srat.snvs) <- "SCT_snn_res.1.4"
srat.snvs.small <- subset(srat.snvs, downsample = 50)
# Which genes to plot
genes_union <- union(top3.de.snvs$gene, top10.pat.snvs$gene)
#Get matrix for Heatmap
dataset <- as.matrix(GetAssayData(srat.snvs.small, assay = "RNA", slot = "scale.data")[genes_union, ])
colnames(dataset) <- gsub('LX279_', '', colnames(dataset))

srat.meta$bcs <- rownames(srat.meta)
srat.meta.sample <- srat.meta[rownames(srat.meta) %in% gsub('LX279_', '',colnames(dataset)),]
srat.meta.sample <- srat.meta.sample[order(srat.meta.sample$reference, srat.meta.sample$cluster),]

ordered_dataset <- dataset[, rownames(srat.meta.sample)]

nuc.srat.annotation <- HeatmapAnnotation(
  Reference = srat.meta.sample$reference,
  Cluster = srat.meta.sample$cluster,
  col = list(Reference = c('p0.healthy' = '#66c2a5', 
                           'p1.healthy' = '#fc8d62',
                           'p0.tumor' = '#8da0cb',
                           'p1.tumor' = '#e78ac3')),
  Cluster = structure(
    circlize::rand_color(17),  # 
    names = as.character(0:16)  # 
  )
)

final.hm <- Heatmap(ordered_dataset,
        top_annotation = nuc.srat.annotation, 
        cluster_columns = FALSE, # Keep order of matrix
        show_column_names = FALSE,
        cluster_rows = FALSE,
        show_row_names = FALSE,
        col = c('black', 'yellow'),
        column_split = srat.meta.sample$reference,
        column_gap = unit(0.5, "mm"))

draw(final.hm)

