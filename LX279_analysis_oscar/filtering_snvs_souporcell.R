remove.packages('Matrix')
devtools::install_version("Matrix", version = "1.5.0", repos = "http://cran.us.r-project.org")

detach("package:Matrix", unload = TRUE)

local.lib <- '/hpc/pmc_holstege/rstudio/oscar_clonetracer/R/x86_64-pc-linux-gnu-library/4.4'

library(Seurat)
library(Matrix)
library(ComplexHeatmap)
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

#Loading data ----
#vafs.merged<-read.table(file = '~/snv_calling/processing.vafs/vafs.merged.tsv', header = TRUE)
vafs.per.cell<-read.table(file = '~/vaf.per.cell.tsv', header = TRUE)
cutoff <- list()
cutoff$QUAL <- 30 #Higher than 30
cutoff$NAs <- 80 #Maximum of n% of NAs, in other words, highly covered
cutoff$cells <- 95 #Cells can have a maximum of n% of NAs for the total of SNVs
cutoff$sds <- 0.3 #SNVs selected must have more than n% SD
cutoff$cumulative.var <- 80 #For selection of PCs, amount of cumulative Var needed 

#String manipulation ----
vafs.per.cell$ID <- NULL
vafs.per.cell<- vafs.per.cell %>% unite(col='SNV',CHROM:ALT, sep=':')
rownames(vafs.per.cell) <- vafs.per.cell$SNV
vafs.per.cell$SNV <- NULL
plot(density(vafs.per.cell$QUAL), main = 'Density distribution of SNV quality', xlab = 'QUAL')

vafs.per.cell <- vafs.per.cell[vafs.per.cell$QUAL > cutoff$QUAL,]
vafs.per.cell$QUAL <- NULL
vafs.per.cell <- as.matrix(vafs.per.cell)

#Split by nuc and mito SNVs
mitochondrial.snvs <- grepl("^chrM", rownames(vafs.per.cell))
mito.matrix <- vafs.per.cell[mitochondrial.snvs, ]
nuc.matrix <- vafs.per.cell[!mitochondrial.snvs, ]

#NUCLEAR SNVS  ----x
##Storing metadata ----
nuc.unfiltered <- rowMeans(nuc.matrix, na.rm = TRUE)
nuc.unfiltered <- as.data.frame(nuc.unfiltered)
colnames(nuc.unfiltered) <- 'mean'
nuc.unfiltered$sd <- rowSds(nuc.matrix, na.rm = TRUE)

informative_value <- function(x) {
  count_zero <- sum(x == 0, na.rm = TRUE)
  count_na <- sum(is.na(x))
  return((ncol(nuc.matrix)-(count_zero + count_na))/ncol(nuc.matrix))
}

nuc.unfiltered$nonzero <- (apply(nuc.matrix, 1, informative_value))*100
nuc.unfiltered$NAs <- ((apply(nuc.matrix, 1, function(x) sum(is.na(x))))/ncol(nuc.matrix))*100

rownames(nuc.matrix[nuc.unfiltered$NAs<cutoff$NAs,])
# [1] "chr11:65502026:G:A" "chr11:65502633:A:G" "chr16:89561052:T:G" "chr16:89561263:C:T" "chr16:89561665:G:A" "chr17:41689883:C:T"
#[7] "chr17:75779744:G:A" "chr6:31353908:A:C"  "chr6:31353975:T:G"  "chr6:31354030:G:A"  "chr6:31354079:A:G"  "chr6:31354105:G:A" 
#[13] "chr6:31354129:A:G"  "chr6:31354138:A:G"  "chr6:31354139:A:G"  "chr6:31354181:G:A"  "chr6:31354249:C:A"  "chr6:31355111:A:G" 
#[19] "chr6:31355134:C:T"  "chr6:31355203:C:T"  "chr6:31355219:C:T"  "chr6:31355456:A:G"  "chr6:31355519:A:G"  "chr6:31355544:G:A" 
#[25] "chr6:31355560:T:C"  "chr6:31355576:G:A"  "chrX:72273841:C:T"  "chrX:72275699:T:C" 

p_nuc_unf_mean <- ggplot(nuc.unfiltered, aes(x = mean)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  geom_density() + 
  labs(title = "Histogram of Mean Values",
       x = "Mean",
       y = "Frequency")
p_nuc_unf_sd <- ggplot(nuc.unfiltered, aes(x = sd)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  geom_density() + 
  labs(title = "Histogram of Standard Deviations",
       x = "SD",
       y = "Frequency")
p_nuc_unf_nonzero <- ggplot(nuc.unfiltered, aes(x = nonzero)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  geom_density() + 
  labs(title = "Histogram of Non-zero values",
       x = "Percentage of Non-zeroes",
       y = "Frequency")
p_nuc_unf_na <- ggplot(nuc.unfiltered, aes(x = NAs)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  geom_density() + 
  labs(title = "Histogram of NAs",
       x = "Percentage of NAs",
       y = "Frequency")
(sum(is.na(nuc.matrix))/(ncol(nuc.matrix)*nrow(nuc.matrix)))*100
#[1] 81.06532 percentage of NAs in our unprocessed data

##Filtering ----
#We have SNVs and barcodes. We first filter barcodes from transcriptomic
#analysis made previously, only taking high quality cells. We also store
#metadata/annotation of these cells for further analysis.
srat.lr <- readRDS(file = "~/nanopore_souporcell/srat-aneuploidy.rds")
srat.meta <- data.frame(row.names = gsub('LX279_','',colnames(srat.lr)),
                        patient = srat.lr@meta.data$patient,
                        is_aneuploid = srat.lr@meta.data$is_aneuploid,
                        cell_type = srat.lr@meta.data$singler)
rownames(srat.meta) <- gsub('-1','', rownames(srat.meta))
length(intersect(rownames(srat.meta), colnames(nuc.matrix)))
#[1] 3691
filtered.barcodes <- rownames(srat.meta)
nuc.selected.cells <- nuc.matrix[, colnames(nuc.matrix) %in% rownames(srat.meta)]
dim(nuc.selected.cells)
#[1] 2633 3691

#filtered.bcs stores information of the selected.bcs matrix, namely, info
#after deleting low quality cells.
nuc.filtered.bcs <- rowMeans(nuc.selected.cells, na.rm = TRUE)
nuc.filtered.bcs <- as.data.frame(nuc.filtered.bcs)
colnames(nuc.filtered.bcs) <- 'mean'
nuc.filtered.bcs$sd <- rowSds(nuc.selected.cells, na.rm = TRUE)
nuc.filtered.bcs$nonzero <- (apply(nuc.selected.cells, 1, informative_value))*100
nuc.filtered.bcs$NAs <- ((apply(nuc.selected.cells, 1, function(x) sum(is.na(x))))/ncol(nuc.selected.cells))*100
(sum(is.na(nuc.selected.cells))/(ncol(nuc.selected.cells)*nrow(nuc.selected.cells)))*100
#[1] 52.36435 percentage of NAs after low qual cell filtering, not a huge change

hist(nuc.filtered.bcs$mean, breaks = 100)
hist(nuc.filtered.bcs$sd, breaks = 100)
hist(nuc.filtered.bcs$nonzero, breaks = 100)
hist(nuc.filtered.bcs$NAs, breaks = 100)


#selected.snvs is our matrix after taking only high coverage SNVs AFTER
#filtering low qual cells
nuc.selected.snvs <- nuc.selected.cells[nuc.filtered.bcs$NAs<cutoff$NAs & nuc.filtered.bcs$sd>cutoff$sds,]

#filtered.snvs object is once again metadata from this matrix, so only high
#quality cells AND SNVs with less than 20% NAs
nuc.filtered.snvs <- rowMeans(nuc.selected.snvs, na.rm = TRUE)
nuc.filtered.snvs <- as.data.frame(nuc.filtered.snvs)
colnames(nuc.filtered.snvs) <- 'mean'
nuc.filtered.snvs$sd <- rowSds(nuc.selected.snvs, na.rm = TRUE)
nuc.filtered.snvs$nonzero <- (apply(nuc.selected.snvs, 1, informative_value))*100
nuc.filtered.snvs$NAs <- ((apply(nuc.selected.snvs, 1, function(x) sum(is.na(x))))/ncol(nuc.selected.snvs))*100

p_nuc_filt_mean <- ggplot(nuc.filtered.snvs, aes(x = mean)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  geom_density() +
  labs(title = "Histogram of Mean Values",
       x = "Mean",
       y = "Frequency")
p_nuc_filt_sd <- ggplot(nuc.filtered.snvs, aes(x = sd)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  geom_density() + 
  labs(title = "Histogram of Standard Deviations",
       x = "SD",
       y = "Frequency")
p_nuc_filt_nonzero <- ggplot(nuc.filtered.snvs, aes(x = nonzero)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  geom_density() + 
  labs(title = "Histogram of Non-zero values",
       x = "Percentage of Non-zeroes",
       y = "Frequency")
p_nuc_filt_na <- ggplot(nuc.filtered.snvs, aes(x = NAs)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  geom_density() + 
  labs(title = "Histogram of NAs",
       x = "Percentage of NAs",
       y = "Frequency")

plot_grid(p_nuc_unf_mean, p_nuc_filt_mean)
plot_grid(p_nuc_unf_sd, p_nuc_filt_sd)
plot_grid(p_nuc_unf_nonzero, p_nuc_filt_nonzero)
plot_grid(p_nuc_unf_na, p_nuc_filt_na)

#Now we look only at cells, not SNVs. We look at columns with high coverage,
#namely, cells that have at least 1 read for more than 60% of the SNVs,
#in other words, less than 40% of non-covered SNVs
nuc.metadata.cells <- ((apply(nuc.selected.snvs, 2, function(x) sum(is.na(x))))/nrow(nuc.selected.snvs))*100
nuc.metadata.cells <- as.data.frame(nuc.metadata.cells)
colnames(nuc.metadata.cells) <- 'NAs'
table(nuc.metadata.cells$NAs < cutoff$cells)
#FALSE  TRUE 
#114  3270

nuc.vafs.filtered <- nuc.selected.snvs[, nuc.metadata.cells$NAs < cutoff$cells]
dim(nuc.vafs.filtered)
#[1]   65 3270

#Percentage of NAs once filtered
(sum(is.na(nuc.vafs.filtered))/(ncol(nuc.vafs.filtered)*nrow(nuc.vafs.filtered)))*100
#[1] 9.670014

##Visualisation of NAs ----
nuc_missing_proportion <- colSums(is.na(nuc.vafs.filtered)) / nrow(nuc.vafs.filtered)
summary(nuc_missing_proportion)
nuc_na_columns <- sum(nuc_missing_proportion > 0)
nuc_na_columns
ggplot(data = melt(is.na(nuc.vafs.filtered)), aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "black"), name = "NA") +
  theme_minimal() +
  labs(title = "Heatmap of Missing Values", x = "Columns", y = "Rows")

##Imputation ----
nuc.vafs.imputed <- impute.knn(nuc.vafs.filtered, colmax = 85)
hist(rowSds(nuc.vafs.imputed$data), breaks = 100)
nuc.srat.meta <- srat.meta[rownames(srat.meta) %in% colnames(nuc.vafs.filtered),]
nuc.srat.meta <- srat.meta[colnames(nuc.vafs.filtered),]
nuc.srat.annotation <- HeatmapAnnotation(df = nuc.srat.meta)

nuc.ht <- ComplexHeatmap::Heatmap(nuc.vafs.imputed$data,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        use_raster = FALSE,
        na_col = '#808080',
        col = brewer.pal(7, "Greens"),
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_title = 'SNVs',
        column_title = 'Barcodes',
        column_km = 2,
        column_km_repeats = 10,
        top_annotation = nuc.srat.annotation)
nuc.ht <- draw(nuc.ht)
dim(nuc.vafs.imputed$data)
column_order(nuc.ht)

##Comparison of annotation souporcell - knn ----
knn_annotation <- rep("NA", nrow(srat.lr@meta.data))
names(knn_annotation) <- rownames(srat.lr@meta.data)

knn_annotation[WhichCells(srat.lr, cells = paste0('LX279_', colnames(nuc.vafs.filtered)[column_order(nuc.ht)[[1]]]))] <- "1"
knn_annotation[WhichCells(srat.lr, cells = paste0('LX279_', colnames(nuc.vafs.filtered)[column_order(nuc.ht)[[2]]]))] <- "0"

srat.lr <- AddMetaData(srat.lr, col.name = "knn_patient", metadata = knn_annotation)

DimPlot(srat.lr, group.by = 'patient', cols = c('#a6e0e9', '#f64941')) + labs(title = 'Annotation using souporcell') |
  DimPlot(srat.lr, group.by = 'knn_patient', cols = c('#a6e0e9', '#f64941', '#222224')) + labs(title = 'Annotation using knn on ONT')

#MITOCHONDRIAL SNVS  ----
##Storing metadata ----
mito.unfiltered <- rowMeans(mito.matrix, na.rm = TRUE)
mito.unfiltered <- as.data.frame(mito.unfiltered)
colnames(mito.unfiltered) <- 'mean'
mito.unfiltered$sd <- rowSds(mito.matrix, na.rm = TRUE)

informative_value <- function(x) {
  count_zero <- sum(x == 0, na.rm = TRUE)
  count_na <- sum(is.na(x))
  return((ncol(mito.matrix)-(count_zero + count_na))/ncol(mito.matrix))
}

mito.unfiltered$nonzero <- (apply(mito.matrix, 1, informative_value))*100
mito.unfiltered$NAs <- ((apply(mito.matrix, 1, function(x) sum(is.na(x))))/ncol(mito.matrix))*100
rownames(mito.matrix[mito.unfiltered$NAs<20,])
#[1] "chrM:663:A:G"   "chrM:750:A:G"   "chrM:1438:A:G"  "chrM:1736:A:G"  "chrM:2706:A:G"  "chrM:2838:A:C"  "chrM:3197:T:C" 
#[8] "chrM:3801:T:C"  "chrM:4248:T:C"  "chrM:4769:A:G"  "chrM:4824:A:G"  "chrM:6488:T:C"  "chrM:7028:C:T"  "chrM:7124:A:G" 
#[15] "chrM:8027:G:A"  "chrM:8158:A:G"  "chrM:8794:C:T"  "chrM:8860:A:G"  "chrM:9477:G:A"  "chrM:11016:G:A" "chrM:11326:C:T"
#[22] "chrM:11467:A:G" "chrM:11719:G:A" "chrM:12007:G:A" "chrM:12101:T:C" "chrM:12705:C:T" "chrM:14766:C:T" "chrM:15232:A:G"
#[29] "chrM:15326:A:G" 

p_mito_unf_mean <- ggplot(mito.unfiltered, aes(x = mean)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Mean Values",
       x = "Mean",
       y = "Frequency")
p_mito_unf_sd <- ggplot(mito.unfiltered, aes(x = sd)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Standard Deviations",
       x = "SD",
       y = "Frequency")
p_mito_unf_nonzero <- ggplot(mito.unfiltered, aes(x = nonzero)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Non-zero values",
       x = "Percentage of Non-zeroes",
       y = "Frequency")
p_mito_unf_na <- ggplot(mito.unfiltered, aes(x = NAs)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of NAs",
       x = "Percentage of NAs",
       y = "Frequency")
(sum(is.na(mito.matrix))/(ncol(mito.matrix)*nrow(mito.matrix)))*100
#[1] 12.01474 percentage of NAs in our unprocessed data. Much lower than nuclear

##Filtering ----
#We have SNVs and barcodes. We first filter barcodes from transcriptomic
#analysis made previously, only taking high quality cells. We also store
#metadata/annotation of these cells for further analysis.
mito.selected.cells <- mito.matrix[, colnames(mito.matrix) %in% rownames(srat.meta)]
dim(mito.selected.cells)
#[1]   34 3691

#filtered.bcs stores information of the selected.bcs matrix, namely, info
#after deleting low quality cells.
mito.filtered.bcs <- rowMeans(mito.selected.cells, na.rm = TRUE)
mito.filtered.bcs <- as.data.frame(mito.filtered.bcs)
colnames(mito.filtered.bcs) <- 'mean'
mito.filtered.bcs$sd <- rowSds(mito.selected.cells, na.rm = TRUE)
mito.filtered.bcs$nonzero <- (apply(mito.selected.cells, 1, informative_value))*100
mito.filtered.bcs$NAs <- ((apply(mito.selected.cells, 1, function(x) sum(is.na(x))))/ncol(mito.selected.cells))*100
(sum(is.na(mito.selected.cells))/(ncol(mito.selected.cells)*nrow(mito.selected.cells)))*100
#[1] 7.445774 percentage of NAs after low qual cell filtering, decent coming form 12%

plot(density(mito.filtered.bcs$mean))
plot(density(mito.filtered.bcs$sd, breaks = 100))
plot(density(mito.filtered.bcs$nonzero, breaks = 100))
plot(density(mito.filtered.bcs$NAs, breaks = 100))


#selected.snvs is our matrix after taking only high coverage SNVs AFTER
#filtering low qual cells
mito.selected.snvs <- mito.selected.cells[mito.filtered.bcs$NAs<20,]
#filtered.snvs object is once again metadata from this matrix, so only high
#quality cells AND SNVs with less than 20% NAs
mito.filtered.snvs <- rowMeans(mito.selected.snvs, na.rm = TRUE)
mito.filtered.snvs <- as.data.frame(mito.filtered.snvs)
colnames(mito.filtered.snvs) <- 'mean'
mito.filtered.snvs$sd <- rowSds(mito.selected.snvs, na.rm = TRUE)
mito.filtered.snvs$nonzero <- (apply(mito.selected.snvs, 1, informative_value))*100
mito.filtered.snvs$NAs <- ((apply(mito.selected.snvs, 1, function(x) sum(is.na(x))))/ncol(mito.selected.snvs))*100

p_mito_filt_mean <- ggplot(mito.filtered.snvs, aes(x = mean)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Mean Values",
       x = "Mean",
       y = "Frequency")
p_mito_filt_sd <- ggplot(mito.filtered.snvs, aes(x = sd)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Standard Deviations",
       x = "SD",
       y = "Frequency")
p_mito_filt_nonzero <- ggplot(mito.filtered.snvs, aes(x = nonzero)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Non-zero values",
       x = "Percentage of Non-zeroes",
       y = "Frequency")
p_mito_filt_na <- ggplot(mito.filtered.snvs, aes(x = NAs)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of NAs",
       x = "Percentage of NAs",
       y = "Frequency")

plot_grid(p_mito_unf_mean, p_mito_filt_mean)
plot_grid(p_mito_unf_sd, p_mito_filt_sd)
plot_grid(p_mito_unf_nonzero, p_mito_filt_nonzero)
plot_grid(p_mito_unf_na, p_mito_filt_na)

#Now we look only at cells, not SNVs. We look at columns with high coverage,
#namely, cells that have at least 1 read for more than 60% of the SNVs,
#in other words, less than 40% of non-covered SNVs
mito.metadata.cells <- ((apply(mito.selected.snvs, 2, function(x) sum(is.na(x))))/nrow(mito.selected.snvs))*100
mito.metadata.cells <- as.data.frame(mito.metadata.cells)
colnames(mito.metadata.cells) <- 'NAs'
table(mito.metadata.cells$NAs < 40)
#FALSE  TRUE 
#27  3357

mito.vafs.filtered <- mito.selected.snvs[, mito.metadata.cells$NAs<40]
dim(mito.vafs.filtered)
#[1]   30 3640

#Percentage of NAs once filtered
(sum(is.na(mito.vafs.filtered))/(ncol(mito.vafs.filtered)*nrow(mito.vafs.filtered)))*100
#[1] 2.234432, much lower from 80%

##Visualisation of NAs ----
mito_missing_proportion <- colSums(is.na(mito.vafs.filtered)) / nrow(mito.vafs.filtered)
summary(mito_missing_proportion)
mito_na_columns <- sum(mito_missing_proportion > 0)
mito_na_columns
ggplot(data = melt(is.na(mito.vafs.filtered)), aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "black"), name = "NA") +
  theme_minimal() +
  labs(title = "Heatmap of Missing Values", x = "Columns", y = "Rows")

##Imputation ----
mito.vafs.imputed <- impute.knn(mito.vafs.filtered)
hist(rowSds(mito.vafs.imputed$data), breaks = 100)
mito.srat.meta <- srat.meta[rownames(srat.meta) %in% colnames(mito.vafs.filtered),]
mito.srat.meta <- srat.meta[colnames(mito.vafs.filtered),]
mito.srat.annotation <- HeatmapAnnotation(df = mito.srat.meta)

mito.ht <- ComplexHeatmap::Heatmap(mito.vafs.imputed$data,
                                  cluster_rows = TRUE,
                                  cluster_columns = TRUE,
                                  use_raster = FALSE,
                                  na_col = '#808080',
                                  col = brewer.pal(7, "Greens"),
                                  show_row_names = TRUE,
                                  show_column_names = FALSE,
                                  row_title = 'SNVs',
                                  column_title = 'Barcodes',
                                  column_km = 2,
                                  column_km_repeats = 10,
                                  top_annotation = mito.srat.annotation)
mito.ht <- draw(mito.ht)
column_order(mito.ht)
#It's a bit confusing but ht[[2]] corresponds to our pat1. To keep nomenclature
#as up until now, we just swap them.
pat1 <- column_order(nuc.ht)[[2]]
pat2 <- column_order(nuc.ht)[[1]]

#PCA ----
#PCA has to be done AFTER imputation because it cannot handle NA values
nuc.pca <- prcomp(t(nuc.vafs.imputed$data))
pc1 <- nuc.pca$rotation[, 1]
pc2 <- nuc.pca$rotation[, 2]
pc1.contributors <- rownames(nuc.vafs.imputed$data)[order(abs(pc1), decreasing = TRUE)]
pc1.contributors <- pc1.contributors[!is.na(pc1.contributors)]
pc1.contributors

pc2.contributors <- rownames(nuc.vafs.imputed$data)[order(abs(pc2), decreasing = TRUE)]
pc2.contributors <- pc2.contributors[!is.na(pc2.contributors)]
pc2.contributors

var.explained <- nuc.pca$sdev^2
var.explained.percent <- var.explained / sum(var.explained) * 100
cumulative.var <- cumsum(var.explained.percent)
num.pcs <- which(cumulative.var >= cutoff$cumulative.var)[1]
num.pcs <- max(num.pcs, 2)

df <- data.frame(PC = 1:length(var.explained.percent),
                 Variance.Explained = var.explained.percent)

ggplot(df, aes(x = PC, y = Variance.Explained)) +
  geom_point() +
  geom_line() +
  xlab("Principal Components") +
  ylab("Var explained (%)") +
  ggtitle("Elbow Plot")

## Projection of patient data on 2D ----
vecp1 <- srat.meta[colnames(nuc.vafs.imputed$data),'patient']
vecp1 <- as.numeric(vecp1)
vecp1[is.na(vecp1)] <- -1 

nuc.pc.data <- as.data.frame(nuc.pca$x)
nuc.pc.data$color <- vecp1 + 1

ggplot(nuc.pc.data, aes(x = PC1, y = PC2, color = factor(color))) +
  geom_point(shape = 20) +
  labs(x = "PC1", y = "PC2", color = "Genotype")


nuc.rotation.matrix <- nuc.pca$rotation

# extracting top snvs from each pc
get_top_snvs <- function(x, pc, top_n = 10) {
  pc_loadings <- x[, pc]
  
  sorted_indices <- order(abs(pc_loadings), decreasing = TRUE)
  
  top_snvs <- rownames(x)[sorted_indices][1:top_n]
  top_loadings <- pc_loadings[sorted_indices][1:top_n]
  
  data.frame(Loading = top_loadings)
}

top_snvs_pc1 <- get_top_snvs(nuc.rotation.matrix, "PC1", top_n = 12)
top_snvs_pc1


nuc.pc.data <- nuc.pc.data[,1:ncol(nuc.pc.data)-1]
nuc.pc.data <- t(as.matrix(nuc.pc.data))
nuc.pc.data <- nuc.pc.data[1:num.pcs,]
nuc.pcs.ht <- ComplexHeatmap::Heatmap(nuc.pc.data,
                                  cluster_rows = TRUE,
                                  cluster_columns = TRUE,
                                  use_raster = FALSE,
                                  na_col = '#808080',
                                  col = brewer.pal(7, "Greens"),
                                  show_row_names = TRUE,
                                  show_column_names = FALSE,
                                  row_title = 'SNVs',
                                  column_title = 'Barcodes',
                                  column_km = 2,
                                  column_km_repeats = 10,
                                  top_annotation = nuc.srat.annotation)
nuc.pcs.ht <- draw(nuc.pcs.ht)
#PCA does not help in demultiplexing. It actually looks a bit worse
column_order(nuc.pcs.ht)


#Patient 0 ----

#Subseting the full matrix (with ALL SNVs and HIGH QUAL cells) for cells
#from patient 0 coming from our clustering of nuclear SNVs
#For some reason, patient number (arbitrary) are swapped. We will follow the annotations
#that we have used so far, so our 'cluster 2' will be patient 0
pat0.selected.bcs <- colnames(nuc.vafs.filtered)[column_order(nuc.ht)[[1]]]
pat0.selected.bcs <- nuc.matrix[,colnames(nuc.matrix) %in% pat0.selected.bcs]

#Nuc patient 0
pat0.filtered.bcs <- rowMeans(pat0.selected.bcs, na.rm = TRUE)
pat0.filtered.bcs <- as.data.frame(pat0.filtered.bcs)
colnames(pat0.filtered.bcs) <- 'mean'
pat0.filtered.bcs$sd <- rowSds(pat0.selected.bcs, na.rm = TRUE)
pat0.filtered.bcs$nonzero <- (apply(pat0.selected.bcs, 1, informative_value))*100
pat0.filtered.bcs$NAs <- ((apply(pat0.selected.bcs, 1, function(x) sum(is.na(x))))/ncol(nuc.matrix))*100
(sum(is.na(pat0.selected.bcs))/(ncol(pat0.selected.bcs)*nrow(pat0.selected.bcs)))*100
#[1] 72.78245 percentage of NAs after low qual cell filtering, not a huge change

hist(pat0.filtered.bcs$mean, breaks = 100)
hist(pat0.filtered.bcs$sd, breaks = 100)
hist(pat0.filtered.bcs$nonzero, breaks = 100)
hist(pat0.filtered.bcs$NAs, breaks = 100)


#pat1.selected.snvs is our matrix after taking only high coverage SNVs AFTER
#filtering low qual cells
pat0.selected.snvs <- pat0.selected.bcs[pat0.filtered.bcs$NAs<cutoff$NAs & pat0.filtered.bcs$sd>cutoff$sds,]
#pat1.filtered.snvs object is once again metadata from this matrix, so only high
#quality cells AND SNVs with less than 20% NAs
pat0.filtered.snvs <- rowMeans(pat0.selected.snvs, na.rm = TRUE)
pat0.filtered.snvs <- as.data.frame(pat0.filtered.snvs)
colnames(pat0.filtered.snvs) <- 'mean'
pat0.filtered.snvs$sd <- rowSds(pat0.selected.snvs, na.rm = TRUE)
pat0.filtered.snvs$nonzero <- (apply(pat0.selected.snvs, 1, informative_value))*100
pat0.filtered.snvs$NAs <- ((apply(pat0.selected.snvs, 1, function(x) sum(is.na(x))))/ncol(pat0.selected.snvs))*100

p0_filt_mean <- ggplot(pat0.filtered.snvs, aes(x = mean)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Mean Values",
       x = "Mean",
       y = "Frequency")
p0_filt_sd <- ggplot(pat0.filtered.snvs, aes(x = sd)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Standard Deviations",
       x = "SD",
       y = "Frequency")
p0_filt_nonzero <- ggplot(pat0.filtered.snvs, aes(x = nonzero)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Non-zero values",
       x = "Percentage of Non-zeroes",
       y = "Frequency")
p0_filt_na <- ggplot(pat0.filtered.snvs, aes(x = NAs)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of NAs",
       x = "Percentage of NAs",
       y = "Frequency")

#Now we look only at cells, not SNVs. We look at columns with high coverage,
#namely, cells that have at least 1 read for more than 60% of the SNVs,
#in other words, less than 40% of non-covered SNVs
pat0.metadata.cells <- ((apply(pat0.selected.snvs, 2, function(x) sum(is.na(x))))/nrow(pat0.selected.snvs))*100
pat0.metadata.cells <- as.data.frame(pat0.metadata.cells)
colnames(pat0.metadata.cells) <- 'NAs'
table(pat0.metadata.cells$NAs < 80)
#Higher threshold, otherwise too many cells are removed
#FALSE  TRUE 
#226  1227
pat0.vafs.filtered <- pat0.selected.snvs[, pat0.metadata.cells$NAs<cutoff$cells]
dim(pat0.vafs.filtered)
#[1]   233 1227
#Percentage of NAs once filtered
(sum(is.na(pat0.vafs.filtered))/(ncol(pat0.vafs.filtered)*nrow(pat0.vafs.filtered)))*100
#[1] 40.68753, pretty high in comparison with other patient

# Trick from Thanasis
# Convert from alternative allele to  uncommon allele ratio-----
#summary(rowMeans(pat0.vafs.filtered, na.rm = TRUE))
#uncommon <- rowMeans(pat0.vafs.filtered, na.rm = TRUE) > 0.5

#UOf <- pat0.vafs.filtered
#UOf[uncommon, ] <- 1 - pat0.vafs.filtered[uncommon, ]

#plot(density(as.vector(pat0.vafs.filtered), na.rm = TRUE))
#lines(density(as.vector(UOf), na.rm = TRUE), col = "red")

#pat0.vafs.imputed <- impute.knn(UOf)

pat0.vafs.imputed <- impute.knn(pat0.vafs.filtered)

hist(rowSds(pat0.vafs.imputed$data), breaks = 100)
pat0.srat.meta <- srat.meta[rownames(srat.meta) %in% colnames(pat0.vafs.imputed$data),]
pat0.srat.meta <- pat0.srat.meta[colnames(pat0.vafs.filtered),]
pat0.srat.annotation <- HeatmapAnnotation(df = pat0.srat.meta)


## PCA patient 0 ----
pcapat0 <- prcomp(t(pat0.vafs.imputed$data))
pc1 <- pcapat0$rotation[, 1]
pc2 <- pcapat0$rotation[, 2]
pat0.pc1.contributors <- rownames(pat0.vafs.imputed$data)[order(abs(pc1), decreasing = TRUE)]
pat0.pc1.contributors <- pat0.pc1.contributors[!is.na(pat0.pc1.contributors)]
pat0.pc1.contributors

pat0.pc2.contributors <- rownames(pat0.vafs.imputed$data)[order(abs(pc2), decreasing = TRUE)]
pat0.pc2.contributors <- pat0.pc2.contributors[!is.na(pat0.pc2.contributors)]
pat0.pc2.contributors

pat0.var.explained <- pcapat0$sdev^2
pat0.var.explained.percent <- pat0.var.explained / sum(pat0.var.explained) * 100
pat0.cumulative.var <- cumsum(pat0.var.explained.percent)
pat0.num.pcs <- which(pat0.cumulative.var >= cutoff$cumulative.var)[1]

pat0.df <- data.frame(PC = 1:length(pat0.var.explained.percent),
                 Variance.Explained = pat0.var.explained.percent)

ggplot(pat0.df, aes(x = PC, y = Variance.Explained)) +
  geom_point() +
  geom_line() +
  xlab("Principal Components") +
  ylab("Var explained (%)") +
  ggtitle("Elbow Plot")

## Projection of aneuploidy data on 2D
vecp1 <- srat.meta[colnames(pat0.vafs.imputed$data),'is_aneuploid']
pat0.pca.data <- as.data.frame(pcapat0$x)
pat0.pca.data$aneuploid <- vecp1

ggplot(pat0.pca.data, aes(x = PC1, y = PC2, color = aneuploid)) +
  geom_point() +
  scale_color_manual(values = c("yes" = "blue", "no" = "red", "NA" = "grey")) +
  labs(title = "PC1 vs PC2 by aneuploidity",
       x = "PC1",
       y = "PC2",
       color = "Aneuploid")

umap <- umap(t(pat0.vafs.imputed$data), n_neighbors = 15)
colnames(umap$layout) <- c('UMAP1', 'UMAP2')
test.umap <- as.data.frame(umap$layout)
test.umap$aneuploid <- vecp1

ggplot(test.umap, aes(x = UMAP1, y = UMAP2, color = aneuploid)) +
  geom_point() +
  scale_color_manual(values = c("yes" = "blue", "no" = "red", "NA" = "grey")) +
  labs(title = "UMAP1 vs UMAP2 by aneuploidity",
       x = "UMAP1",
       y = "UMAP2",
       color = "Aneuploid")

pat0.rotation.matrix <- pcapat0$rotation

# extracting top snvs from each pc
get_top_snvs <- function(x, pc, top_n = 10) {
  pc_loadings <- x[, pc]
  
  sorted_indices <- order(abs(pc_loadings), decreasing = TRUE)
  
  top_snvs <- rownames(x)[sorted_indices][1:top_n]
  top_loadings <- pc_loadings[sorted_indices][1:top_n]
  
  data.frame(Loading = top_loadings)
}

pat0.top_snvs_pc1 <- get_top_snvs(pat0.rotation.matrix, "PC1", top_n = 12)
pat0.top_snvs_pc1

#eliminate last column, because it has info used for the plot
pat0.pca.data <- pat0.pca.data[,1:ncol(pat0.pca.data)-1]

pat0.pca.data <- t(as.matrix(pat0.pca.data))
pat0.pca.data <- pat0.pca.data[1:pat0.num.pcs,]
pat0.nuc.pcs.ht <- ComplexHeatmap::Heatmap(pat0.pca.data,
                                      cluster_rows = TRUE,
                                      cluster_columns = TRUE,
                                      use_raster = FALSE,
                                      na_col = '#808080',
                                      col = brewer.pal(7, "Greens"),
                                      show_row_names = TRUE,
                                      show_column_names = FALSE,
                                      row_title = 'SNVs',
                                      column_title = 'Barcodes',
                                      column_km = 2,
                                      column_km_repeats = 10,
                                      top_annotation = pat0.srat.annotation)
pat0.nuc.pcs.ht <- draw(pat0.nuc.pcs.ht)
column_order(pat0.nuc.pcs.ht)



p0.ht <- ComplexHeatmap::Heatmap(pat0.vafs.imputed$data,
                        cluster_rows = TRUE,
                        cluster_columns = TRUE,
                        use_raster = FALSE,
                        na_col = '#808080',
                        col = brewer.pal(7, "Greens"),
                        show_row_names = TRUE,
                        show_column_names = FALSE,
                        row_title = 'SNVs',
                        column_title = 'Barcodes Patient 1',
                        column_km = 2,
                        column_km_repeats = 10,
                        top_annotation = pat0.srat.annotation)
p0.ht <- draw(p0.ht)
dim(pat0.vafs.imputed$data)

p0.healthy <- colnames(pat0.vafs.filtered)[column_order(pat0.nuc.pcs.ht)[[2]]]
p0.tumor <- colnames(pat0.vafs.filtered)[column_order(pat0.nuc.pcs.ht)[[1]]]

## Adding metadata to srat ----
p0.healthy <- paste0('LX279_', p0.healthy)
p0.tumor <- paste0('LX279_', p0.tumor)

test <- rep("NA", nrow(srat.lr@meta.data))
names(test) <- rownames(srat.lr@meta.data)

#In the reference the patients are inversed
test[WhichCells(srat.lr, cells = p0.healthy)] <- "p1.tumor"
test[WhichCells(srat.lr, cells = p0.tumor)] <- "p1.healthy"



#Patient 1 ----
#Subseting the full matrix (with ALL SNVs and HIGH QUAL cells) for cells
#from patient 1 coming from our clustering of SNVs
pat1.selected.bcs <- colnames(nuc.vafs.filtered)[column_order(nuc.ht)[[2]]]
pat1.selected.bcs <- nuc.matrix[,colnames(nuc.matrix) %in% pat1.selected.bcs]

#Nuc patient 1
pat1.filtered.bcs <- rowMeans(pat1.selected.bcs, na.rm = TRUE)
pat1.filtered.bcs <- as.data.frame(pat1.filtered.bcs)
colnames(pat1.filtered.bcs) <- 'mean'
pat1.filtered.bcs$sd <- rowSds(pat1.selected.bcs, na.rm = TRUE)
pat1.filtered.bcs$nonzero <- (apply(pat1.selected.bcs, 1, informative_value))*100
pat1.filtered.bcs$NAs <- ((apply(pat1.selected.bcs, 1, function(x) sum(is.na(x))))/ncol(vafs.per.cell))*100
(sum(is.na(pat1.selected.bcs))/(ncol(pat1.selected.bcs)*nrow(pat1.selected.bcs)))*100
#[1] 82.26402 percentage of NAs after low qual cell filtering, not a huge change

hist(pat1.filtered.bcs$mean, breaks = 100)
hist(pat1.filtered.bcs$sd, breaks = 100)
hist(pat1.filtered.bcs$nonzero, breaks = 100)
hist(pat1.filtered.bcs$NAs, breaks = 100)


#pat1.selected.snvs is our matrix after taking only high coverage SNVs AFTER
#filtering low qual cells
pat1.selected.snvs <- pat1.selected.bcs[pat1.filtered.bcs$NAs<cutoff$NAs & pat1.filtered.bcs$sd>cutoff$sds,]
#pat1.filtered.snvs object is once again metadata from this matrix, so only high
#quality cells AND SNVs with less than 20% NAs
pat1.filtered.snvs <- rowMeans(pat1.selected.snvs, na.rm = TRUE)
pat1.filtered.snvs <- as.data.frame(pat1.filtered.snvs)
colnames(pat1.filtered.snvs) <- 'mean'
pat1.filtered.snvs$sd <- rowSds(pat1.selected.snvs, na.rm = TRUE)
pat1.filtered.snvs$nonzero <- (apply(pat1.selected.snvs, 1, informative_value))*100
pat1.filtered.snvs$NAs <- ((apply(pat1.selected.snvs, 1, function(x) sum(is.na(x))))/ncol(pat1.selected.snvs))*100

p1_filt_mean <- ggplot(pat1.filtered.snvs, aes(x = mean)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Mean Values",
       x = "Mean",
       y = "Frequency")
p1_filt_sd <- ggplot(pat1.filtered.snvs, aes(x = sd)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Standard Deviations",
       x = "SD",
       y = "Frequency")
p1_filt_nonzero <- ggplot(pat1.filtered.snvs, aes(x = nonzero)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Non-zero values",
       x = "Percentage of Non-zeroes",
       y = "Frequency")
p1_filt_na <- ggplot(pat1.filtered.snvs, aes(x = NAs)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of NAs",
       x = "Percentage of NAs",
       y = "Frequency")

#Now we look only at cells, not SNVs. We look at columns with high coverage,
#namely, cells that have at least 1 read for more than 60% of the SNVs,
#in other words, less than 40% of non-covered SNVs
pat1.metadata.cells <- ((apply(pat1.selected.snvs, 2, function(x) sum(is.na(x))))/nrow(pat1.selected.snvs))*100
pat1.metadata.cells <- as.data.frame(pat1.metadata.cells)
colnames(pat1.metadata.cells) <- 'NAs'
table(pat1.metadata.cells$NAs < 80)
#FALSE  TRUE 
#233  1584
pat1.vafs.filtered <- pat1.selected.snvs[, pat1.metadata.cells$NAs < cutoff$cells]

#Percentage of NAs once filtered
(sum(is.na(pat1.vafs.filtered))/(ncol(pat1.vafs.filtered)*nrow(pat1.vafs.filtered)))*100
#[1] 30.97718, much lower from 80%

#Visualisation of NAs
pat1.missing_proportion <- colSums(is.na(pat1.vafs.filtered)) / nrow(pat1.vafs.filtered)
summary(pat1.missing_proportion)
pat1.na_columns <- sum(pat1.missing_proportion > 0)
pat1.na_columns
ggplot(data = melt(is.na(nuc.vafs.filtered)), aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "black"), name = "NA") +
  theme_minimal() +
  labs(title = "Heatmap of Missing Values", x = "Columns", y = "Rows")

pat1.vafs.imputed <- impute.knn(pat1.vafs.filtered)
hist(rowSds(pat1.vafs.imputed$data), breaks = 100)
pat1.srat.meta <- srat.meta[rownames(srat.meta) %in% colnames(pat1.vafs.imputed$data),]
pat1.srat.meta <- pat1.srat.meta[colnames(pat1.vafs.filtered),]
pat1.srat.annotation <- HeatmapAnnotation(df = pat1.srat.meta)

p1.ht <- ComplexHeatmap::Heatmap(pat1.vafs.imputed$data,
                        cluster_rows = TRUE,
                        cluster_columns = TRUE,
                        use_raster = FALSE,
                        na_col = '#808080',
                        col = brewer.pal(7, "Greens"),
                        show_row_names = TRUE,
                        show_column_names = FALSE,
                        row_title = 'SNVs',
                        column_title = 'Barcodes Patient 2',
                        column_km = 2,
                        column_km_repeats = 10,
                        top_annotation = pat1.srat.annotation)

p1.ht <- draw(p1.ht)

##Heatmap based on PCA ----
pcapat1 <- prcomp(t(pat1.vafs.imputed$data))
pc1 <- pcapat1$rotation[, 1]
pc2 <- pcapat1$rotation[, 2]
pat1.pc1.contributors <- rownames(pat1.vafs.imputed$data)[order(abs(pc1), decreasing = TRUE)]
pat1.pc1.contributors <- pat1.pc1.contributors[!is.na(pat1.pc1.contributors)]
pat1.pc1.contributors

pat1.pc2.contributors <- rownames(pat1.vafs.imputed$data)[order(abs(pc2), decreasing = TRUE)]
pat1.pc2.contributors <- pat1.pc2.contributors[!is.na(pat1.pc2.contributors)]
pat1.pc2.contributors

pat1.var.explained <- pcapat1$sdev^2
pat1.var.explained.percent <- pat1.var.explained / sum(pat1.var.explained) * 100
pat1.cumulative.var <- cumsum(pat1.var.explained.percent)
pat1.num.pcs <- which(pat1.cumulative.var >= cutoff$cumulative.var)[1]

pat1.df <- data.frame(PC = 1:length(pat1.var.explained.percent),
                 Variance.Explained = pat1.var.explained.percent)

ggplot(pat1.df, aes(x = PC, y = Variance.Explained)) +
  geom_point() +
  geom_line() +
  xlab("Principal Components") +
  ylab("Var explained (%)") +
  ggtitle("Elbow Plot")

## Projection of aneuploidy data on 2D
vecp2 <- srat.meta[colnames(pat1.vafs.imputed$data),'is_aneuploid']
pat1.pca.data <- as.data.frame(pcapat1$x)
pat1.pca.data$aneuploid <- vecp2

ggplot(pat1.pca.data, aes(x = PC1, y = PC2, color = aneuploid)) +
  geom_point() +
  scale_color_manual(values = c("yes" = "blue", "no" = "red", "NA" = "grey")) +
  labs(title = "PC1 vs PC2 by aneuploidity",
       x = "PC1",
       y = "PC2",
       color = "Aneuploid")


p1.umap <- umap(t(pat1.vafs.imputed$data))
colnames(p1.umap$layout) <- c('UMAP1', 'UMAP2')
test.p1.umap <- as.data.frame(p1.umap$layout)
test.p1.umap$aneuploid <- vecp2

ggplot(test.p1.umap, aes(x = UMAP1, y = UMAP2, color = aneuploid)) +
  geom_point() +
  scale_color_manual(values = c("yes" = "blue", "no" = "red", "NA" = "grey")) +
  labs(title = "UMAP1 vs UMAP2 by aneuploidity",
       x = "UMAP1",
       y = "UMAP2",
       color = "Aneuploid")

pat1.rotation.matrix <- pcapat1$rotation

# extracting top snvs from each pc
top_snvs_pc1 <- get_top_snvs(pat1.rotation.matrix, "PC1", top_n = 12)
top_snvs_pc1
#eliminate last column, because it has info used for the plot
pat1.pca.data <- pat1.pca.data[,1:ncol(pat1.pca.data)-1]

pat1.pca.data <- t(as.matrix(pat1.pca.data))
pat1.pca.data <- pat1.pca.data[1:pat1.num.pcs,]
pat1.nuc.pcs.ht <- ComplexHeatmap::Heatmap(pat1.pca.data,
                                           cluster_rows = TRUE,
                                           cluster_columns = TRUE,
                                           use_raster = FALSE,
                                           na_col = '#808080',
                                           col = brewer.pal(7, "Greens"),
                                           show_row_names = TRUE,
                                           show_column_names = FALSE,
                                           row_title = 'SNVs',
                                           column_title = 'Barcodes',
                                           column_km = 2,
                                           column_km_repeats = 10,
                                           top_annotation = pat1.srat.annotation)
pat1.nuc.pcs.ht <- draw(pat1.nuc.pcs.ht)
column_order(pat0.nuc.pcs.ht)


## Adding metadata to srat ----
p1.healthy <- colnames(pat1.vafs.filtered)[column_order(p1.ht)[[1]]]
p1.tumor <- colnames(pat1.vafs.filtered)[column_order(p1.ht)[[2]]]

p1.healthy <- paste0('LX279_', p1.healthy)
p1.tumor <- paste0('LX279_', p1.tumor)

#In the reference the patients are inversed
test[WhichCells(srat.lr, cells = p1.healthy)] <- "p0.tumor"
test[WhichCells(srat.lr, cells = p1.tumor)] <- "p0.healthy"

srat.lr <- AddMetaData(srat.lr, col.name = "test", metadata = test)

DimPlot(srat.lr,
        reduction = "umap",
        pt.size = 0.5,
        group.by = "test",
        cols = c('#AC92EB', '#4FC1E8', '#A0D568', '#FFCE54', '#ED5564')
        ) + labs(title = 'Annotation using knn') |
  FeaturePlot(srat.lr,
              reduction = "umap",
              pt.size = 0.5,
              features = "nCount_RNA")

DimPlot(
  srat.lr,
  group.by = 'test',
  cols = c('#AC92EB', '#4FC1E8', '#A0D568', '#FFCE54', '#ED5564')
) + labs(title = 'Patient and cell condition annotation using knn')|
  DimPlot(
    srat.lr,
    group.by = 'reference',
    cols = c('#4FC1E8', '#A0D568', '#FFCE54', '#ED5564')
  )+ labs(title = 'Patient and cell condition annotation using souporcell + InferCNV')

VlnPlot(srat.lr,
        features = 'nCount_RNA',
        group.by = 'test') +
  labs(title = 'N RNA counts profile using SNV annotation')

srat.lr@meta.data <- srat.lr@meta.data %>%
  mutate(matching_cases = case_when(
    test == 'p0.healthy' & is_aneuploid == 'yes' ~ 'p0.healthy_yes',
    test == 'p0.healthy' & is_aneuploid == 'no' ~ 'p0.healthy_no',
    test == 'p0.tumor' & is_aneuploid == 'yes' ~ 'p0.tumor_yes',
    test == 'p0.tumor' & is_aneuploid == 'no' ~ 'p0.tumor_no',
    test == 'p1.healthy' & is_aneuploid == 'yes' ~ 'p1.healthy_yes',
    test == 'p1.healthy' & is_aneuploid == 'no' ~ 'p1.healthy_no',
    test == 'p1.tumor' & is_aneuploid == 'yes' ~ 'p1.tumor_yes',
    test == 'p1.tumor' & is_aneuploid == 'no' ~ 'p1.tumor_no',
    TRUE ~ 'NA'
  ))

VlnPlot(srat.lr,
        features = 'nCount_RNA',
        group.by = 'matching_cases') +
  labs(title = 'Matching cases between transcriptomic and SNV aneuploidy calling')
VlnPlot(srat.lr,
        features = 'nCount_RNA',
        group.by = 'patient') +
  labs(title = 'Number of counts per cell split by patient using Souporcell annotation')

srat.lr@meta.data <- srat.lr@meta.data %>%
  mutate(matching_cases_sop = case_when(
    patient == '0' & is_aneuploid == 'yes' ~ '0_tumor',
    patient == '0' & is_aneuploid == 'no' ~ '0_healthy',
    patient == '1' & is_aneuploid == 'yes' ~ '1_tumor',
    patient == '1' & is_aneuploid == 'no' ~ '1_healthy',
    TRUE ~ 'NA'
  ))
VlnPlot(srat.lr,
        features = 'nCount_RNA',
        group.by = 'matching_cases_sop') +
  labs(title = 'Number of counts per cell split by patient using Souporcell and InferCNV annotation')

srat.lr@meta.data <- srat.lr@meta.data %>%
  mutate(reference = case_when(
    patient == '0' & is_aneuploid == 'yes' ~ 'p0.tumor',
    patient == '0' & is_aneuploid == 'no' ~ 'p0.healthy',
    patient == '1' & is_aneuploid == 'yes' ~ 'p1.tumor',
    patient == '1' & is_aneuploid == 'no' ~ 'p1.healthy',
    TRUE ~ 'NA'
  ))

table(srat.lr$test, srat.lr$reference)
#           p0.healthy p0.tumor p1.healthy p1.tumor
#NA                 83      186         71      225
#p0.healthy          0        0         73      817
#p0.tumor            0        0        271       88
#p1.healthy        634      186          0        0
#p1.tumor          597      171          0        0

#Souporcell annotations ----
##Patient demultiplexing
sop <- './sop_input/'
patient0 <- 'p0/'
patient1 <- 'p1/'

k2_nuc <- read.table(file = paste0(sop,'nuc/','clusters.tsv'), header = TRUE)
rownames(k2_nuc) <- paste0('LX279_', k2_nuc$barcode)
k2_nuc$barcode <- NULL

goodBc <- rownames(k2_nuc)[rownames(k2_nuc) %in% colnames(srat.lr)]
srat.lr@meta.data[goodBc,'patient_ont'] <- k2_nuc[goodBc,'assignment']

DimPlot(srat.lr, group.by = 'patient_ont', cols = c('#AC92EB', '#A0D568')) + labs(title = 'Annotation using souporcell on ONT') |
  DimPlot(srat.lr, group.by = 'patient_ilmn', cols = c('#AC92EB', '#4FC1E8', '#A0D568', '#FFCE54', '#ED5564')) + labs(title = 'Annotation using souporcell on Ilmn')

##k2 ----
k2_tsv.p0 <- read.table(file = paste0(sop, patient0,'clusters_k2.tsv'), header = TRUE)
k2_tsv.p1 <- read.table(file = paste0(sop,patient1,'clusters_k2.tsv'), header = TRUE)

rownames(k2_tsv.p0) <- paste0('LX279_', k2_tsv.p0$barcode)
k2_tsv.p0$barcode <- NULL

rownames(k2_tsv.p1) <- paste0('LX279_', k2_tsv.p1$barcode)
k2_tsv.p1$barcode <- NULL

k2 <- rep("NA", nrow(srat.lr@meta.data))
names(k2) <- rownames(srat.lr@meta.data)

k2_condition0 <- rownames(k2_tsv.p0)[which(k2_tsv.p0$assignment == '0')]
k2_condition1 <- rownames(k2_tsv.p0)[which(k2_tsv.p0$assignment == '1')]
k2_condition2 <- rownames(k2_tsv.p1)[which(k2_tsv.p1$assignment == '0')]
k2_condition3 <- rownames(k2_tsv.p1)[which(k2_tsv.p1$assignment == '1')]

k2[WhichCells(srat.lr, cells = k2_condition0)] <- "p1.tumor"
k2[WhichCells(srat.lr, cells = k2_condition1)] <- "p1.healthy"
k2[WhichCells(srat.lr, cells = k2_condition2)] <- "p0.healthy"
k2[WhichCells(srat.lr, cells = k2_condition3)] <- "p0.tumor"

srat.lr <- AddMetaData(srat.lr, col.name = "k2", metadata = k2)

DimPlot(srat.lr, group.by = 'k2', cols = c('#AC92EB', '#4FC1E8', '#A0D568', '#FFCE54', '#ED5564')) + labs(title = 'Annotation based on Souporcell only') |
  DimPlot(srat.lr, group.by = 'reference', cols = c('#4FC1E8', '#A0D568', '#FFCE54', '#ED5564')) + labs(title = 'Annotation based on Souporcell + InferCNV')

DimPlot(srat.lr, group.by = 'k2', cols = c('#AC92EB', '#4FC1E8', '#A0D568', '#FFCE54', '#ED5564')) + labs(title = 'Annotation based on Souporcell only') |
  DimPlot(srat.lr, group.by = 'test', cols = c('#AC92EB', '#4FC1E8', '#A0D568', '#FFCE54', '#ED5564')) + labs(title = 'Annotation based on knn')

table(srat.lr$k2, srat.lr$reference)
#           p0.healthy p0.tumor p1.healthy p1.tumor
#NA                321      131         89      187
#p0.healthy          0        0        295      177
#p0.tumor            0        0         31      766
#p1.healthy        602      243          0        0
#p1.tumor          391      169          0        0

##k3 ----
k3_tsv.p0 <- read.table(file = paste0(sop,patient0,'clusters_k3.tsv'), header = TRUE)
k3_tsv.p1 <- read.table(file = paste0(sop,patient1,'clusters_k3.tsv'), header = TRUE)

rownames(k3_tsv.p0) <- paste0('LX279_', k3_tsv.p0$barcode)
k3_tsv.p0$barcode <- NULL

rownames(k3_tsv.p1) <- paste0('LX279_', k3_tsv.p1$barcode)
k3_tsv.p1$barcode <- NULL

k3 <- rep("NA", nrow(srat.lr@meta.data))
names(k3) <- rownames(srat.lr@meta.data)

k3_condition0 <- rownames(k3_tsv.p0)[which(k3_tsv.p0$assignment == '0')]
k3_condition1 <- rownames(k3_tsv.p0)[which(k3_tsv.p0$assignment == '1')]
k3_condition2 <- rownames(k3_tsv.p0)[which(k3_tsv.p0$assignment == '2')]
k3_condition3 <- rownames(k3_tsv.p1)[which(k3_tsv.p1$assignment == '0')]
k3_condition4 <- rownames(k3_tsv.p1)[which(k3_tsv.p1$assignment == '1')]
k3_condition5 <- rownames(k3_tsv.p1)[which(k3_tsv.p1$assignment == '2')]

k3[WhichCells(srat.lr, cells = k3_condition0)] <- "p0.0"
k3[WhichCells(srat.lr, cells = k3_condition1)] <- "p0.1"
k3[WhichCells(srat.lr, cells = k3_condition2)] <- "p0.2"
k3[WhichCells(srat.lr, cells = k3_condition3)] <- "p1.0"
k3[WhichCells(srat.lr, cells = k3_condition4)] <- "p1.1"
k3[WhichCells(srat.lr, cells = k3_condition5)] <- "p1.2"

srat.lr <- AddMetaData(srat.lr, col.name = "k3", metadata = k3)

DimPlot(srat.lr, group.by = 'k3') | DimPlot(srat.lr, group.by = 'is_aneuploid') + labs(title = 'Annotation comparison Souporcell/InferCNV')



##k4 ----
k4_tsv.p0 <- read.table(file = paste0(sop,patient0,'clusters_k4.tsv'), header = TRUE)
k4_tsv.p1 <- read.table(file = paste0(sop,patient1,'clusters_k4.tsv'), header = TRUE)

rownames(k4_tsv.p0) <- paste0('LX279_', k4_tsv.p0$barcode)
k4_tsv.p0$barcode <- NULL

rownames(k4_tsv.p1) <- paste0('LX279_', k4_tsv.p1$barcode)
k4_tsv.p1$barcode <- NULL

k4 <- rep("NA", nrow(srat.lr@meta.data))
names(k4) <- rownames(srat.lr@meta.data)

k4_condition0 <- rownames(k4_tsv.p0)[which(k4_tsv.p0$assignment == '0')]
k4_condition1 <- rownames(k4_tsv.p0)[which(k4_tsv.p0$assignment == '1')]
k4_condition2 <- rownames(k4_tsv.p0)[which(k4_tsv.p0$assignment == '2')]
k4_condition3 <- rownames(k4_tsv.p0)[which(k4_tsv.p0$assignment == '3')]
k4_condition4 <- rownames(k4_tsv.p1)[which(k4_tsv.p1$assignment == '0')]
k4_condition5 <- rownames(k4_tsv.p1)[which(k4_tsv.p1$assignment == '1')]
k4_condition6 <- rownames(k4_tsv.p1)[which(k4_tsv.p1$assignment == '2')]
k4_condition7 <- rownames(k4_tsv.p1)[which(k4_tsv.p1$assignment == '3')]

k4[WhichCells(srat.lr, cells = k4_condition0)] <- "p0.0"
k4[WhichCells(srat.lr, cells = k4_condition1)] <- "p0.1"
k4[WhichCells(srat.lr, cells = k4_condition2)] <- "p0.2"
k4[WhichCells(srat.lr, cells = k4_condition3)] <- "p0.3"
k4[WhichCells(srat.lr, cells = k4_condition4)] <- "p1.0"
k4[WhichCells(srat.lr, cells = k4_condition5)] <- "p1.1"
k4[WhichCells(srat.lr, cells = k4_condition6)] <- "p1.2"
k4[WhichCells(srat.lr, cells = k4_condition7)] <- "p1.3"

srat.lr <- AddMetaData(srat.lr, col.name = "k4", metadata = k4)

DimPlot(srat.lr, group.by = 'k4') + labs(title = 'Annotation based on Souporcell only') |
  DimPlot(srat.lr, group.by = 'reference') + labs(title = 'Annotation based on Souporcell + InferCNV')




#Overlap calculation ----
srat.lr@meta.data <- srat.lr@meta.data %>%
  mutate(reference = case_when(
    patient == '0' & is_aneuploid == 'yes' ~ 'p0.tumor',
    patient == '0' & is_aneuploid == 'no' ~ 'p0.healthy',
    patient == '1' & is_aneuploid == 'yes' ~ 'p1.tumor',
    patient == '1' & is_aneuploid == 'no' ~ 'p1.healthy',
    TRUE ~ 'NA'
  ))
table(srat.lr$k2, srat.lr$reference)
#           p0.healthy p0.tumor p1.healthy p1.tumor
#NA                127      180        327      180
#p0.healthy        319      165          0        0
#p0.tumor           64      735          0        0
#p1.healthy          0        0        365      179
#p1.tumor            2        0        561      264

table(srat.lr$k2)
#                    NA p0.healthy   p0.tumor p1.healthy   p1.tumor 
# With filtering    862        466        784        710        646
# Without filtering 915        484        799        492        778




#TRASH ----
dim(illumina)
dim(ont)

shared.bcs <- intersect(colnames(illumina), colnames(ont))
illumina <- illumina[, shared.bcs]
ont <- ont[, shared.bcs]

overlapping.snvs <- intersect((rownames(illumina)), rownames(ont))
#Visualisation of NAs
pat0.missing_proportion <- colSums(is.na(pat0.vafs.filtered)) / nrow(pat0.vafs.filtered)
summary(pat0.missing_proportion)
pat0.na_columns <- sum(pat0.missing_proportion > 0)
pat0.na_columns
ggplot(data = melt(is.na(pat0.vafs.filtered)), aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "black"), name = "NA") +
  theme_minimal() +
  labs(title = "Heatmap of Missing Values", x = "Columns", y = "Rows")

#For patient ----
sop <- read.table("clusters_k2.tsv", header=T, row.names = 1)
rownames(sop) <- gsub("-1","", rownames(sop))
#colnames(denseIll) <- gsub(".1","", colnames(denseIll))
assVec <- sop[colnames(illumina),'assignment']

df = data.frame(c("Illumina", "ONT", "Intersect"),
                c(length(rownames(illumina)),
                  length(rownames(ont)), 
                  length(overlapping.snvs)))
colnames(df) <- c("Technique", "SNVs_called")
ggplot(df, aes(x = reorder(Technique,-SNVs_called), y = SNVs_called)) +
  geom_col() + scale_y_log10()

venn.cols <- brewer.pal(3, "Pastel1")
#[1] "#FBB4AE" "#B3CDE3" "#CCEBC5"
ggvenn(list(Illumina = rownames(illumina), ONT = rownames(ont)),
       fill_color = c("#FBB4AE", "#B3CDE3"))

#PCA ----
sop <- read.table("~/analysis_snv/clusters_k2.tsv", header=T, row.names = 1)
#colnames(denseIll) <- gsub(".1","", colnames(denseIll))
assVec <- sop[colnames(illumina),'assignment']
assVec <- as.numeric(assVec)
assVec[is.na(assVec)] <- -1 

pcill <- prcomp(illumina,
                center = TRUE,
                scale. = F)

pcont <- prcomp(ont,
                center = TRUE,
                scale. = F)

pc_ill_patient <-plot(pcill$rotation[,'PC1'], pcill$rotation[,'PC2'], col=assVec+1, frame=F, pch=20)
plot(pcill$rotation[assVec==0,'PC1'], pcill$rotation[assVec==0,'PC2'])
plot(pcill$rotation[assVec==1,'PC1'], pcill$rotation[assVec==1,'PC2'])

assVecOnt <- sop[colnames(ont),'assignment']
assVecOnt <- as.numeric(assVecOnt)
assVecOnt[is.na(assVecOnt)] <- -1 

plot(pcont$rotation[,'PC1'], pcont$rotation[,'PC2'], col=assVecOnt+1,pch=20)
pc_ont_patient<-plot(pcont$rotation[,'PC2'], pcont$rotation[,'PC3'], col=assVecOnt+1, pch=20, frame=F)
plot(pcont$rotation[assVecOnt==0,'PC2'], pcont$rotation[assVecOnt==0,'PC3'])
plot(pcont$rotation[assVecOnt==1,'PC1'], pcont$rotation[assVecOnt==1,'PC2'])

#For aneuploidy ----
srat.sr <- readRDS("~/analysis_illumina_LX279/srat-aneuploidy.rds")
aneupVec <- srat.sr@meta.data[colnames(srat.sr),'is_aneuploid']
bcs.sr <- gsub("-1","", colnames(srat.sr))
aneupVec <- ifelse(aneupVec == "yes", "1", ifelse(aneupVec == "no", "0", aneupVec))
aneupVec <- as.numeric(aneupVec)
aneupVec[is.na(aneupVec)] <- -1

ggplot(as.data.frame(pcill$rotation),
       aes(pcill$rotation[, 'PC1'], pcill$rotation[, 'PC2'])) + geom_point() + scale_color_manual(c("red", "green"))
plot(pcill$rotation[,'PC1'], pcill$rotation[,'PC2'], col=aneupVec+1,pch = 20, frame=F)
plot(pcill$rotation[aneupVec==0,'PC1'], pcill$rotation[aneupVec==0,'PC2'], pch = 20)
plot(pcill$rotation[aneupVec==1,'PC1'], pcill$rotation[aneupVec==1,'PC2'], pch = 20) +legend("topright",
                                                                                             legend = c("Healthy", "Tumor"),
                                                                                             fill = assVecOnt+1,       # Color of the squares
                                                                                             border = "black",
                                                                                             bty='n') # Color of the border of the squares


aneupVecOnt <- srat.sr@meta.data[colnames(srat.gene),'is_aneuploid']
aneupVecOnt <- ifelse(aneupVecOnt == "yes", "1", ifelse(aneupVecOnt == "no", "0", aneupVecOnt))
aneupVecOnt <- as.numeric(aneupVecOnt)
aneupVecOnt[is.na(aneupVecOnt)] <- -1

plot(pcont$rotation[,'PC1'], pcont$rotation[,'PC2'], col=aneupVecOnt+1, pch=20)
plot(pcont$rotation[,'PC2'], pcont$rotation[,'PC3'], col=aneupVecOnt+1, pch=20)
plot(pcont$rotation[aneupVecOnt==0,'PC2'], pcont$rotation[aneupVecOnt==0,'PC3'])
plot(pcont$rotation[aneupVecOnt==1,'PC2'], pcont$rotation[aneupVecOnt==1,'PC3'])
