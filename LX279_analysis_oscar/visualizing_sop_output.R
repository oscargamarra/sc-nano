k2_tsv.p0 <- read.table(file = '~/sop_input/input_filtered/p0/clusters.tsv', header = T)
k2_tsv.p1 <- read.table(file = '~/sop_input/input_filtered/p1/clusters.tsv', header = T)

cluster_values.p0 <- data.frame(cluster0 = k2_tsv.p0$cluster0, cluster1 = k2_tsv.p0$cluster1)
cluster_values.p1 <- data.frame(cluster0 = k2_tsv.p1$cluster0, cluster1 = k2_tsv.p1$cluster1)

cluster_values.p0 <- BBmisc::normalize(
  cluster_values.p0,
  method = "standardize",
  range = c(0, 1),
  margin = 1L,
  on.constant = "quiet"
)
cluster_values.p1 <- BBmisc::normalize(
  cluster_values.p1,
  method = "standardize",
  range = c(0, 1),
  margin = 1L,
  on.constant = "quiet"
)

cluster_values.p0$assignment <- k2_tsv.p0$assignment
cluster_values.p0$status <- k2_tsv.p0$status

cluster_values.p1$assignment <- k2_tsv.p1$assignment
cluster_values.p1$status <- k2_tsv.p1$status

ggplot(data = cluster_values.p0, aes(x = cluster0, y=cluster1)) + geom_point(aes(colour = assignment))
ggplot(data = cluster_values.p1, aes(x = cluster0, y=cluster1)) + geom_point(aes(colour = assignment))

srat.lr <- readRDS(file = './nanopore_souporcell/srat-aneuploidy.rds')

##k2 ----
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

k2[WhichCells(srat.lr, cells = k2_condition0)] <- "p0.healthy"
k2[WhichCells(srat.lr, cells = k2_condition1)] <- "p0.tumor"
k2[WhichCells(srat.lr, cells = k2_condition2)] <- "p1.healthy"
k2[WhichCells(srat.lr, cells = k2_condition3)] <- "p1.tumor"


srat.lr <- AddMetaData(srat.lr, col.name = "k2", metadata = k2)

DimPlot(srat.lr, group.by = 'k2', cols = c('#AC92EB', '#4FC1E8', '#A0D568', '#FFCE54', '#ED5564')) + labs(title = 'Souporcell only - Filtered by SD') |
  DimPlot(srat.lr, group.by = 'reference', cols = c('#4FC1E8', '#A0D568', '#FFCE54', '#ED5564')) + labs(title = 'Reference (souporcell + InferCNV Ilmn)')

table(srat.lr$k2, srat.lr$reference)


srat.sr <- readRDS(file = './srat-cellsidentified.rds')

VlnPlot(srat.sr, features = 'percent_mito', group.by = 'patient')
VlnPlot(srat.lr, features = 'percent_mito', group.by = 'patient')