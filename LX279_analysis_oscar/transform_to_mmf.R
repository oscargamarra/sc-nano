#First 'cycle'/step/level. Matrix including all cells to split by patient
sop.input <- './sop_input/'
patient0 <- 'p0/'
patient1 <- 'p1/'

alt.mtx <- read.table(paste0(sop.input, 'alt.mtx'), header = TRUE)
ref.mtx <- read.table(paste0(sop.input, 'ref.mtx'), header = TRUE)

library(tidyverse)
library(Matrix)

#ALT ----
alt.mtx$ID <- NULL
alt.mtx<- alt.mtx %>% unite(col='SNV',CHROM:ALT, sep=':')
rownames(alt.mtx) <- alt.mtx$SNV
alt.mtx$SNV <- NULL
alt.mtx$QUAL <- NULL
View(alt.mtx)

barcodes <- colnames(alt.mtx)
mitochondrial.snvs <- grepl("^chrM", rownames(alt.mtx))
mito.alt.mtx <- alt.mtx[mitochondrial.snvs, ]
nuc.alt.mtx <- alt.mtx[!mitochondrial.snvs, ]

rownames(nuc.alt.mtx) <- NULL
colnames(nuc.alt.mtx) <- NULL

#REF ----
ref.mtx$ID <- NULL
ref.mtx<- ref.mtx %>% unite(col='SNV',CHROM:ALT, sep=':')
rownames(ref.mtx) <- ref.mtx$SNV
ref.mtx$SNV <- NULL
ref.mtx$QUAL <- NULL
View(ref.mtx)

mito.ref.mtx <- ref.mtx[mitochondrial.snvs, ]
nuc.ref.mtx <- ref.mtx[!mitochondrial.snvs, ]

#Saving matrices ----
#Little trick is needed here, and then change manually the -1
##Nuclear ----
nuc.ref.mtx <- as.matrix(nuc.ref.mtx)
nuc.alt.mtx <- as.matrix(nuc.alt.mtx)

index.nuc.alt <- nuc.alt.mtx != 0
index.nuc.ref <- nuc.ref.mtx != 0

nuc.ref.mtx[index.nuc.alt] <- -1
nuc.alt.mtx[index.nuc.ref] <- -1

nuc.sparse.alt.mtx <- as(as.matrix(nuc.alt.mtx), "sparseMatrix")
nuc.sparse.ref.mtx <- as(as.matrix(nuc.ref.mtx), "sparseMatrix")

writeMM(nuc.sparse.alt.mtx, file = paste0(sop.input, 'nuc/market.alt.mtx'))
writeMM(nuc.sparse.ref.mtx, file = paste0(sop.input, 'nuc/market.ref.mtx'))
writeLines(barcodes, paste0(sop.input, "nuc/barcodes.txt"))

##Mito ----
mito.ref.mtx <- as.matrix(mito.alt.mtx)
mito.alt.mtx <- as.matrix(mito.ref.mtx)

index.mito.alt <- mito.alt.mtx != 0
index.mito.ref <- mito.ref.mtx != 0

mito.ref.mtx[index.mito.alt] <- -1
mito.alt.mtx[index.mito.ref] <- -1

mito.sparse.alt.mtx <- as(as.matrix(mito.alt.mtx), "sparseMatrix")
mito.sparse.ref.mtx <- as(as.matrix(mito.ref.mtx), "sparseMatrix")

writeMM(mito.sparse.alt.mtx, file = paste0(sop.input, 'mito/market.alt.mtx'))
writeMM(mito.sparse.ref.mtx, file = paste0(sop.input, 'mito/market.ref.mtx'))
writeLines(barcodes, paste0(sop.input, "mito/barcodes.txt"))

# 34 MITO SNVS SEEM INSUFFICIENT TO CLUSTER USING SOP

#With the output from sop, we now do the patient demultiplexing
#Let's add some sanity check code to see how the demult looks like.
k2_nuc <- read.table(file = paste0(sop.input, 'nuc/clusters.tsv'), header = TRUE)
srat.lr <- readRDS(file = './srat-aneuploidy.rds')

rownames(k2_nuc) <- paste0('LX279_', k2_nuc$barcode)
k2_nuc$barcode <- NULL

ont_sop <- rep("NA", nrow(srat.lr@meta.data))
names(ont_sop) <- rownames(srat.lr@meta.data)

patient0 <- rownames(k2_nuc)[k2_nuc$assignment == '0']
patient1 <- rownames(k2_nuc)[k2_nuc$assignment == '1']

ont_sop[WhichCells(srat.lr, cells = patient0)] <- "0"
ont_sop[WhichCells(srat.lr, cells = patient1)] <- "1"

srat.lr <- AddMetaData(srat.lr, col.name = "ont_sop", metadata = ont_sop)

DimPlot(srat.lr, group.by = 'ont_sop')

table(k2_nuc$assignment)
#0  0/1    1  1/0 
#2208   93 2021  101 

#Split by patient ----
#ALT ----
##Patient 0 ----
p0.alt.mtx <- as.matrix(nuc.alt.mtx[, colnames(alt.mtx) %in% gsub('LX279_','',patient0)])
p0.barcodes <- colnames(p0.alt.mtx)
colnames(p0.alt.mtx) <- NULL
rownames(p0.alt.mtx) <- NULL

##Patient 1 ----
p1.alt.mtx <- as.matrix(nuc.alt.mtx[, colnames(alt.mtx) %in% gsub('LX279_','',patient1)])
p1.barcodes <- colnames(p1.alt.mtx)

p1.snvs <- rownames(p1.alt.mtx)
colnames(p1.alt.mtx) <- NULL
rownames(p1.alt.mtx) <- NULL

#REF ----
##Patient 0 ----
p0.ref.mtx <- as.matrix(nuc.ref.mtx[, colnames(ref.mtx) %in% gsub('LX279_','',patient0)])
p0.barcodes <- colnames(p0.ref.mtx)
colnames(p0.ref.mtx) <- NULL
rownames(p0.ref.mtx) <- NULL

##Patient 1 ----
p1.ref.mtx <- as.matrix(nuc.ref.mtx[, colnames(ref.mtx) %in% gsub('LX279_','',patient1)])
p1.barcodes <- colnames(p1.ref.mtx)
colnames(p1.ref.mtx) <- NULL
rownames(p1.ref.mtx) <- NULL

#Saving matrices ----
#Little trick is needed here, and then change manually the -1
##Patient 0 ----
p0.alt.mtx <- as.matrix(p0.alt.mtx)
p0.ref.mtx <- as.matrix(p0.ref.mtx)

index.p0.alt <- p0.alt.mtx != 0
index.p0.ref <- p0.ref.mtx != 0

p0.ref.mtx[index.p0.alt] <- -1
p0.alt.mtx[index.p0.ref] <- -1

p0.sparse.alt.mtx <- as(p0.alt.mtx, "sparseMatrix")
p0.sparse.ref.mtx <- as(p0.ref.mtx, "sparseMatrix")

writeMM(p0.sparse.alt.mtx, file = paste0(sop.input, 'p0/market.alt.mtx'))
writeMM(p0.sparse.ref.mtx, file = paste0(sop.input, 'p0/market.ref.mtx'))
writeLines(p0.barcodes, paste0(sop.input, "p0/barcodes.txt"))

##Patient 1 ----
p1.alt.mtx <- as.matrix(p1.alt.mtx)
p1.ref.mtx <- as.matrix(p1.ref.mtx)

index.p1.alt <- p1.alt.mtx != 0
index.p1.ref <- p1.ref.mtx != 0

p1.ref.mtx[index.p1.alt] <- -1
p1.alt.mtx[index.p1.ref] <- -1

p1.sparse.alt.mtx <- as(p1.alt.mtx, "sparseMatrix")
p1.sparse.ref.mtx <- as(p1.ref.mtx, "sparseMatrix")

writeMM(p1.sparse.alt.mtx, file = paste0(sop.input, 'p1/market.alt.mtx'))
writeMM(p1.sparse.ref.mtx, file = paste0(sop.input, 'p1/market.ref.mtx'))
writeLines(p1.barcodes, paste0(sop.input, "p1/barcodes.txt"))
