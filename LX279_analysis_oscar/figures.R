srat <- readRDS(file = '~/analysis_illumina_LX279/srat-aneuploidy.rds')
srat.lr <- readRDS(file = '~/analysis_nanopore_LX279/genes/srat-aneuploidy.rds')
srat.transcript <- readRDS(file = '~/analysis_nanopore_LX279/transcripts/srat-enriched.rds')
srat.snvs <- readRDS(file = '~/srat.snvs-clusters.rds')

#We will also include the counts coming from the downsampled Illumina
counts.sr.ds <- Read10X(data.dir = "~/DS_LX279_cellranger_outs/filtered_feature_bc_matrix")
counts.sr <- Read10X(data.dir = "filtered_matrix_LX279_short_reads")
counts <- read.delim("~/analysis_nanopore_LX279/genes/LX279_sup_complete.gene_expression.counts.tsv",
                     row.names=1)

#Change colnames, as they have different format in Ilmn
colnames(counts.sr) <- gsub('-1', '', colnames(counts.sr))
colnames(counts.sr.ds) <- gsub('-1', '', colnames(counts.sr.ds))

#Subset count matrices to bcs only found in both modalities
intersect.ds.ilmn.matrix <- counts.sr.ds[, colnames(counts.sr.ds) %in% gsub('LX279_', '', colnames(srat.lr))]
intersect.ilmn.matrix <- counts.sr[, colnames(counts.sr) %in% colnames(counts.sr.ds)]
intersect.ont.matrix <- counts[, colnames(counts) %in% colnames(counts.sr.ds)]

#To subset srats also to bcs only present in both srats after analysis
intersect.bcs <- intersect(colnames(srat), colnames(srat.lr))

#Use paste0 because it is easier for downsampling later the srat object
intersect.subset.bcs <- intersect(paste0('LX279_', colnames(counts.sr.ds)), colnames(srat.lr))

##And subset based on them
intersect.srat.ilmn <- subset(x = srat, cells = intersect.bcs)
intersect.srat.ONT <- subset(x = srat.lr, cells = intersect.bcs)

intersect.subset.srat.ONT <- subset(x = srat.lr, cells = intersect.subset.bcs)

#Filtered cells
df_filtered <- data.frame(
  Ilmn_counts = Matrix::colSums(intersect.srat.ilmn@assays$RNA$counts),
  ONT_counts = Matrix::colSums(intersect.srat.ONT@assays$RNA$counts),
  Ilmn_features = colSums(intersect.srat.ilmn@assays$RNA$counts > 0),
  ONT_features = colSums(intersect.srat.ONT@assays$RNA$counts > 0)
) 

#Subset Ilmn
df_subset <- data.frame(
  DS_Ilmn_counts = Matrix::colSums(intersect.ds.ilmn.matrix),
  ONT_counts = Matrix::colSums(intersect.subset.srat.ONT@assays$RNA$counts),
  DS_Ilmn_features = colSums(intersect.ds.ilmn.matrix > 0),
  ONT_features = colSums(intersect.subset.srat.ONT@assays$RNA$counts > 0)
)

rownames(df_filtered) <- gsub('LX279_', '', rownames(df_filtered))

df_filtered$singler <- srat.lr$singler[rownames(df_filtered)]
df_filtered$is_aneuploid <- srat.lr$is_aneuploid[rownames(df_filtered)]
df_filtered$patient_ont <- srat.lr$patient_ont[rownames(df_filtered)]
df_filtered$normalized_Ilmn <- df_filtered$Ilmn_counts/df_filtered$Ilmn_features
df_filtered$normalized_ONT <- df_filtered$ONT_counts/df_filtered$ONT_features

#Now for subset too, we ignore metadata coming from Seurat pipeline
df_subset$normalized_DS_Ilmn <- df_subset$DS_Ilmn_counts/df_subset$DS_Ilmn_features
df_subset$normalized_ONT <- df_subset$ONT_counts/df_subset$ONT_features


#Raw reads
rownames(df_unfiltered) <- gsub('LX279_', '', rownames(df_unfiltered))

df_unfiltered$singler <- srat.lr$singler[rownames(df_unfiltered)]
df_unfiltered$is_aneuploid <- srat.lr$is_aneuploid[rownames(df_unfiltered)]
df_unfiltered$patient_ont <- srat.lr$patient_ont[rownames(df_unfiltered)]
df_unfiltered$normalized_Ilmn <- df_unfiltered$Ilmn_counts/df_unfiltered$Ilmn_features
df_unfiltered$normalized_ONT <- df_unfiltered$ONT_counts/df_unfiltered$ONT_features

colors <- RColorBrewer::brewer.pal(n = 12, name = "Paired")

model <- lm(ONT_features ~ Ilmn_features, data = df_unfiltered)

#For R squared value
r_squared <- summary(model)$r.squared


ggplot(df_unfiltered, aes(x = normalized_Ilmn, y = normalized_ONT, color = singler)) +
  geom_point(size = 1) +
  flatuicoloRs::scale_color_spanish() +
  theme_minimal() +
  labs(title = 'Correlation between transcript reads using Ilmn and ONT grouped by Cell Type', x = "Ilmn Features", y = "ONT Features", color = "Cell Type") +
  geom_smooth(
  method = lm ,
  color = "black",
  fill = "#69b3a2") +
  theme_minimal() +
  annotate("text", x = 4, y = 4, label = paste("R² = ", round(r_squared, 3)), size = 5)

#Scatter of ratio UMI/Genes per platform
ggscatter(df_filtered, x = "normalized_Ilmn", y = "normalized_ONT",size = 0.5, color = 'black',
          add = "reg.line", conf.int = FALSE, 
          cor.coef = TRUE, cor.method = "spearman",
          add.params = list(color = "blue")) +
  labs(title = 'Comparison of UMI to Gene ratios between Illumina and ONT Platforms', x = "UMI/Gene Ratio (Illumina)", y = "UMI/Gene Ratio (ONT)") +
  theme_minimal()

#Scatter of ratio UMI/Genes per platform, now downsampling Ilmn
ggscatter(df_subset, x = "normalized_DS_Ilmn", y = "normalized_ONT",size = 0.5, color = 'black',
          add = "reg.line", conf.int = FALSE, 
          cor.coef = TRUE, cor.method = "spearman",
          add.params = list(color = "blue")) +
  labs(title = 'Comparison of UMI to Gene ratios between Illumina and ONT Platforms', x = "UMI/Gene Ratio (Illumina)", y = "UMI/Gene Ratio (ONT)") +
  theme_minimal()

#Comparison of cell annotation Ilmn - ONT
##First rename colnames from srats so they are equal
colnames(srat) <- gsub('-1', '', colnames(srat))

##Get matching barcodes
intersect.bcs <- intersect(colnames(srat), colnames(srat.lr))
##And subset based on them
intersect.srat.ilmn <- subset(x = srat, cells = intersect.bcs)
intersect.srat.ONT <- subset(x = srat.lr, cells = intersect.bcs)

##Now we see matching cases for cell annotation
table(intersect.srat.ONT$singler, intersect.srat.ilmn$singler)

## Change annotation of underrepresented cell types ----
# Total amount of cells
cell_percentage_ilmn <- sum(table(intersect.srat.ilmn$singler))

# Renaming
intersect.srat.ilmn$singler <- ifelse(
  table(intersect.srat.ilmn$singler)[intersect.srat.ilmn$singler] / cell_percentage_ilmn < 0.01,
  'Other',
  intersect.srat.ilmn$singler
)

table(intersect.srat.ilmn$singler)

#Same for ONT
cell_percentage_ont <- sum(table(intersect.srat.ONT$singler))

# Renaming
intersect.srat.ONT$singler <- ifelse(
  table(intersect.srat.ONT$singler)[intersect.srat.ONT$singler] / cell_percentage_ont < 0.01,
  'Other',
  intersect.srat.ONT$singler
)

table(intersect.srat.ONT$singler)

#Less confusing table
table(intersect.srat.ONT$singler, intersect.srat.ilmn$singler)

#UMAPs ----
##Cell typing Illumina
DimPlot(intersect.srat.ilmn, group.by = 'singler', pt.size = 0.5) +
  ggsci::scale_color_simpsons() +
  labs(title = 'Cell annotation in Ilmn') +
  xlab('UMAP_1') + ylab('UMAP_2')

#Genotype deconvolution and aneuploidy inference Ilmn
DimPlot(intersect.srat.ilmn, group.by = 'reference', pt.size = 0.5) +
  ggsci::scale_color_jco() +
  labs(title = 'Genotype and CNV inference') +
  xlab('UMAP_1') + ylab('UMAP_2')

##Cell typing in ONT
DimPlot(intersect.srat.ONT, group.by = 'singler', pt.size = 0.5) +
  ggsci::scale_color_simpsons() +
  labs(title = 'Cell annotation in ONT') +
  xlab('UMAP_1') + ylab('UMAP_2')

#Genotype deconvolution and aneuploidy inference ONT
DimPlot(intersect.srat.ONT, group.by = 'reference', pt.size = 0.5) +
  ggsci::scale_color_jco() +
  labs(title = 'Genotype and CNV inference ONT') +
  xlab('UMAP_1') + ylab('UMAP_2')

DimPlot(intersect.srat.ONT, group.by = 'singler', pt.size = 0.5) +
  ggsci::scale_color_simpsons()

#SNV based UMAP
#Simplify annotation
cell_percentage_snv <- sum(table(srat.snvs$singler))

srat.snvs$singler <- ifelse(
  table(srat.snvs$singler)[srat.snvs$singler] / cell_percentage_snv < 0.01,
  'Other',
  srat.snvs$singler
) + 
  
  table(srat.snvs$singler)

DimPlot(srat.snvs, group.by = 'singler', pt.size = 0.5) +
  ggsci::scale_color_simpsons() +
  labs(title = 'Cell annotation transfered to SNV-based UMAP') +
  xlab('UMAP_1') + ylab('UMAP_2')

DimPlot(srat.snvs, group.by = 'reference', pt.size = 0.5) +
  ggsci::scale_color_jco() +
  labs(title = 'Genotype and aneuploidy inference transferred to SNV-based UMAP') +
  xlab('UMAP_1') + ylab('UMAP_2')

#Log2 fold change of Ilmn against ONT cell typing ----
matching_cell_typing <- as.vector(log2(table(intersect.srat.ilmn$singler)) - log2(table(intersect.srat.ONT$singler)))

df_celltype <- data.frame(FC = matching_cell_typing,
                  cell_type = names(table(intersect.srat.ONT$singler)))

#Ilmn against ONT, >0 values means overrepresentation in Ilmn in comparison to ONT
ggplot(df_celltype, aes(x = FC, y = cell_type, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() + labs(title = 'Log 2 Fold change of Cells identified in Ilmn against ONT') +
  xlab("log2FC Ilmn vs ONT") + ylab("Cell type") + scale_y_discrete(limits =
                                                                      rev) +
  ggsci::scale_fill_simpsons() 

#UMAP transcript with Singler from ONT
srat.transcript <- AddMetaData(srat.transcript, colnames(srat.transcript), col.name = 'bcs')
srat.lr <-
  AddMetaData(srat.lr, colnames(srat.lr), col.name = 'bcs')
merged.metadata <- merge(srat.transcript@meta.data,
                         srat.lr@meta.data[, c('singler', 'bcs')], by =
                           'bcs')
srat.transcript <- srat.transcript[, colnames(srat.transcript) %in% merged.metadata$bcs]
srat.transcript@meta.data$singler <- merged.metadata$singler.x

#Simplify cell annotation
cell_percentage_trans <- sum(table(srat.transcript$singler))

srat.transcript$singler <- ifelse(
  table(srat.transcript$singler)[srat.transcript$singler] / cell_percentage_trans < 0.01,
  'Other',
  srat.transcript$singler
)

table(srat.transcript$singler)

DimPlot(srat.transcript, group.by = 'singler', pt.size = 0.5) +
  ggsci::scale_color_simpsons() +
  labs(title = 'Cell annotation from ONT transcript') +
  xlab('UMAP_1') + ylab('UMAP_2')

#Log2 fold change of Ilmn against ONT cell typing ----
intersect.srat.ONT$reference <- gsub("p0", "temp", intersect.srat.ONT$reference)  # Replace 'p0' with a temporary placeholder
intersect.srat.ONT$reference <- gsub("p1", "p0", intersect.srat.ONT$reference)    # Replace 'p1' with 'p0'
intersect.srat.ONT$reference <- gsub("temp", "p1", intersect.srat.ONT$reference)  # Replace temporary placeholder with 'p1'

matching_aneuploidy <- as.vector(log2(table(intersect.srat.ilmn$reference)) - log2(table(intersect.srat.ONT$reference)))

df_aneuploidy <- data.frame(FC = matching_aneuploidy,
                          cell_type = names(table(intersect.srat.ONT$reference)))

#Ilmn against ONT for aneuploidy, >0 values means overrepresentation in Ilmn in comparison to ONT
ggplot(df_aneuploidy, aes(x = FC, y = cell_type, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() + labs(title = 'Log 2 Fold change of Aneuploidy detection in Ilmn against ONT') +
  xlab("log2FC Ilmn vs ONT") + ylab("Patient and Cell State") + scale_y_discrete(limits =
                                                                      rev) +
  ggsci::scale_fill_jco() 


#Read length ONT ----
#Ilmn has 90bp
read_length <- read.table("~/analysis_nanopore_LX279/hist_readlength", header = FALSE, sep = "", col.names = c("length", "count"))

# rearrange df (1st col is read length, 2nd col is n of counts with such read length)
expanded_lengths <- rep(read_length$length, read_length$count)
expanded_lengths<- log10(expanded_lengths)

# subsetting because there are too many values
set.seed(123)  # Para reproducibilidad
subset_lengths <- sample(expanded_lengths, 4000)

# density + rug plot
ggplot(data.frame(length = expanded_lengths), aes(x = length)) +
  geom_density(fill = "blue", alpha = 0.5) +  # Plot de densidad con todos los valores
  geom_rug(data = data.frame(length = subset_lengths), aes(x = length), sides = "b", color = "black", alpha = 0.5) +  # Rug plot solo con el subset
  labs(title = "Density Plot of Read Lengths in ONT", x = "Log10 Base Pair Length", y = "Density") +
  theme_minimal()

#Heatmap cell typing ----
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(reshape2)

# Create a confusion matrix to calculate overlaps
confusion_matrix <- table(intersect.srat.ilmn$singler, intersect.srat.ONT$singler)

# Calculate percentage overlaps by dividing each cell by row sums
percentage_overlap <- prop.table(confusion_matrix, 1) * 100

# Convert to a data frame for ggplot
percentage_overlap_df <- as.data.frame(as.table(percentage_overlap))

# Plot heatmap
ggplot(percentage_overlap_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "indianred") +
  labs(
    title = ' Overlap of cell type identification between ONT and Illumina',
    x = "Cell Type in ONT",
    y = "Cell Type in Illumina",
    fill = "Overlap (%)"
  ) +
  geom_text(aes(label = sprintf("%.1f", Freq)), color = "black", size = 4) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(limits=rev)

#AUROC genotypes ----
#Loading genotypes from souporcell

k2_standard <- read.table(file = '~/sop_input/clusters.tsv', header = T)
cluster_values <- data.frame(cluster0 = k2_standard$cluster0, cluster1 = k2_standard$cluster1)
cluster_values$assignment <- k2_standard$assignment


