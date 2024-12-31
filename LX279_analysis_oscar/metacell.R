library(metacell)
library(tidyverse)
library(tidyr)
library(Seurat)

if(!dir.exists("testdb")) dir.create("testdb/")
scdb_init("testdb/", force_reinit=T)

alt <- read.table(file = '~/sop_input/alt.mtx', header = TRUE)
alt$ID <- NULL
alt<- alt %>% unite(col='SNV',CHROM:ALT, sep=':')
rownames(alt) <- alt$SNV
rownames(alt) <- paste0('alt:', rownames(alt))
alt$SNV <- NULL
alt$QUAL <- NULL
alt <- as.matrix(alt)

ref <- read.table(file = '~/sop_input/ref.mtx', header = TRUE)
ref$ID <- NULL
ref<- ref %>% unite(col='SNV',CHROM:ALT, sep=':')
rownames(ref) <- ref$SNV
rownames(ref) <- paste0('ref:', rownames(ref))
ref$SNV <- NULL
ref$QUAL <- NULL
ref <- as.matrix(ref)

reads <- (rbind(alt, ref))

write.table(reads, file="./testdb/reads.tsv", sep="\t", row.names=TRUE, col.names=TRUE)
srat.lr <- readRDS(file = "./nanopore_souporcell/srat-aneuploidy.rds")

mcell_import_scmat_tsv(mat_nm = 'rbind', fn = './testdb/reads.tsv', dset_nm = 'rbind')
mat = scdb_mat("rbind")
print(dim(mat@mat))

if(!dir.exists("./testdb/figs")) dir.create("./testdb/figs/")
scfigs_init("./testdb/figs/")

#plot is in /testdb/figs
nms <- c(rownames(mat@mat), rownames(mat@ignore_gmat))
bad_genes = unique(c(grep("^alt:chrM", nms, v=T), grep("^ref:chrM", nms, v=T)))
mcell_mat_ignore_genes(new_mat_id="rbind", mat_id="rbind", bad_genes, reverse=F)

mcell_plot_umis_per_cell("rbind", min_umis_cutoff = NA, bin_for_cutoff = 40)

#Filtering based on previous plot
not_in_srat <- setdiff(colnames(mat@mat),(gsub('LX279_', '', colnames(srat.lr))))
mcell_mat_ignore_cells("rbind_filt","rbind",ig_cells=c(not_in_srat))
mat <- scdb_mat("rbind_filt")
print(dim(mat@mat))


#This creates a gstat object named rbind. To retrieve, use gstats=scdb_gstat("rbind")
mcell_add_gene_stat(gstat_id="rbind", mat_id="rbind_filt", force=T)


#This func creates a similarity graph
mcell_add_cgraph_from_mat_bknn(mat_id="rbind_filt", 
                               gset_id = "rbind_feats", 
                               graph_id="rbind_graph",
                               K=100, #number of neighbours per metacell
                               dsamp=F)


mcell_coclust_from_graph_resamp(
  coc_id="rbind_coc500", 
  graph_id="rbind_graph",
  min_mc_size=20, 
  p_resamp=0.75, n_resamp=500)

mcell_mc_from_coclust_balanced(
  coc_id="rbind_coc500", 
  mat_id= "rbind",
  mc_id= "rbind_mc", 
  K=30, min_mc_size=30, alpha=2)

mcell_plot_outlier_heatmap(mc_id="rbind_mc", mat_id = "rbind", T_lfc=3)

mcell_mc_split_filt(new_mc_id="rbind_mc_f", 
                    mc_id="rbind_mc", 
                    mat_id="rbind",
                    T_lfc=3, plot_mats=F)

mcell_gset_from_mc_markers(gset_id="rbind_markers", mc_id="rbind_mc_f")

mcell_mc_plot_marks(mc_id="rbind_mc_f", gset_id="rbind_markers", mat_id="rbind")

mcell_mc2d_force_knn(mc2d_id="rbind_2dproj",mc_id="rbind_mc_f", graph_id="rbind_graph")

tgconfig::set_param("mcell_mc2d_height",1000, "metacell")

tgconfig::set_param("mcell_mc2d_width",1000, "metacell")

mcell_mc2d_plot(mc2d_id="rbind_2dproj")