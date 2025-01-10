The scripts available are meant to:
- Perform Transcriptomic analysis (ONT/ILMN)
- Run scNanoGPS
- Reformat scNanoGPS' VCF
- Unsupervised clustering from SNV VAFS/counts
- Reformat scNanoGPS' VCF to a Souporcell-friendly format
- Integration of RNA + SNV assays with Seurat.

TLDR of each file:
- LX279_ONT: transcriptomic analysis from ONT
- LX279_ONT_transcripts: isoform analysis using FULL LENGTH transcripts.
- analysis_snvs_counts: unsupervised clustering with SNVs counts to deconvolve cells (demultiplexing and malignancy inference)
- analysis_snvs_vafs: unsupervised clustering using SNVs VAFS (demultiplexing and malignancy inference)
- changing.matrix: AFTER RUNNING TRANSFORM_TO_MMF. Reformats MMF to Souporcell's format (changes -1 to 0 and adds extra info needed in the matrix)
- figures.R: replicates figures used in Report
- filtering_snvs_souporcell: filters SNVs for Souporcell's filtered alternative
- ge+snv.pipeline: unified transcriptomic + SNV analysis
- integration.Rmd: integrates RNA + SNV assays from ONT
- metacell: attempt at running metacell, UNUSED
- run_clonetracer: attempt at running clonetracer, SUCCESSFUL but UNUSED
- run_scnano: example of how scNanoGPS was run
- seurat_snvs: SNV analysis using Seurat
- subsampling: downsample of ONT data for report
- transform_to_mmf: script to turn from SNV COUNT data to Market Matrix Format (souporcell's format)
- vcf_to_allele_counts: creates VAF matrix & ALT and REF counts matrices from scNanoGPS VCF
- visualizing_sop_output: script to look at Souporcell's clustering quality
