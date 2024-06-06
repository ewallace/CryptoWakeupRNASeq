# Rmd

R markdown files for analysis following read processing.

## QC_PCA.Rmd

Quality Control by Scatter plots and PCA.

Also organises the read counts and outputs organised summaries, `counts_bysamplecode.txt` and `tpms_bysamplecode.txt`, in `results/summaries` folder.

## `load_count_data.Rmd`

To load and organise the raw count data in R, as well as loading libraries for data analysis and plotting.

Other scripts just run this one instead of duplicating code, especially scripts for differential gene expression analysis.

## `Heatmap.Rmd`

Makes clustered heatmap summarising the dataset.


## TPM_plots.Rmd

Plot RNA abundances in transcripts per million.
Uses data in `tpms_bysamplecode.txt`, and the sample sheet.


## DGE...Rmd - differential gene expression analyses.

Focus on pairwise differential gene expression analyses, with DESeq2, mostly at the 1 hr timepoint.

Comparing temperature 37C vs 25C in each media, also comparing media RPMI+ to YPD at each temperature.
Another file compares RPMI+ 37C 2hr to 0 min.

Results go in `results/DEGs` - deseqdf with all results, and deglist with just lists of up- and down-regulated genes.


### `DGE_37Cvs25C_RPMI1hr.Rmd`

### `DGE_37Cvs25C_YPD1hr.Rmd`

### `DGE_RPMIvsYPD_25C1hr.Rmd`

### `DGE_RPMIvsYPD_37C1hr.Rmd`

### `DGE_RPMI37C2h_vs_0min.Rmd`

RPMI 37C 2hr samples compared to pre-wakeup 0 min (= 5 days in YPD).

### `Combine_DGE_tables.Rmd`

Combines DGE tables into a single .xlsx file for output: `/results/geneclusters/CryptoWakeup_geneclusters.xlsx`



