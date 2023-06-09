# Rmd

R markdown files for analysis following read processing.

## QC_PCA.Rmd

Quality Control by Scatter plots and PCA.

Also organises the read counts and outputs organised summaries, `counts_bysamplecode.txt` and `tpms_bysamplecode.txt`, in `results/summaries` folder.

## TPM_plots.Rmd

Plot RNA abundances in transcripts per million.
Uses data in `tpms_bysamplecode.txt`, and the sample sheet.

Under construction as we haven't chosen which genes to focus on.

## DGE...Rmd - differential gene expression analyses.

Focus on pairwise differential gene expression analyses at the 1 hr timepoint.

Comparing temperature 37C vs 25C in each media, also comparing media RPMI+ to YPD at each temperature.

Results go in `results/DEGs` - deseqdf with all results, and deglist with just lists of up- and down-regulated genes.

### `DGE_37Cvs25C_RPMI1hr.Rmd`

### `DGE_37Cvs25C_YPD1hr.Rmd`

### `DGE_RPMIvsYPD_25C1hr.Rmd`

### `DGE_RPMIvsYPD_37C1hr.Rmd`

### `Combine_DGE_tables.Rmd`

Combines tables into a single .xlsx file for output: `/results/geneclusters/CryptoWakeup_geneclusters.xlsx`
## MORE TO ADD
