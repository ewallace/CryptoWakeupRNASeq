# CryptoWakeupRNASeq/results

Results of data analysis in various formats

# counts

RNA-seq data analysis directly output from the nextflow pipeline `process_tagseq.nf`, including:

- Each sample has a folder containing FastQC quality control, hisat2 alignment summary, bedgraphs for genome browser plots on both plus and minus strands. These folders have names described in the sample sheet `/input/experiment/CryptoWakeUpSampleSheetPlusTagsAndRawFilenames.txt`
- `multiqc_report.html` MultiQC report summarising statistics across all samples
- `nextflow_report.html` nextflow report summarising the programs that were run by `process_tagseq.nf`
- `flowchart.png` is a picture describing the programs run by `process_tagseq.nf`
- `counts.txt` counts per gene and per sample as collected by featureCounts. A cleaned version of this data, easier for humans to read, is found in `/results/summaries/counts_bysamplecode.txt`
- `counts.txt.summary` is the summary of FeatureCounts output, assigned counts per sample.


# counts_test100000

As for folder `counts`, but for a test run of 100000 reads.

Ignore this, unless you are troubleshooting re-running the analysis pipeline.

# summaries

Cleaned and nicely labeled summaries per gene and per sample, created by script `/Rmd/QC_PCA.Rmd` using the sample codes (labels) found in `/results/summaries/counts_bysamplecode.txt`

- `counts_bysamplecode.txt` - raw counts per gene
- `tpms_bysamplecode.txt` - transcripts per million (tpm) per gene


# geneclusters

Clusters of co-regulated genes, output by `Rmd/Heatmap.Rmd` - see the script for full details. This contains just the genelists

- `genecluster_GAT201_k700.txt` - cluster containing GAT201, induced in RPMI
- `genecluster_heat_k300.txt` - cluster induced by heat
- `genecluster_stationary_k300.txt` - cluster up in stationary phase
- `genecluster_translation_k300.txt` - cluster of translation-related genes
- `genecluster_TSA3_k500.txt` - cluster containing TSA3, strongly induced in RPMI



# geneclusters_GOKEGG

Functional analysis by gene ontology and KEGG of the gene clusters above.
This analysis was performed on FungiDB in March 2023.

# DEGs

Differentially expressed gene analysis. Both data frames of statistics for all genes `deseq_df....txt` and genelists of differentially up or down-regulated genes 2-fold at 5% false discovery rate `deglist_...2x_FDR5.txt`.

These lists are made by scripts in the `Rmd` folder, see the scripts for details:

- `DGE_RPMIvsYPD_37C1hr.Rmd`
- `DGE_37Cvs25C_RPMI1hr.Rmd`
- `DGE_37Cvs25C_YPD1hr.Rmd`
- `DGE_RPMI37C2h_vs_0min.Rmd`
- `DGE_RPMIvsYPD_25C1hr.Rmd`
