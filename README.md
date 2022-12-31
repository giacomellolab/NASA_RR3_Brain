# NASA_RR3_Brain
Here, one can find the scripts used in the analysis of the mouse brain data used in the study titled "**_Neurological effects of spaceflight on mouse brains revealed by spatial and single-cell multiomics_**"

### Spatial Transcriptomics data analysis
The code for spatial transcriptomics data analysis can be found under the [code](spatial-transcriptomics/code/RR3_brain_ST_clustering.Rmd) section of the spatial-transcriptomics folder. Corresponding metadata and geneslist used for filtering is also provided in the [data](spatial-transcriptomics/data) section.
Processed count matrices from spaceranger used in the analysis, original HnE tissue images, and the final seurat object are available under this dataset in Mendeley Data.

### Multiomics data analysis
Combined code for snRNAseq, ATACseq, integrated analysis, differential expression analysis, motif analysis is available in rmarkdown format under the [code](multiomics/code/real_data_B_multiomics_20221222.rmd) section of multiomics folder. The corresponding html document with results and plots, and final seurat object with corresponding metadata can be downloaded from Mendeley Dataset [here].

### Deconvolution
ST spot decomposition was performed using Stereoscope. The steps include:
1. Prepare data: First, create h5ad objects of both ST and SN data. For SN data, subsample by cluster to similar number of cells per cluster and run DEG detection. Select a geneset for deconvolution. Was done with either top20, top50 or top100 genes per cluster, or by taking top 2000 variable genes. Also, only select genes that are expressed in the ST data as well.
2. Run deconvolution: [Bash](deconvolution-stereoscope/code/run_stereoscope_brain.bash) script (have only changed the line with SC_LIST to define different gene sets)
3. Summarize results: Finally, summarize the results using this [script](deconvolution-stereoscope/code/summary_stereoscope_brain.Rmd).

The deconvolution results saved as an R object can be downloaded from the Mendeley Dataset [here].

### Ligand Receptor interactions analysis
L-R interactions were annalysed via this [code](L-R_interactions/code/L_R_Brain.Rmd). Resultant plots are shown and discussed in the manuscript.



