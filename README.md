# NASA_RR3_Brain
Here, one can find the scripts used in the analysis of the mouse brain data used in the study titled "**_Spatially resolved multiomics on the neuronal effects induced by spaceflight_**"

*Impairment of the central nervous system (CNS) functions in astronauts is a major health risk for long-duration space missions. Here, for the first time, we combine single-cell multiomics (transcriptomics and chromatin accessibility) and spatial transcriptomics analyses to discover spaceflight-mediated changes in the mouse brain. By comparing ground control and spaceflight animals, we found that the main processes affected by spaceflight include neurogenesis, synaptogenesis and synaptic transmission in cortex, hippocampus, striatum and neuroendocrine structures as well as astrocyte activation and immune dysfunction. At the pathway level, spaceflight resembles neurodegenerative diseases with oxidative stress and protein misfolding components, such as in Parkinson’s disease. Our integrated spatial multiomics approach reveals both widespread and localized brain impairments and suggests key structures and mechanisms to be targeted for countermeasure development. All datasets can be mined through an interactive data portal as well as the National Aeronautics and Space Administration (NASA) GeneLab repositories.*

## Analysis and scripts

### Spatial Transcriptomics data analysis
The code for spatial transcriptomics data analysis can be found under the [code](spatial-transcriptomics/code/RR3_brain_ST_clustering.Rmd) section of the spatial-transcriptomics folder. Corresponding metadata and geneslist used for filtering is also provided in the [data](spatial-transcriptomics/data) section.
Processed count matrices from spaceranger used in the analysis, original HnE tissue images, and the final seurat object are available under this dataset in Mendeley Data.

### Multiomics data analysis
Combined code for snRNAseq, ATACseq, integrated analysis, differential expression analysis, motif analysis is available in rmarkdown format under the [code](multiomics/code/real_data_B_multiomics_20230404.rmd) section of multiomics folder. The corresponding html document with results and plots, and final seurat object with corresponding metadata can be downloaded from Mendeley Dataset https://data.mendeley.com/datasets/fjxrcbh672/draft?a=69394d54-235c-436e-be60-520cd2899517.

Genecount matrices and ATAC peaks matrices for each sample were generated using cellranger (details of parameters available in *Methods* section of the manuscript. The results were merged again to get aggregated counts and peak matrices for all samples. The [script](multiomics/data/run_merge_brain.sh) to merge the sample and the corresponding [input csv](multiomics/data/libraries_brain.csv) file with input paths are under the data section of multiomics folder.

### Deconvolution
ST spot decomposition was performed using Stereoscope. The steps include:
1. Prepare data: First, create h5ad objects of both [ST](prep_st_brain_220512_for_deconv.Rmd) and [SN](deconvolution-stereoscope/code/prep_sc_data_for_deconv_brain_220609.Rmd) data. For SN data, subsample by cluster to similar number of cells per cluster and run DEG detection. Select a geneset for deconvolution. Was done with either top20, top50 or top100 genes per cluster, or by taking top 2000 variable genes. Also, only select genes that are expressed in the ST data as well.
2. Run deconvolution: [Bash](deconvolution-stereoscope/code/run_stereoscope_brain_220609.bash) script (have only changed the line with SC_LIST to define different gene sets)
3. Summarize results: Finally, summarize the results using this [script](deconvolution-stereoscope/code/summary_stereoscope_brain_220609.Rmd).

The deconvolution results saved as an R object can be downloaded from the Mendeley Dataset https://data.mendeley.com/datasets/fjxrcbh672/draft?a=69394d54-235c-436e-be60-520cd2899517.

### Ligand Receptor interactions analysis
L-R interactions were annalysed via this [code](L-R_interactions/code/L_R_Brain.Rmd). Resultant plots are shown and discussed in the manuscript.

### Celltype interactions analysis using Misty and Progeny
More details on this analysis can be found in our manuscript and all the code required to run this analysis is available as a Snakemake pipeline in the github repository: **https://github.com/saezlab/astromouse**

### Metabolic analysis
This collection of programs generates pathway enrichment analyses for metabolic collections of genes featured in the Recon3D model and found under [metabolic-analysis](metabolic_analysis/code)

Recon3D is described in *Brunk E, Sahoo S, Zielinski DC, Altunkaya A, Dräger A, Mih N, Gatto F, Nilsson A, Preciat Gonzalez GA, Aurich MK, Prlić A, Sastry A, Danielsdottir AD, Heinken A, Noronha A, Rose PW, Burley SK, Fleming RMT, Nielsen J, Thiele I, Palsson BO. Recon3D enables a three-dimensional view of gene variation in human metabolism. Nat Biotechnol. 2018 Mar;36(3):272-281. doi: 10.1038/nbt.4072. Epub 2018 Feb 19. PMID: 29457794; PMCID: PMC5840010.*

#### Process to generate gene enrichment analyses:

1. Extract out the metabolic subsystem ("pathway") memberships of all genes in Recon3D human metabolic model. This is done with the programs and files in the folder, [RECON3D_Hs_genes_subsystems_to_Mm](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm)
	a) Download (or use our copy) of [Recon3D.json](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D.json) from the BIGG models database: https://bigg.ucsd.edu/models/Recon3D
	b) Use [parse_recon3d_for_id_subsystem.py](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/parse_recon3d_for_id_subsystem.py) on Recon3D.json to generate the reaction-to-subsystem file, [Recon3D_id_subsystem.txt](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D_id_subsystem.txt)
	c) Use [parse_recon3d.py](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/parse_recon3d.py) to generate a list of subsystems (pathways) and associated genes as [Recon3D_pathways_gene_symbol_tbl.csv](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D_pathways_gene_symbol_tbl.csv)
2. Translate the human Recon3D genes to mouse IDs using the Hugo Gene Nomenclature Committee (HGNC) HCOP (HGNC Comparison of Orthology Predictions) from https://www.genenames.org/tools/hcop/ -- we used a very specific structure of HCOP included here, in the file [HCOP_mouse_human_15col.txt](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/HCOP_mouse_human_15col.txt)

	a) In the directory, [GSEA_with_Recon3D_GMT_files](metabolic_analysis/code/GSEA_with_Recon3D_GMT_files) we have files for creating GMT files. Use [replace_human_with_mouse_ids.py](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/replace_human_with_mouse_ids.py) to use files [Recon3D_pathways_gene_symbol_tbl.csv](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D_pathways_gene_symbol_tbl.csv) and [HCOP_mouse_human_15col.txt](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/HCOP_mouse_human_15col.txt) to output a mouse-translated form of the Recon3D pathways, output as "Recon3D_pathway_Mm_gene_symbol_tbl.tab". You MUST double-check the translations as a new download of HCOP does not necessarily do a good job in translating certain genes from human to mouse. The Oxidative Phosphorylation pathway is especially affected. Manual review is always necessary. Our review of this output with corrected OxPhos genes is given here as [Recon3D_pathway_Mm_gene_symbol_withATP_tbl.tab](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D_pathway_Mm_gene_symbol_withATP_tbl.tab)
3. Make the GMT files from the human and mouse gene-to-subsystem files for downstream use.
	a) The python program  [make_mouse_gmt.py](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/make_mouse_gmt.py) is used with our manual review output file from (2): [Recon3D_pathway_Mm_gene_symbol_withATP_tbl.tab](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D_pathway_Mm_gene_symbol_withATP_tbl.tab) to generate a .gmt file.
4. Analysis was executed in R (RStudio) and all files and output are in the directory [R_analyses_from_Seurat_4_output](metabolic_analysis/code/R_analyses_from_Seurat_4_output). Analyses were done using the R Markdown File [brain_Seurat_4_FGSEA_analysis.Rmd](metabolic_analysis/code/R_analyses_from_Seurat_4_output/brain_Seurat_4_FGSEA_analysis.Rmd) which is well-annotated. All steps above generate the [Recon3D_pathway_Mm_gene_symbol_withATP_tbl.tab](metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D_pathway_Mm_gene_symbol_withATP_tbl.tab) file used in the analysis. The first part of the script is annotated as configurable and is very specific to this dataset. This particular analysis just calculates the pathway enrichment from the Recon3D gmt files created in steps 1-3.
This script also generates the figures from the paper.

The input files and output data are available via Mendeley dataset https://data.mendeley.com/datasets/fjxrcbh672/draft?a=69394d54-235c-436e-be60-520cd2899517 as also mentioned in the 'Data availability' section of the paper. 

### Metabolic analysis
This collection of programs generates pathway enrichment analyses for metabolic collections of genes featured in the Recon3D model and found under [metabolic-analysis] (metabolic_analysis/code)

Recon3D is described in *Brunk E, Sahoo S, Zielinski DC, Altunkaya A, Dräger A, Mih N, Gatto F, Nilsson A, Preciat Gonzalez GA, Aurich MK, Prlić A, Sastry A, Danielsdottir AD, Heinken A, Noronha A, Rose PW, Burley SK, Fleming RMT, Nielsen J, Thiele I, Palsson BO. Recon3D enables a three-dimensional view of gene variation in human metabolism. Nat Biotechnol. 2018 Mar;36(3):272-281. doi: 10.1038/nbt.4072. Epub 2018 Feb 19. PMID: 29457794; PMCID: PMC5840010.*

#### Process to generate gene enrichment analyses:

1. Extract out the metabolic subsystem ("pathway") memberships of all genes in Recon3D human metabolic model. This is done with the programs and files in the folder, [RECON3D_Hs_genes_subsystems_to_Mm] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm)
	a) Download (or use our copy) of [Recon3D.json] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D.json) from the BIGG models database: https://bigg.ucsd.edu/models/Recon3D
	b) Use [parse_recon3d_for_id_subsystem.py] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/parse_recon3d_for_id_subsystem.py) on Recon3D.json to generate the reaction-to-subsystem file, [Recon3D_id_subsystem.txt] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D_id_subsystem.txt)
	c) Use [parse_recon3d.py] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/parse_recon3d.py) to generate a list of subsystems (pathways) and associated genes as [Recon3D_pathways_gene_symbol_tbl.csv] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D_pathways_gene_symbol_tbl.csv)
2. Translate the human Recon3D genes to mouse IDs using the Hugo Gene Nomenclature Committee (HGNC) HCOP (HGNC Comparison of Orthology Predictions) from https://www.genenames.org/tools/hcop/ -- we used a very specific structure of HCOP included here, in the file [HCOP_mouse_human_15col.txt] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/HCOP_mouse_human_15col.txt)

	a) In the directory, [GSEA_with_Recon3D_GMT_files] (metabolic_analysis/code/GSEA_with_Recon3D_GMT_files) we have files for creating GMT files. Use [replace_human_with_mouse_ids.py] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/replace_human_with_mouse_ids.py) to use files [Recon3D_pathways_gene_symbol_tbl.csv] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D_pathways_gene_symbol_tbl.csv) and [HCOP_mouse_human_15col.txt] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/HCOP_mouse_human_15col.txt) to output a mouse-translated form of the Recon3D pathways, output as "Recon3D_pathway_Mm_gene_symbol_tbl.tab". You MUST double-check the translations as a new download of HCOP does not necessarily do a good job in translating certain genes from human to mouse. The Oxidative Phosphorylation pathway is especially affected. Manual review is always necessary. Our review of this output with corrected OxPhos genes is given here as [Recon3D_pathway_Mm_gene_symbol_withATP_tbl.tab] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D_pathway_Mm_gene_symbol_withATP_tbl.tab)
3. Make the GMT files from the human and mouse gene-to-subsystem files for downstream use.
	a) The python program  [make_mouse_gmt.py] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/make_mouse_gmt.py) is used with our manual review output file from (2): [Recon3D_pathway_Mm_gene_symbol_withATP_tbl.tab] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D_pathway_Mm_gene_symbol_withATP_tbl.tab) to generate a .gmt file.
4. Analysis was executed in R (RStudio) and all files and output are in the directory [R_analyses_from_Seurat_4_output](metabolic_analysis/code/R_analyses_from_Seurat_4_output). Analyses were done using the R Markdown File [brain_Seurat_4_FGSEA_analysis.Rmd] (metabolic_analysis/code/R_analyses_from_Seurat_4_output/brain_Seurat_4_FGSEA_analysis.Rmd) which is well-annotated. All steps above generate the [Recon3D_pathway_Mm_gene_symbol_withATP_tbl.tab] (metabolic_analysis/code/RECON3D_Hs_genes_subsystems_to_Mm/Recon3D_pathway_Mm_gene_symbol_withATP_tbl.tab) file used in the analysis. The first part of the script is annotated as configurable and is very specific to this dataset. This particular analysis just calculates the pathway enrichment from the Recon3D gmt files created in steps 1-3.
This script also generates the figures from the paper.

The input files and output data are available via Mendeley dataset link of which is available in the 'Data availability' section of the paper. 

### Validation_RNAscope_quantification
This section contains the code and data for the validation experiments performed in this study. The code folder contains the code to generate the two boxplots that show the RNAscope quantification for two genes of interest (Adcy1 and Gpc5) in this study. The data folder contains the csv files with quantified RNAscope signal for each gene as a separate folder. Each csv filename starts with the samplename used in the study.
