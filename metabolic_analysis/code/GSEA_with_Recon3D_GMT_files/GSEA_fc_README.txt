
to calculate the FC in files in flight_vs_ground_exp_fc_tbls_Nov_2022:

FindMarkers(
  xseu,
  ident.1 = "flight", ident.2 =  "ground",
  group.by = condition_col, subset.ident = xclust,
  logfc.threshold = 0,  # No filter on logfc, because low logfc genes are also informative in GSEA
  min.pct = 0.01  # Remove genes that are expressed in < 1% of cells in either group, to avoid division by 0
)
