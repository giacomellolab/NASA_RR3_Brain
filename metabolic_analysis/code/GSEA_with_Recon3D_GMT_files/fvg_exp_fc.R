library(Signac)
library(Seurat)
library(tidyverse)


# Ref: https://stackoverflow.com/a/5577647/4638182
load_multiomics_rdata <- function(file_path) {
  mo_env <- new.env()
  mo_obj_names <- load(file_path, mo_env)
  stopifnot(identical(mo_obj_names, c("seurat_multi")))
  mo_seu <- mo_env[["seurat_multi"]]

  return(mo_seu)
}

sid_seu_list <- list(
  mo_brain = load_multiomics_rdata(
    "data/sc_rna/seurat_data_brain_cc_2022-05-17.RData"),

  st_brain = readRDS(
    "data/sc_rna/220512_RR3_brain_SCTpersection_harmony_annot.rds")
)


iwalk(sid_seu_list, function(xseu, xname) {
  print(xname)
  print(colSums((xseu@assays$RNA@data[, 1:3])))
  print(colSums((xseu@assays$RNA@counts[, 1:3])))
})

pp_sid_seu_list <- map(sid_seu_list, function(xseu) {
  DefaultAssay(xseu) <- "RNA"

  xseu <- NormalizeData(
    xseu, normalization.method = "LogNormalize", scale.factor = 10000,
    verbose = FALSE)

  return(xseu)
})

iwalk(pp_sid_seu_list, function(xseu, xname) {
  print(xname)
  print(colSums((xseu@assays$RNA@data[, 1:3])))
  print(colSums(expm1(xseu@assays$RNA@data[, 1:3])))
  print(colSums((xseu@assays$RNA@counts[, 1:3])))
})



fvg_efc_tbls_list <- imap(pp_sid_seu_list, function(xseu, xname) {
  if (xname == "mo_brain") {
    condition_col <- "cell_flight_status"

  } else if (xname == "st_brain") {
    condition_col <- "sample_condition"

  } else {
    stop(paste("Uknown sample id", xname))
  }

  stopifnot(condition_col %in% colnames(xseu@meta.data))

  stopifnot(identical(sum(is.na(xseu@meta.data[, condition_col])), 0L))
  stopifnot(is.character(xseu@meta.data[, condition_col]))

  stopifnot(identical(
    sort(unique(xseu@meta.data[, condition_col])),
    c("flight",  "ground")
  ))

  stopifnot(identical(
    rownames(xseu@meta.data),
    names(Idents(xseu))
  ))

  print(xname)
  print(
    table(
      data.frame(
        cluster = as.character(Idents(xseu)),
        condition = xseu@meta.data[, condition_col])
    )
  )

  fvg_efc_tbls <- map(unique(as.character(Idents(xseu))), function(xclust) {
    print(xclust)
    x_ss_seu <- subset(xseu, idents = c(xclust))
    if (sum(x_ss_seu@meta.data[, condition_col] == "flight") < 3) {
      print(paste(xname, "cluster", xclust, "has < 3 flight cells."))
      return(NULL)
    }

    if (sum(x_ss_seu@meta.data[, condition_col] ==  "ground") < 3) {
      print(paste(xname, "cluster", xclust, "has < 3 ground cells."))
      return(NULL)
    }

    x_dge_tbl <- FindMarkers(
      xseu,
      ident.1 = "flight", ident.2 =  "ground",
      group.by = condition_col, subset.ident = xclust,
      logfc.threshold = 0, min.pct = 0.01)

    res <- list(
      cluster = xclust,
      fvg_efc_tbl = x_dge_tbl
    )

    return(res)
  })

  return(fvg_efc_tbls)
})

out_res_dir <- "results/flight_vs_ground_exp_fc_tbls_Nov_2022"
if (!dir.exists(out_res_dir)) {
  dir.create(out_res_dir)
}

iwalk(fvg_efc_tbls_list, function(xl, x_sample) {
  print(x_sample)
  walk(xl, function(x_clust_list) {
    print(x_clust_list$cluster)
    print(
      sum(is.na(rownames_to_column(x_clust_list$fvg_efc_tbl, "gene_symbol"))))
  })
})

iwalk(fvg_efc_tbls_list, function(xl, x_sample) {
  print(x_sample)

  walk(xl, function(x_clust_list) {
    print(x_clust_list$cluster)

    stopifnot(identical(
      sum(is.na(rownames_to_column(x_clust_list$fvg_efc_tbl, "gene_symbol"))),
      0L
    ))

    write_tsv(
      rownames_to_column(x_clust_list$fvg_efc_tbl, "gene_symbol") %>%
        arrange(desc(avg_log2FC)) %>%
        mutate(sample_id = x_sample, cluster = x_clust_list$cluster),
      file.path(
        out_res_dir,
        paste0(x_sample, "_cluster", x_clust_list$cluster,
        "_flight_vs_ground_log_norm_count_wilcox_tbl.tsv"))
    )

  })
})
