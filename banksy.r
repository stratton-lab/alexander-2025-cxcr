options(future.globals.maxSize = 10^12)
suppressPackageStartupMessages({
  library(SeuratObject)
  library(Seurat)
  library(Matrix)
  library(SeuratWrappers)
  library(purrr)
  library(reticulate)
})

from_pseudospot <- function(ad.path, id) {
  message(sprintf("%s: loading anndata", id))
  ad <- import("anndata", convert = FALSE)$read_h5ad(ad.path)

  message(sprintf("%s: creating seurat object", id))
  so <- py_to_r(ad$T$X) |>
    `dimnames<-`(list(
      py_to_r(ad$var_names$to_list()),
      py_to_r(ad$obs_names$to_list())
    )) |>
    CreateSeuratObject(
      assay = "MERSCOPE",
      project = id
    )
  message(sprintf("%s: loading blanks assay", id))
  so[["blanks"]] <- ad$obsm["blanks"] |>
    py_to_r() |>
    t() |>
    Matrix() |>
    CreateAssay5Object()

  so <- AddMetaData(so, py_to_r(ad$obs))

  so
}

# patch get_locs function to properly add z coordinates
# adapted from: https://github.com/satijalab/seurat-wrappers/blob/a1eb0d8b039ad6d5ef0ff1332fd8eb1c0c223553/R/banksy.R#L159
assignInNamespace(
  "get_locs",
  function(object, dimx, dimy, dimz, ndim, data_own, group, verbose) {
    if (!is.null(dimx) & !is.null(dimy)) {
      # Extract locations from metadata
      locs <- data.frame(
        sdimx = unlist(object[[dimx]]),
        sdimy = unlist(object[[dimy]])
      )
      rownames(locs) <- colnames(object)

      # Add z-dim if present
      # CHANGE FROM SOURCE: `object[[dimz]]` to `unlist(object[[dimz]])`
      if (!is.null(dimz)) locs[["sdimz"]] <- unlist(object[[dimz]])

      # Check locations
      obj_samples <- colnames(data_own)
      locs_samples <- rownames(locs)
      if (any(is.na(match(obj_samples, locs_samples)))) {
        na_id <- which(is.na(match(obj_samples, locs_samples)))
        warning(
          "No centroids found for samples: ",
          paste(obj_samples[na_id], collapse = ", "), ". Dropping samples."
        )
        data_own <- data_own[, -na_id, drop = FALSE]
      }
      locs <- locs[match(obj_samples, locs_samples), , drop = FALSE]
    } else {
      # Extract locations with Seurat accessor
      locs <- Seurat::GetTissueCoordinates(object)[, seq_len(ndim)]
    }

    dim_names <- paste0("sdim", c("x", "y", "z"))
    colnames(locs) <- dim_names[seq_len(ncol(locs))]

    if (!is.null(group)) {
      # Stagger locations by group
      if (verbose) message("Staggering locations by ", group)
      locs[, 1] <- locs[, 1] + abs(min(locs[, 1]))
      max_x <- max(locs[, 1]) * 2
      n_groups <- length(unique(unlist(object[[group]])))
      shift <- seq(from = 0, length.out = n_groups, by = max_x)
      locs[, 1] <- locs[, 1] + rep(shift, table(object[[group]]))
    }

    locs
  },
  "SeuratWrappers"
)

so.list <- c("s5r0", "s5r1", "s11r0", "s11r1") |>
  map(function(ident) {
    sprintf("data/ad/%s.h5ad", ident) |>
      from_pseudospot(ident) |>
      subset(nCount_MERSCOPE >= 8) |>
      NormalizeData() |>
      AddMetaData(ident, "orig.ident")
  })

so <- JoinLayers(merge(so.list[[1]], so.list[[2:length(so.list)]]))

so <- RunBanksy(so,
  lambda = 0.2,
  assay = "MERSCOPE", slot = "data",
  dimx = "spot_x", dimy = "spot_y", dimz = "spot_z", ndim = 3,
  features = "all", group = "orig.ident", split.scale = TRUE
)

so <- AddMetaData(so, colSums(GetAssayData(so)), "nCount_BANKSY")
so <- AddMetaData(so, colSums(GetAssayData(so) > 0), "nFeature_BANKSY")

so <- subset(so, nCount_BANKSY > 70)
VariableFeatures(so) <- Features(so)

so <- so |>
  RunPCA() |>
  FindNeighbors(dims = 1:15) |>
  FindClusters(algorithm = "leiden")

saveRDS(so, "data/rds/all_banksy.seurat.rds")

iwalk(
  SplitObject(so, "orig.ident"),
  function(so, id) {
    so |>
      RunPCA() |>
      FindNeighbors(dims = 1:15) |>
      RunUMAP(dims = 1:15) |>
      saveRDS(sprintf("data/rds/%s_banksy.seurat.rds", id))
  }
)
