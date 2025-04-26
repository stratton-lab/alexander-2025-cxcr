suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SeuratWrappers)
  library(Matrix)
  library(reticulate)
  library(ggplot2)
  library(data.table)
  library(purrr)
  library(ggridges)
  library(patchwork)
  library(glmGamPoi)
  library(forcats)
  library(ggrepel)
})

gcs <- function(clusts, meta.dt) {
  map(c("s5r0", "s5r1", "s11r0", "s11r1"), function(i) {
    so <- readRDS(sprintf("data/rds/%s_banksy.seurat.rds", i))

    so <- so[
      ,
      meta.dt[Cells(so), seurat_clusters %in% clusts, on = .(cell)]
    ]

    so <- AddMetaData(
      so,
      meta.dt[
        Cells(so),
        .(cell, condition, mapped.ident),
        on = .(cell)
      ] |> data.frame(row.names = "cell")
    )

    so
  }) |>
    purrr::reduce(merge) |>
    `DefaultAssay<-`(value = "MERSCOPE") |>
    JoinLayers()
}

dbc <- function(so, drop.features = NULL) {
  mtx <- as.matrix(GetAssayData(so, assay = "MERSCOPE", layer = "counts"))
  if (!is.null(drop.features)) {
    mtx <- mtx[!(rownames(mtx) %in% unique(drop.features)), ]
  }
  pseudobulk(
    mtx,
    vars(mapped.ident, condition),
    col_data = so[[]]
  ) |>
    glm_gp(design = ~condition) |>
    test_de(cond(condition = "inf") - cond(condition = "ctl")) |>
    as.data.table()
}


meta.all <- as.data.table(
  readRDS("data/rds/all_banksy.seurat.rds")[[]],
  keep.rownames = "cell"
)[, `:=`(
  condition = fifelse(
    endsWith(orig.ident, "r0"),
    "inf",
    "ctl"
  ) |> fct_relevel("ctl", "inf"),
  mapped.ident = fcase(
    orig.ident == "s5r1", "sham:1",
    orig.ident == "s11r1", "sham:2",
    orig.ident == "s5r0", "inf:1",
    orig.ident == "s11r0", "inf:2"
  ) |> fct_relevel("sham:1", "sham:2", "inf:1", "inf:2")
)][]

data.sel <- rbindlist(map(c("s5r0", "s5r1", "s11r0", "s11r1"), function(i) {
  so <- readRDS(sprintf("data/rds/%s_banksy.seurat.rds", i))

  as.data.table(
    FetchData(
      AddMetaData(
        so,
        meta.all[
          Cells(so),
          .(cell, merged.clust = seurat_clusters, condition, mapped.ident),
          on = .(cell)
        ] |> data.frame(row.names = "cell")
      ),
      c(
        "condition", "merged.clust", "mapped.ident", # metadata
        "LTag MuPyV", "VP1 MuPyV", # viral expression
        "Mbp", "Mobp", # oligo
        "Pdgfra", "Pcdh15", # OPC
        "Csf1r", "Cx3cr1", # myeloid
        "Cd8a", "Cd8b1", # CD8 T cell
        "Foxj1", "Cfap65", # ependyma
        "Cxcl12", "Cxcl16", "Cxcr4", "Cxcr6", # genes of interest
        "spot_x", "spot_y", "spot_z" # spatial
      )
    ),
    keep.rownames = "cell"
  )
}))[
  data.table(
    mapped.ident = c("inf:1", "sham:1", "inf:2", "sham:2"),
    rot.theta = c(-80, -80, 180, 180) * pi / 180
  ),
  `:=`(
    x = spot_x * cos(rot.theta) - spot_y * sin(rot.theta),
    y = spot_x * sin(rot.theta) + spot_y * cos(rot.theta)
  ),
  on = .(mapped.ident)
][]

tcd8.so <- gcs(c(24), meta.all)
ep.so <- gcs(c(28), meta.all)
ol.so <- gcs(c(23), meta.all)
opc.so <- gcs(c(25), meta.all)
myeloid.so <- gcs(c(42, 17, 27, 3, 34), meta.all)

ep.so.s <- ep.so[
  ,
  rbindlist(
    map(
      c("s5r0", "s5r1", "s11r0", "s11r1"),
      function(x) {
        fread(sprintf("data/scs/%s.tsv", x))
      }
    )
  )[within_bb == TRUE][
    !rbindlist(
      map(
        c("s5r0_choroid", "s5r1_choroid"),
        function(x) {
          fread(sprintf("data/scs/%s.tsv", x))
        }
      )
    )[within_bb == TRUE],
    cell,
    on = .(cell)
  ]
] |>
  subset(
    Cd8a < 1 & Cd8b1 < 1 &
      Csf1r < 1 & Cx3cr1 < 1 &
      Foxj1 > 3 & Cfap65 > 3
  )

tcd8.so.s <- subset(
  tcd8.so,
  Foxj1 < 1 & Cfap65 < 1 &
    Csf1r < 1 & Cx3cr1 < 1 &
    Cd8a > 3 & Cd8b1 > 3
)

myeloid.so.s <- subset(
  myeloid.so,
  Foxj1 < 1 & Cfap65 < 1 &
    Cd8a < 1 & Cd8b1 < 1 &
    Csf1r > 3 & Cx3cr1 > 3
)


data.sel[, `:=`(cell_type = as.character(NA))]
data.sel[
  Cells(ep.so.s), `:=`(cell_type = "ependyma"),
  on = .(cell)
][
  Cells(tcd8.so.s), `:=`(cell_type = "Cd8 T"),
  on = .(cell)
][
  Cells(myeloid.so.s), `:=`(cell_type = "myeloid"),
  on = .(cell)
][
  Cells(ol.so), `:=`(cell_type = "oligo"),
  on = .(cell)
][
  Cells(opc.so), `:=`(cell_type = "OPC"),
  on = .(cell)
]

fwrite(data.sel, "out/data.tsv", sep = "\t")

hec <- function(so, threshold) {
  mtx <- GetAssayData(so)
  exp <- rowSums(mtx > 0) / ncol(mtx)

  names(exp[exp > threshold])
}

deg.list <- list(
  "ependyma" = dbc(
    ep.so.s, c(
      hec(subset(tcd8.so.s, condition == "inf"), 0.35),
      hec(subset(myeloid.so.s, condition == "inf"), 0.35)
    )
  ),
  "oligo" = dbc(ol.so),
  "opc" = dbc(opc.so),
  "myeloid" = dbc(myeloid.so.s, c(
    hec(subset(tcd8.so.s, condition == "inf"), 0.35),
    hec(subset(ep.so.s, condition == "inf"), 0.35)
  ))
)

iwalk(list(
  "ependyma" = ep.so.s,
  "tcd8" = tcd8.so.s,
  "oligo" = ol.so,
  "opc" = opc.so,
  "myeloid" = myeloid.so.s
), function(x, i) {
  saveRDS(x, sprintf("out/rds/%s.rds", i))
})

iwalk(deg.list, function(x, i) {
  fwrite(x, sprintf("out/deg/%s.tsv", i), sep = "\t")
})
