suppressPackageStartupMessages({
  library(ggplot2)
  library(SeuratObject)
  library(forcats)
  library(ggrepel)
  library(ggridges)
  library(patchwork)
  library(data.table)
  library(purrr)
})

vlc <- function(deg, pval.thresh = 0.05, abs.lfc.thresh = 1) {
  ggplot(deg) +
    aes(
      x = lfc,
      y = -log10(adj_pval),
      label = fifelse(adj_pval < pval.thresh & abs(lfc) > abs.lfc.thresh, name, NA),
      colour = fct_relevel(
        fcase(
          adj_pval >= 0.05, "ns",
          lfc > 0, "up",
          lfc < 0, "down"
        ),
        "up", "down", "ns"
      )
    ) +
    geom_point() +
    geom_label_repel(show.legend = FALSE) +
    labs(
      x = "log fold change",
      y = expression(-log[10](adj_pval)),
      colour = NULL,
    )
}

sep <- function(gene, DT, filter_type = NA, lims = NULL) {
  ggplot(
    if (!is.na(filter_type)) {
      DT[order(fcoalesce(cell_type == filter_type, FALSE), decreasing = TRUE)]
    } else {
      DT
    }
  ) +
    aes(
      x = x,
      y = y,
      colour = if (!is.na(filter_type)) {
        fifelse(fcoalesce(cell_type == filter_type, FALSE), get(gene), NA)
      } else {
        get(gene)
      },
      alpha = if (!is.na(filter_type)) {
        fifelse(fcoalesce(cell_type == filter_type, FALSE), 1, 0.006)
      } else {
        1
      }
    ) +
    geom_point(shape = ".") +
    scale_alpha_identity() +
    theme_void() +
    scale_colour_viridis_c(
      option = "cividis",
      limits = lims
    ) +
    coord_fixed() +
    theme(
      plot.background = element_rect(fill = "black", colour = NA),
      panel.background = element_rect(fill = "black", colour = NA),
      text = element_text(colour = "white"),
      plot.margin = margin()
    ) +
    labs(colour = gene)
}
gss <- function(s.dt, ct, padding = 100) {
  vls <- s.dt[cell_type == ct, c(
    (max(max(x) - min(x), max(y) - min(y)) / 2) + padding,
    (max(x) + min(x)) / 2,
    (max(y) + min(y)) / 2
  )]
  s.dt[
    x %between% c(vls[[2]] - vls[[1]], vls[[2]] + vls[[1]]) &
      y %between% c(vls[[3]] - vls[[1]], vls[[3]] + vls[[1]])
  ]
}

gce <- function(x, y, DT) {
  ggplot() +
    aes(x = get(x), y = get(y)) +
    geom_text(
      aes(label = merged.clust, colour = condition),
      DT[
        , lapply(.SD, mean),
        .SDcols = c(x, y), by = .(merged.clust, condition)
      ],
      key_glyph = draw_key_point
    ) +
    scale_fill_viridis_c(option = "cividis") +
    labs(x = x, y = y, fill = NULL)
}

main.fig <- function(DT, tcd8.so.s, ep.deg, interaction.dt) {
  split.by.id <- split(
    DT[order(fct_relevel(
      mapped.ident,
      "sham:1", "sham:2", "inf:1", "inf:2"
    ))],
    by = "mapped.ident", drop = TRUE
  )

  ltag.pl <- split.by.id |>
    imap(function(x, i) {
      sep(
        "LTag MuPyV", x,
        lims = c(0, DT[, max(`LTag MuPyV`)])
      ) + labs(title = i)
    }) |>
    wrap_plots(... = _, widths = 1, heights = 1, guides = "collect") &
    theme(plot.margin = margin())

  cluster.ids <- wrap_plots(
    gce("Mbp", "Mobp", DT) + labs(title = "oligodendrocyte"),
    gce("Pdgfra", "Pcdh15", DT) + labs(title = "OPC"),
    gce("Csf1r", "Cx3cr1", DT) + labs(title = "myeloid"),
    gce("Cd8a", "Cd8b1", DT) + labs(title = "Cd8 T cell"),
    gce("Foxj1", "Cfap65", DT) + labs(title = "ependymal"),
    gce("LTag MuPyV", "VP1 MuPyV", DT) + labs(title = "MuPyV"),
    guides = "collect"
  )

  ep.vlc <- vlc(ep.deg)

  cxcl.pl <- map(c("Cxcl12", "Cxcl16"), function(gene) {
    split.by.id |>
      imap(function(s.dt, i) {
        wrap_plots(
          sep(
            gene, s.dt,
            "ependyma", c(0, DT[cell_type == "ependyma", max(get(gene))])
          ) + labs(title = i),
          sep(
            gene, gss(s.dt, "ependyma"),
            "ependyma", c(0, DT[cell_type == "ependyma", max(get(gene))])
          ),
          design = "A####\n#BBBB\n#BBBB\n#BBBB"
        )
      }) |>
      wrap_plots(... = _, widths = 1, heights = 1, guides = "collect")
  }) |>
    wrap_plots(... = _, design = "AB", widths = 1) &
    plot_annotation(
      theme = theme(plot.background = element_rect(fill = "black"))
    ) &
    theme(plot.margin = margin())

  loi <- c(
    "Ccr1", "Ccr2", "Ccr3", "Ccr5", "Ccr6",
    "Cxcr1", "Cxcr2", "Cxcr3", "Cxcr4", "Cxcr5", "Cxcr6",
    "Ifnar1", "Ifngr1", "Ifnlr1",
    "Il23r", "Il12rb2",
    "Ifng", "Gzmb"
  )
  cd8.data <- subset(tcd8.so.s, condition == "inf") |>
    FetchData(c(loi, "mapped.ident")) |>
    as.data.table(keep.rownames = "cell") |>
    melt(
      id.vars = c("cell", "mapped.ident"),
      value.name = "expr",
      variable.name = "gene"
    ) |>
    _[, `:=`(expr.avg = mean(expr)), by = .(gene, mapped.ident)][]
  cd8.pl <- ggplot() +
    aes(x = expr, y = fct_rev(gene), fill = mapped.ident) +
    geom_ridgeline(
      aes(height = after_stat(scaled)),
      cd8.data[mapped.ident == "inf:1"],
      stat = "density", trim = TRUE, scale = 0.4
    ) +
    geom_ridgeline(
      aes(height = -after_stat(scaled)),
      cd8.data[mapped.ident == "inf:2"],
      stat = "density", trim = TRUE, scale = 0.4, min_height = -1
    ) +
    labs(x = "normalized expression", y = NULL, fill = "identifier") +
    guides(fill = guide_legend(override.aes = list(alpha = 1)))

  cxcr.pl <- map(c("Cxcr4", "Cxcr6"), function(gene) {
    wrap_elements(
      split(DT[condition == "inf"], by = "mapped.ident", drop = TRUE) |>
        imap(function(x, i) {
          sep(gene, x, "Cd8 T") + labs(title = i)
        }) |>
        wrap_plots(... = _, widths = 1, heights = 1, guides = "collect") &
        theme(plot.margin = margin())
    )
  }) |>
    wrap_plots(... = _, design = "A\nB", widths = 1) &
    theme(plot.margin = margin())

  tcd8.mtx <- GetAssayData(subset(tcd8.so.s, condition == "inf"))
  tcd8.exp <- rowSums(tcd8.mtx > 0) / ncol(tcd8.mtx)
  lr.pl <- interaction.dt[
    ep.deg[adj_pval < 0.05 & lfc > 0, name],
    on = .(from),
    nomatch = NULL
  ][
    names(tcd8.exp[tcd8.exp > 0.1]),
    on = .(to),
    nomatch = NULL
  ] |>
    dcast(from ~ to, value.var = "weight", fill = 0) |>
    melt(id.vars = "from", variable.name = "to", value.name = "interaction") |>
    ggplot() +
    aes(
      x = fct_reorder(to, interaction, mean, .desc = TRUE),
      y = fct_reorder(from, interaction, mean),
      fill = interaction
    ) +
    geom_tile() +
    scale_fill_viridis_c(option = "turbo") +
    scale_x_discrete(position = "top") +
    labs(
      x = "receptors present in Cd8 T cells following infection",
      y = "ligands upregulated in ependyma following infection",
      fill = "interaction\npotential"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0))


  wrap_plots(
    ltag.pl,
    cluster.ids,
    free(ep.vlc),
    free(wrap_elements(cxcl.pl)),
    free(cd8.pl),
    free(cxcr.pl),
    lr.pl,
    # blanks added added to provide
    # space for figure annotations
    design = "
      #AAAAAAAAAAAA#BBBBBBBBBBBBB
      #CCCCCCCCCC#DDDDDDDDDDDDDDD
      #EEEEEEEE#FFFFFFFF#GGGGGGGG
    "
  )
}

supp.fig <- function(ol.deg, opc.deg, mu.deg, DT) {
  append(
    imap(
      list(oligo = ol.deg, OPC = opc.deg, myeloid = mu.deg),
      function(x, i) {
        vlc(x, abs.lfc.thresh = 0) + labs(title = i)
      }
    ),
    map(c("oligo", "OPC", "myeloid"), function(ct) {
      map(c("Cxcl12", "Cxcl16"), function(gene) {
        wrap_elements(
          split(
            DT[order(fct_relevel(
              mapped.ident,
              "sham:1", "sham:2", "inf:1", "inf:2"
            ))],
            by = "mapped.ident", drop = TRUE
          ) |>
            imap(function(x, i) {
              sep(gene, x, ct, c(0, DT[, max(get(gene))])) + labs(title = i)
            }) |>
            wrap_plots(... = _, widths = 1, heights = 1, guides = "collect") &
            theme(plot.margin = margin())
        )
      }) |>
        wrap_plots(... = _, design = "AB", widths = 1)
    })
  ) |>
    wrap_plots(
      ... = _,
      design = "
        AADDDD
        BBEEEE
        CCFFFF
      "
    ) & theme(plot.margin = margin())
}

# load data
DT <- fread("out/data.tsv")[, `:=`(
  condition = fct_relevel(
    fifelse(condition == "ctl", "sham", condition),
    "sham", "inf"
  )
)][]
so.list <- map(
  setNames(nm = c(
    "ependyma",
    "oligo",
    "opc",
    "myeloid",
    "tcd8"
  )), function(i) {
    readRDS(sprintf("out/rds/%s.rds", i))
  }
)
deg.list <- map(
  setNames(nm = c(
    "ependyma",
    "oligo",
    "opc",
    "myeloid"
  )), function(i) {
    fread(sprintf("out/deg/%s.tsv", i))
  }
)
nn.dt <- if (file.exists("data/nn/wn.rds")) {
  as.data.table(readRDS("data/nn/wn.rds")[["lr_sig"]])
} else {
  wn <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  saveRDS(wn, "data/nn/wn.rds")
  as.data.table(wn[["lr_sig"]])
}


# generate figures
main <- main.fig(
  DT, so.list[["tcd8"]], deg.list[["ependyma"]], nn.dt
) & theme(text = element_text(size = 14))
supp <- supp.fig(
  deg.list[["oligo"]], deg.list[["opc"]], deg.list[["myeloid"]], DT
) & theme(text = element_text(size = 16))

# utility function to save plots
dff <- function(path, plot, w, h, rec = FALSE, png = TRUE, pdf = TRUE, rds = TRUE) {
  if (png) {
    ggsave(
      sprintf("%s/fig.png", path), plot,
      w = w, h = h, units = "px",
      create.dir = TRUE
    )
  }
  if (pdf) {
    ggsave(sprintf("%s/fig.pdf", path), plot, w = w, h = h, units = "px")
  }
  if (rds) {
    saveRDS(plot, sprintf("%s/fig.rds", path))
  }

  if (
    ("patchwork" %in% class(plot)) &&
      rec > 0 &&
      !is.null(plot$patches[["layout"]][["design"]])
  ) {
    ncol <- max(plot$patches[["layout"]][["design"]][["r"]])
    nrow <- max(plot$patches[["layout"]][["design"]][["b"]])
    pwalk(list(
      as.list(plot),
      seq_along(plot),
      plot$patches[["layout"]][["design"]][["t"]],
      plot$patches[["layout"]][["design"]][["l"]],
      plot$patches[["layout"]][["design"]][["b"]],
      plot$patches[["layout"]][["design"]][["r"]]
    ), function(x, i, t, l, b, r) {
      dff(
        sprintf("%s/e%s", path, i),
        x,
        w * (((r - l) + 1) / ncol),
        h * (((b - t) + 1) / nrow),
        if (is.logical(rec)) {
          rec
        } else {
          rec - 1
        }
      )
    })
  }
}

# if script is run non-interactively, save all figures
if (!interactive()) {
  dff("out/fig/main", main, 8.5 * 600, 11 * 600, 1)
  dff("out/fig/supp", supp, 8.5 * 600, 11 * 600, 1)
}
