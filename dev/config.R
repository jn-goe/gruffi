# Configuration for the Package
DESCRIPTION <- list(
  package.name = "gruffi",
  version = "1.5.5",
  title = "Gruffi identifies and removes stressed cells from brain organoid single-cell datasets.",
  description = "The Gruffi R package helps you (1) to identify stressed cells in single-cell
    RNA-seq datasets using  *granular funcitonal filtering*, and (2) you can use it to calculate any
    GO-term defined gene set's pathway activity. Gruffi integrates into single-cell analysis with
    Seurat and comes with a graphical user interface.",

  # THIS PASRT NEEDS TO BE CORRECTED MANUALLY TO INCLUDE JULIA BC PACKAGE TOOLS BREAKS WITH 2 AUTHORS
  author.given = "Abel",
  author.family = "Vertesy",
  author.email = "abel.vertesy@imba.oeaw.ac.at",
  github.user = "vertesy",

  license = "GPL-3 + file LICENSE",
  depends = "Seurat, magrittr, Stringendo (>= 0.5.0), MarkdownReports, DOSE",
  imports = "cowplot, dplyr, ggplot2, raster, rgl
                    , MarkdownHelpers (>= 1.0.1), CodeAndRoll2, Seurat.utils, ggExpress, stringr
                    , IRanges, Matrix, biomaRt, clipr, rlang, shiny, tictoc, viridis, org.Hs.eg.db, AnnotationDbi",
  suggests = "htmlwidgets, sm",
  bug.reports = "https://github.com/jn-goe/gruffi/issues/"
)
