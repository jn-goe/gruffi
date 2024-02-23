# Configuration for the Package
DESCRIPTION <- list(
  package.name = "gruffi"
  , version = "1.3.0"
  , title = "Gruffi identifies and removes stressed cells from brain organoid single-cell datasets."
  , "Author" = c(
    person(given = "Julia", family = "Naas", email = "julia.naas@meduniwien.ac.at", role =  c("aut", "cre") ),
    person(given = "Abel", family = "Vertesy", email = "abel.vertesy@imba.oeaw.ac.at", role =  c("aut", "cre") )
  )
  , "Description" = "The Gruffi R package helps you (1) to identify stressed cells in single-cell
    RNA-seq datasets using  *granular funcitonal filtering*, and (2) you can use it to calculate any
    GO-term defined gene set's pathway activity. Gruffi integrates into single-cell analysis with
    Seurat and comes with a graphical user interface."
  , "License" = "GPL-3 + file LICENSE"
  , "Packaged" =  Sys.time()
  , "Depends" =  "Seurat, magrittr, Stringendo (>= 0.5.0), MarkdownReports"
  , "Imports" = "cowplot, dplyr, ggplot2, DOSE
                    , MarkdownHelpers (>= 1.0.1), CodeAndRoll2, Seurat.utils, ggExpress, stringr, AnnotationDbi
                    , IRanges, Matrix, biomaRt, clipr, org.Hs.eg.db, rlang, shiny, tictoc, viridis"
  , "Suggests" = "raster, viridis, rgl, htmlwidgets, sm"
  , "BugReports"= "https://github.com/jn-goe/gruffi/issues/"
)


