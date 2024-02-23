##########################################################################################
# Gruffi - extended functionality for stress removal
##########################################################################################




# _____________________________________________________________________________________________ ----
# 1. Granules ---------------------------------------------------------------------------


# _____________________________________________________________________________________________ ----
# 5. Visualization ---------------------------------------------------------------------------

# _________________________________________________________________________________________________
#' @title UMAP3Dcubes
#'
#' @description UMAP with 3D cubes.
#' @param obj Seurat single cell object, Default: combined.obj
#' @param n_bin Number of bins, Default: 10
#' @param plot Plot results? Default: FALSE
#' @param save.plot Save plot into a file. Default: FALSE
#' @param reduction Which dim. reduction to use? Default: c("umap", "pca", "tsne")[1]
#' @param ident Identity (metadata), Default: Seurat.utils::GetClusteringRuns(obj)[1]
#' @param PCA_dim PCA dimensions, Default: 50
#' @seealso
#'  \code{\link[Seurat]{RunUMAP}}, \code{\link[Seurat]{RunTSNE}}, \code{\link[Seurat]{reexports}}
#'  \code{\link[sm]{binning}}
#'  \code{\link[raster]{modal}}
#'  \code{\link[viridis]{reexports}}
#'  \code{\link[rgl]{open3d}}, \code{\link[rgl]{shade3d}}, \code{\link[rgl]{matrices}}, \code{\link[rgl]{cube3d}}, \code{\link[rgl]{rglwidget}}
#'  \code{\link[htmlwidgets]{saveWidget}}
#' @export
#' @importFrom Seurat RunUMAP RunTSNE Cells
#' @importFrom sm binning
#' @importFrom raster modal
#' @importFrom viridis cividis
#' @importFrom rgl open3d shade3d translate3d cube3d rglwidget
#' @importFrom htmlwidgets saveWidget

UMAP3Dcubes <- function(obj = combined.obj,
                        n_bin = 10, plot = FALSE, save.plot = FALSE,
                        reduction = c("umap", "pca", "tsne")[1],
                        ident = Seurat.utils::GetClusteringRuns(obj)[1],
                        PCA_dim = 50) {

  warning("This function is not part of the standard workflow, and thus not fully supported.\n
          It requires the following additional packages: raster, viridis, rgl, htmlwidgets, sm")

  if (dim(obj[[reduction]]@cell.embeddings)[2] < 3) {
    if (reduction == "umap") {
      obj.3d <- Seurat::RunUMAP(obj, dims = 1:PCA_dim, n.components = 3)
      obj@misc$"reductions.backup"$umap3d <- obj.3d@reductions$umap
    }
    if (reduction == "tsne") {
      obj.3d <- Seurat::RunTSNE(obj, dims = 1:PCA_dim, n.components = 3L)
      obj@misc$"reductions.backup"$tsne3d <- obj.3d@reductions$tsne
    }
  }

  reduction <- paste0(reduction, "3d")

  umap_1 <- obj@misc$"reductions.backup"[[reduction]]@cell.embeddings[, 1]
  umap_2 <- obj@misc$"reductions.backup"[[reduction]]@cell.embeddings[, 2]
  umap_3 <- obj@misc$"reductions.backup"[[reduction]]@cell.embeddings[, 3]

  xzy <- sm::binning(cbind(umap_1, umap_2, umap_3), nbins = n_bin)
  xzy$cube_ID <- 1:dim(xzy$x)[1]

  x_width <- abs(sort(unique(xzy$x[, 1]))[2] - sort(unique(xzy$x[, 1]))[1])
  y_width <- abs(sort(unique(xzy$x[, 2]))[2] - sort(unique(xzy$x[, 2]))[1])
  z_width <- abs(sort(unique(xzy$x[, 3]))[2] - sort(unique(xzy$x[, 3]))[1])

  maj_l <- logical(length = dim(xzy$x)[1])
  cube_ID <- rep(NA, length(Seurat::Cells(obj)))

  for (l in 1:dim(xzy$x)[1]) {
    cells_in_cube <- which(umap_1 <= (xzy$x[l, 1] + x_width / 2) &
                             umap_2 <= (xzy$x[l, 2] + y_width / 2) &
                             umap_3 <= (xzy$x[l, 3] + z_width / 2) &
                             umap_1 > (xzy$x[l, 1] - x_width / 2) &
                             umap_2 > (xzy$x[l, 2] - y_width / 2) &
                             umap_3 > (xzy$x[l, 3] - z_width / 2))

    cube_ID[cells_in_cube] <- xzy$cube_ID[l]

    val <- raster::modal(obj@meta.data[ident][cells_in_cube, ])
    if (length(val) > 1) {
      val <- sample(val, 1)
    }
    maj_l[l] <- val
  }

  cols <- viridis::cividis(length(levels(factor(maj_l))))[factor(maj_l)]

  positions <- data.frame(xzy$x, freq = xzy$x.freq / max(xzy$x.freq))

  positions$umap_1 <- positions$umap_1 / x_width * 2
  positions$umap_2 <- positions$umap_2 / y_width * 2
  positions$umap_3 <- positions$umap_3 / z_width * 2

  if (plot) {
    rgl::open3d()
    for (i in 1:dim(positions)[1]) {
      rgl::shade3d(
        rgl::translate3d(
          rgl::cube3d(col = cols[i]),
          positions$umap_1[i],
          positions$umap_2[i],
          positions$umap_3[i]
        ),
        alpha = positions[i, 4]
      )
    }
    if (save.plot) {
      htmlwidgets::saveWidget(rgl::rglwidget(), "3d_UMAP_cube.html")
    }
  }
  cube_ID <- as.factor(cube_ID)
  obj@meta.data[paste0("cube_ID_", n_bin)] <- cube_ID
  return(obj)
}

