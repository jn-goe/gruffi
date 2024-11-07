#' @title Wrapper for clUMAP with Thresholding and Plotting
#'
#' @description This function applies thresholding to granule scores within a specified object and column,
#' inverts the thresholding if a specific miscellaneous name is provided, and then calls clUMAP
#' for plotting based on the thresholded values.
#'
#' @param obj An object containing the data and miscellaneous threshold information.
#' Default is `combined.obj`.
#' @param colname Name of the column containing granule scores.
#' Default is `'RNA_snn_res.6.reassigned_cl.av_GO:0042063'`.
#' @param miscname Name of the miscellaneous threshold information within `obj@misc$gruffi`.
#' Default is `'thresh.stress.ident1'`.
#' @param auto Automatically flip filtering
#' @param ... Additional arguments to be passed to `clUMAP`.
#'
#' @importFrom Seurat.utils clUMAP
#' @importFrom CodeAndRoll2 as.numeric.wNames.character
#'
#' @examples
#' # Assuming combined.obj is available and properly formatted
#' customClUMAPWrapper(
#'   combined.obj,
#'   "RNA_snn_res.6.reassigned_cl.av_GO:0042063",
#'   "thresh.stress.ident1"
#' )
#'
#' @export
StressUMAP <- function(obj = combined.obj,
                        colname = "RNA_snn_res.6.reassigned_cl.av_GO:0042063",
                        miscname = "thresh.stress.ident1",
                        auto = TRUE,
                        ...) {

  if (is.null(colname)) {
    message("Score is NULL")
    return(NULL)
  }

  stopifnot(
    inherits(obj, "Seurat"), !is.null(obj@misc$gruffi),
    miscname %in% names(obj@misc$gruffi),
    colname %in% names(obj@meta.data)
  )

  thr <- obj@misc$gruffi[[miscname]] # Retrieve threshold from object

  # Convert granule scores to numeric with names
  granule_scores <- CodeAndRoll2::as.numeric.wNames.character(obj@meta.data[[colname]], verbose = FALSE)

  # Apply threshold to create a logical vector for identification
  obj[["pass"]] <- (granule_scores < thr)

  # Invert the logical vector for a specific miscname
  if (miscname == "thresh.notstress.ident3" & auto) {
    message("thresh.notstress.ident3 <> filtering is flipped")
    obj@meta.data[["pass"]] <- !(obj@meta.data[["pass"]])
  }

  components <- strsplit(colname, "_cl\\.av_")[[1]]

  nr.clust <- NaN
  try(nr.clust <- length(unique(obj@meta.data[[components[1]]])), silent = TRUE)

  subt <- paste(
    "threshold:", thr, " | Pass: ", pc_TRUE(obj@meta.data[["pass"]]), " | ",
    ncol(obj), "cells |", nr.clust, "clusters."
  )
  # Call clUMAP with the thresholded identification and any additional arguments
  Seurat.utils::clUMAP(
    ident = "pass", obj = obj, label = FALSE,
    title = paste("Granule Thresholding", make.names(components[2])),
    prefix = colname, sub = subt, caption = colname, ...
  )
}


