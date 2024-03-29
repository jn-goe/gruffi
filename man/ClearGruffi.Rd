% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gruffi.R
\name{ClearGruffi}
\alias{ClearGruffi}
\title{Remove gruffi results from the Seurat object}
\usage{
ClearGruffi(
  obj,
  pattern.GO.score = "^Score\\\\.GO\\\\.[0-9]+",
  pattern.granule.score = ".*_snn_res\\\\..*\\\\_cl\\\\.av_GO\\\\.[0-9]+",
  stress.assignment = "is.Stressed",
  remove.granules = FALSE
)
}
\arguments{
\item{obj}{A Seurat object.}

\item{pattern.GO.score}{A regex pattern to match GO score columns that should be removed.
Default is "^Score\\.GO\\.\link{0-9}+", which matches columns starting with "Score.GO."
followed by any number of digits.}

\item{pattern.granule.score}{A regex pattern to match granule score columns that should be
removed. Default is ".\emph{_snn_res\\..}\\_cl\\.av_GO\\.\link{0-9}+", matching a specific pattern
related to granule scores.}

\item{stress.assignment}{The name of the stress assignment column. Default: "is.Stressed".}

\item{remove.granules}{Option to also remove the optimal granule resolution. Default: FALSE,
because it can take very long time to recompute.}
}
\value{
The modified Seurat object with specified columns removed from @meta.data
and the @misc$gruffi slot cleared.
}
\description{
Removes corresponding metadata columns and clear @misc$gruffi slot
}
\examples{
obj <- removeColumnsAndClearGruffi(obj)
}
