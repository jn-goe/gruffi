% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gruffi.R
\name{StressUMAP}
\alias{StressUMAP}
\title{Plot Stress Assignment UMAP}
\usage{
StressUMAP(
  obj = combined.obj,
  colname = "is.Stressed",
  GO.terms = names(obj@misc$gruffi$GO),
  suffix = NULL,
  ...
)
}
\arguments{
\item{obj}{An object containing the data and miscellaneous threshold information.
Default is \code{combined.obj}.}

\item{colname}{Name of the column containing granule scores.
Default is \code{'RNA_snn_res.6.reassigned_cl.av_GO:0042063'}.}

\item{GO.terms}{used for the computation}

\item{suffix}{Suffix to append to the plot/file name. Default: NULL.}

\item{...}{Additional arguments to be passed to \code{clUMAP}.}
}
\description{
Plot Stress Assignment UMAP after running the gruffi pipeline
}
\examples{
StressUMAP(combined.obj)
}
