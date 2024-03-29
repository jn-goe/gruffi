% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gruffi.R
\name{IntersectWithExpressed}
\alias{IntersectWithExpressed}
\title{IntersectWithExpressed}
\usage{
IntersectWithExpressed(genes, obj = combined.obj, genes.shown = 10)
}
\arguments{
\item{genes}{Charcter vector of genes}

\item{obj}{Seurat single cell object, Default: combined.obj}

\item{genes.shown}{Number of genes shown, Default: 10}
}
\description{
Intersect a gene set with the list of expressed genes.
}
\seealso{
\code{\link[Stringendo]{iprint}}
}
