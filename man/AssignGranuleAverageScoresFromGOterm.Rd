% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gruffi.R
\name{AssignGranuleAverageScoresFromGOterm}
\alias{AssignGranuleAverageScoresFromGOterm}
\title{Assign granule average scores from Gene Ontology (GO) terms.}
\usage{
AssignGranuleAverageScoresFromGOterm(
  obj = combined.obj,
  GO_term = "GO:0034976",
  new_GO_term_computation = FALSE,
  clustering = GetGruffiClusteringName(obj),
  save.UMAP = FALSE,
  plot.each.gene = FALSE,
  assay = "RNA",
  description = (if (!is.null(names(GO_term))) names(GO_term) else NULL),
  mirror = NULL,
  stat.av = c("mean", "median", "normalized.mean", "normalized.median")[3],
  overwrite.misc.GO_genes = FALSE,
  ...
)
}
\arguments{
\item{obj}{Seurat single cell object, Default: combined.obj}

\item{GO_term}{The GO term identifier for evaluation, default: 'GO:0034976'.}

\item{new_GO_term_computation}{Calculate new GO-term score? (or use am existing one?) Default: FALSE}

\item{clustering}{Which clustering to use (from metadata)? Default: GetGruffiClusteringName(obj)
e.g. "integrated_snn_res.48.reassigned"}

\item{save.UMAP}{Logical indicating whether to save UMAP plots, default: \code{FALSE}.}

\item{plot.each.gene}{Logical indicating whether to plot each gene's expression on UMAP, default: \code{FALSE}.}

\item{assay}{The assay type to use, default: 'RNA'.}

\item{description}{Description to added to plot title, e.g. GO-terms name, Default: 'NULL'}

\item{mirror}{Which Ensembl mirror to use in biomaRt::useEnsembl()? If connection to Ensembl
fails with default settings, mirror can be specified. Default: 'NULL'}

\item{stat.av}{How to caluclate the central tendency? Default: c("mean", "median", "normalized.mean", "normalized.median")\link{3}}

\item{overwrite.misc.GO_genes}{Overwrite preexisting gene set in \code{misc$gruffi$GO[[make.names(GO)]]}.}

\item{...}{Additional parameters to be passed to the CalculateAndPlotGoTermScores function.}
}
\description{
Calculates and assigns granule average scores from Gene Ontology (GO) term in a
Seurat object. It can compute new GO term scores, plot gene expressions, and save UMAP visualizations.
The function updates the Seurat object with cluster average scores for the specified GO term.
}
