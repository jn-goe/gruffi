##########################################################################################
# Gruffi - finding stressed cells.
##########################################################################################




# _____________________________________________________________________________________________ ----
# 1. Granules ---------------------------------------------------------------------------


# _________________________________________________________________________________________________
#' @title AutoFindGranuleResolution
#'
#' @details The function iterates over clustering resolutions within the specified bounds to find
#' an optimal resolution that results in median cluster sizes within the target range. It uses
#' Seurat's clustering functionality and adjusts based on the median cell count per cluster.
#'
#' @param obj Seurat object to perform clustering on. Defaults to 'combined.obj'.
#' @param min.med.granule.size Minimum acceptable median cell count per cluster. Default is 100.
#' @param max.med.granule.size Maximum acceptable median cell count per cluster. Default is 200.
#' @param min.res Starting lower bound for clustering resolution search.
#' Defaults to: `round(ncol(obj)/10000)`.
#' @param max.res Starting upper bound for clustering resolution search.
#' Defaults to: `round(ncol(obj)/1000)`.
#' @param assay Specifies the assay to use for clustering. Default is 'integrated, with
#' "RNA" as an alternative taken from on `Seurat::DefaultAssay(obj)`.
#' @param max.loop Maximum number of iterations for adjusting clustering resolution. Default is 20.
# #' @param n.threads Number of threads to use for parallel computation. Default is 1.
#'
#' @return A Seurat object with updated clustering at the optimal resolution found, including
#' modifications to `obj@meta.data` and `Idents(obj)` to reflect the new clustering.
#' @examples
#' # Assuming `combined.obj` is a pre-loaded Seurat object
#' optimal.obj <- AutoFindGranuleResolution(obj = combined.obj)
#'
#' @seealso
#'  \code{\link[Seurat]{reexports}}, \code{\link[Seurat]{FindClusters}}
#'  \code{\link[Stringendo]{iprint}}
#' @export
#' @importFrom Stringendo iprint
#' @importFrom tictoc tic toc
#' @importFrom future plan nbrOfWorkers
AutoFindGranuleResolution <- function(obj = combined.obj,
                                      min.med.granule.size = 100,
                                      max.med.granule.size = 200,
                                      min.res = round(ncol(obj) / 1e4),
                                      max.res = round(ncol(obj) / 1e3),
                                      assay = c("integrated", "RNA")[1], # Seurat::DefaultAssay(obj),
                                      max.loop = 20,
                                      n.threads = 1) {
  message(
    ncol(combined.obj), " cells in object. Searching for optimal resolution between res: ",
    min.res, " & ", max.res, "...\nGranule Size Requirements: ",
    min.med.granule.size, " >X> ", max.med.granule.size
  )

  assays.present <- names(obj@assays)
  if (assay %in% assays.present) {
    Seurat::DefaultAssay(obj) <- assay
  } else {
    message("Assay: (", assay, ") is not found in the object. Replaced by: (", assays.present[1], ")")
    assay <- assays.present[1]
  }

  r.current <- min.res
  r.lower <- min.res
  r.upper <- max.res

  if (n.threads > 1) {
    warning("Using multicore for parallel computation is not debgugged.", immediate. = TRUE)
    z <- require(future)
    stopifnot("future package needed for multicore" = z)
    future::plan("multisession", workers = n.threads)
    message("Threads: ", future::nbrOfWorkers())
  }

  message("\nClustering at res.: ", r.current)
  tictoc::tic(); obj <- Seurat::FindClusters(obj, resolution = r.current, verbose = FALSE); tictoc::toc()
  m <- CalculateMedianClusterSize(obj, assay, res = r.current)
  stopifnot("Define a lower `min.res` or decrease `min.med.granule.size` if necessary. Granule size too small" = m > min.med.granule.size)


  message("\nClustering at res.: ", r.upper)
  tictoc::tic(); obj <- Seurat::FindClusters(obj, resolution = r.upper, verbose = FALSE); tictoc::toc()
  m.up <- CalculateMedianClusterSize(obj, assay, res = r.upper)
  stopifnot("Define a higher `max.res`. Granule size too big" = m.up < max.med.granule.size)



  loop.count <- 1
  while (m > max.med.granule.size | m < min.med.granule.size) {
    # message("Current clustering resolution: ", r.current)
    # message("Current median cluster size: ", m)

    obj@meta.data[[paste0(assay, "_snn_res.", r.current)]] <- NULL

    if (m > max.med.granule.size) {
      r.lower <- r.current
      r.current <- round((max.res - r.current) / 2 + r.current)
    }

    if (m < min.med.granule.size) {
      r.upper <- r.current
      r.current <- round((r.current - min.res) / 2 + min.res)
    }

    message("\nSearch between resolutions: ", r.lower, " & ", r.upper, ".")
    tictoc::tic(); obj <- Seurat::FindClusters(obj, resolution = r.current, verbose = FALSE); tictoc::toc()

    res.name <- paste0(assay, "_snn_res.", r.current)
    m <- CalculateMedianClusterSize(obj, assay, res = r.current)

    loop.count <- loop.count + 1
    if (loop.count > max.loop) {
      next
    }
  }
  res.name <- paste0(assay, "_snn_res.", r.current)
  clustering <- obj@meta.data[[res.name]]

  obj@misc$gruffi$"optimal.granule.res" <- res.name
  message("Suggested Granule resolution (meta.data colname) is stored in: obj@misc$gruffi$optimal.granule.res:\n",
          res.name)

  Seurat::Idents(obj) <- clustering

  nr.granules <- length(unique(clustering))
  if (nr.granules < 100) {
    warning("Too few granules is too small. Try to define a higher `min.res`.", immediate. = TRUE)
    message(" -> nr.granules: ", nr.granules)
  }

  cv.final <- signif(cv(table(clustering)))
  if (cv.final > 0.5) {
    warning("Granule sizes vary significantly. Try to define a higher `min.res`.", immediate. = TRUE)
    message(" -> CV: ", Stringendo::percentage_formatter(cv.final))
  }

  print(paste0(
    "\nSuggested resolution is ", r.current,
    " and the respective clustering is stored in Idents(obj) and in @meta.data.",
    " The median numbers of cells per cluster is ", m,
    ". The number of clusters is ", length(unique(clustering))
  ))
  return(obj)
}


# _________________________________________________________________________________________________
#' @title ReassignSmallClusters
#'
#' @description Reassign granules (clusters) smaller than X to nearby granules based on 3D UMAP coordinates.
#' @param obj Seurat single cell object, Default: combined.obj
#' @param ident Identity (granules / clustering) to reassign, Default: obj@misc$gruffi$optimal.granule.res
#' @param cell.num Minimum number of cells per granule / cluster, Default: 30
#' @param reduction Dim. reduction used to estimate closest cluster center when reassigning. Default: c("pca", "3d_umap")[2]
#' @param nr.of.PCs.for.3D.umap.computation If 3D umap need to be recomputed, how many PCA-dimensions to be used?, Default: ncol(obj@reductions$pca@cell.embeddings)
#' @seealso
#'  \code{\link[Stringendo]{iprint}}
#'  \code{\link[Seurat]{RunUMAP}}, \code{\link[Seurat]{reexports}}
#' @export
#' @importFrom Stringendo iprint
#' @importFrom Seurat RunUMAP Idents RenameIdents

ReassignSmallClusters <- function(obj = combined.obj,
                                  ident = obj@misc$gruffi$"optimal.granule.res",
                                  cell.num = 30,
                                  reduction = c("pca", "3d_umap")[2],
                                  nr.of.PCs.for.3D.umap.computation = ncol(obj@reductions$pca@cell.embeddings)) {

  # Create a new identifier for reassigned granules
  name.id.reassigned <- paste0(ident, ".reassigned")
  message("Reassigned granule res. column: ", name.id.reassigned)

  granules <- obj@meta.data[[name.id.reassigned]] <- obj@meta.data[[ident]]
  granules.sizes <- table(granules) # Calculate sizes of all granules
  granules.too.small.bool <- (granules.sizes < cell.num) # Identify granules smaller than the specified minimum cell number

  # If there are small granules, proceed with reassignment
  if (sum(granules.too.small.bool) > 0) {

    # Choose the dimensionality reduction technique for calculating embeddings
    if (reduction == "pca") {
      message("PCA cell.embeddings used.")
      embedding <- obj@reductions$pca@cell.embeddings

    } else if (reduction == "3d_umap") {
      message("3D umap cell.embeddings used.")

      if (is.null(obj@misc$"reductions.backup")) {
        message("The 3D UMAP will now be computed. PC dimensions used: ", nr.of.PCs.for.3D.umap.computation, "\n")

        obj.3d <- Seurat::RunUMAP(obj, dims = 1:nr.of.PCs.for.3D.umap.computation, n.components = 3)

        obj@misc$reductions.backup$"umap3d" <- obj.3d@reductions$umap
        message("3D UMAP embedding is now saved to: obj@misc$reductions.backup")

      } else {
        message("A 3D UMAP in @misc$'reductions.backup is found and will be used.\n")
      } # if exists / calculate

      embedding <- obj@misc$"reductions.backup"$umap3d@cell.embeddings
      stopifnot(ncol(embedding) == 3)
    }  # if reduction == "3d_umap"
  } # if any granules too small


  # Iterate over reassignment process until no small granules remain
  while (sum(table(obj@meta.data[[name.id.reassigned]]) < cell.num) > 0) {

    cl.reassign <- names(which(table(obj@meta.data[[name.id.reassigned]]) < cell.num))
    message("\n", length(cl.reassign), " clusters have <",
            cell.num, " cells, and need to be reassigned.")

    # Calculate the center (mean position) of each cluster in the embedding space
    clustering <- obj@meta.data[[name.id.reassigned]]
    cluster.means <- matrix(NA, nrow = length(unique(clustering)), ncol = dim(embedding)[2])
    rownames(cluster.means) <- levels(clustering)
    for (d in 1:dim(embedding)[2]) {
      cluster.means[, d] <- by(embedding[, d], clustering, mean)
    }

    # Prepare for reassignment by calculating distances between cluster centers
    dist.mat <- as.matrix(stats::dist(cluster.means))
    diag(dist.mat) <- Inf

    # Identify and reassign small clusters based on proximity to larger ones
    names(new_names) <- new_names <- levels(clustering)
    for (r in cl.reassign) {
      Stringendo::iprint("cluster:", r)
      new_names[which(new_names == r)] <- names(which.min(dist.mat[r, ]))
    }

    # Update identities in the Seurat object with reassigned cluster names
    Seurat::Idents(obj) <- obj@meta.data[[name.id.reassigned]]
    obj <- Seurat::RenameIdents(obj, new_names)
    obj@meta.data[[name.id.reassigned]] <- Seurat::Idents(obj)
  } # while

  message("\nReassignment resulted in " , length(levels(obj@meta.data[[name.id.reassigned]])),
          " granules down from ", length(levels(obj@meta.data[[ident]])), " original granules.")

  # Update the optimal granule resolution in the Seurat object's @misc slot
  obj@misc$gruffi$"optimal.granule.res" <- name.id.reassigned
  message("Suggested Granule Resolution (meta.data colname) is stored in:")
  message("obj@misc$gruffi$optimal.granule.res:\n", name.id.reassigned)
  message("\nCall PlotClustSizeDistr() to check results.\n")

  CalculateMedianClusterSize(columnName = ident, obj = obj)
  CalculateMedianClusterSize(columnName = name.id.reassigned, obj = obj)
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Calculate Median Cluster Size
#'
#' @description This function calculates and returns the median size of clusters for a given assay
#' and resolution in a Seurat object. It provides an option to display the median size via message.
#'
#' @param obj A Seurat object containing metadata with clustering information.
#' @param assay The name of the assay to which the resolution parameter is applied.
#' @param res The resolution level for which to calculate the median cluster size.
#' @param suffix String suffix for clustering name
#' @param columnName Directly provide a meta.data column name instead of parsing it.
#' @param q quantile. Default: 0.5, i.e: the median.
#' @param verbose Logical; if TRUE, prints the median cluster size to the console. Defaults to TRUE.
#'
#' @return The median size of the clusters at the specified resolution.
#'
#' @details The function generates a column name based on the assay and resolution parameters,
#' checks for the existence of this column in the object's metadata, calculates the median size
#' of the clusters, and optionally outputs this information as a message.
#'
#' @examples
#' # Assuming `seuratObj` is a Seurat object with appropriate metadata
#' medianSize <- CalculateMedianClusterSize(obj = seuratObj, assay = "RNA", res = 0.8)
#'
#' @export
CalculateMedianClusterSize <- function(
    obj,
    assay = NULL, res = NULL, suffix = "_snn_res.",
    columnName = NULL,
    q = 0.5,
    verbose = TRUE) {
  # Generate the resolution-specific metadata column name

  if(is.null(columnName)) columnName <- paste0(assay, suffix, res)

  # Ensure the column exists in the metadata
  stopifnot(columnName %in% colnames(obj@meta.data))

  # Calculate the median size of clusters at the specified resolution
  cluster.sizes <- table(obj@meta.data[[columnName]])
  median.size <- quantile(cluster.sizes, probs = q)
  cv.size <- signif(cv(cluster.sizes))

  # Output the median cluster size
  if (verbose) {
    prefix <- if (q == 0.5) "Median" else paste("quantile", q)
    if (is.null(res)) res <- columnName
    message("Resolution: ", res, " | ", length(cluster.sizes), " clusters | ",
            prefix, " size: ", median.size,
            " | CV of sizes: ", Stringendo::percentage_formatter(cv.size))
  }

  return(median.size)
}


# _____________________________________________________________________________________________ ----
# 2. Obtain and prepare GO Terms ---------------------------------------------------------------------------

#' @title Retrieve Gene Ontology (GO) Terms and Associated Genes
#'
#' @description Fetches genes associated with a specified Gene Ontology (GO) term and optionally
#' opens the GO term page. The genes are retrieved using either the Ensembl database via biomaRt or
#'  from precomputed GO_genes results stored in the Seurat object.
#'
#' @param obj A Seurat object potentially containing precomputed GO_genes results. Default: combined.obj
#' @param GO The GO term identifier for which to fetch associated genes. Default: 'GO:0034976'
#' @param assay Assay to consider when retrieving precomputed GO_genes results.
#' @param data.base.access Which tool to use to access gene list databases? Options:
#' biomaRt (deault) or AnnotationDbi.
#' @param version Ensembl version, useful if wanting to connect to an archived Ensembl version.
#' @param mirror Ensembl mirror to use, helpful if the default connection fails with
#' default settings, mirror can be specified. Default: 'NULL'
#' @param GRCh GRCh version if not using the current GRCh38, currently this can only be 37 (via useEnsembl())
#' @param web.open Logical; if TRUE, opens the web page for the specified GO term.
#' @param genes.shown Number of genes shown, Default: 10
#'
#' @return The input Seurat object with the retrieved genes associated with the
#' specified GO term stored in `obj@misc$gruffi$GO`.
#'
#' @seealso
#'  \code{\link[biomaRt]{useEnsembl}}, \code{\link[biomaRt]{getBM}}
#'  \code{\link[AnnotationDbi]{AnnotationDb-objects}}
#'  \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}
#'  \code{\link[Stringendo]{kollapse}}
#' @export
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom AnnotationDbi select
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom Stringendo kollapse
#' @importFrom DOSE geneInCategory

GetGOTerms <- function(obj = combined.obj,
                       GO = "GO:0034976",
                       assay = "RNA",
                       data.base.access = c("biomaRt", "AnnotationDbi")[1],
                       version = NULL,
                       GRCh = NULL,
                       genes.shown = 10,
                       web.open = FALSE,
                       overwrite.misc.GO_genes = FALSE,
                       mirror = NULL) {

  message("Running GetGOTerms()")

  GO.gene.symbol.slot <- obj@misc$gruffi$GO[[make.names(GO)]]
  # GO.gene.symbol.slot <- obj@misc$gruffi$"GO_genes"[[assay]]

  if(overwrite.misc.GO_genes) GO.gene.symbol.slot <- NULL # so it will be recomputed and overwritten

  if (is.null(GO.gene.symbol.slot)) {
    message("No pre-existing GO_genes results found in obj@misc$gruffi$GO[[make.names(GO)]].")

    if (data.base.access == "biomaRt") { # If using biomaRt to access Ensembl database

      # Check if Ensembl mart object already exists, if not, create one
      if (!exists("ensembl")) { message("Running biomaRt::useEnsembl()")
        ensembl <<- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = version, GRCh = GRCh, mirror = mirror)
      }

      # Retrieve gene symbols associated with the GO term using biomaRt
      genes <- biomaRt::getBM(
        attributes = c("hgnc_symbol"), # 'ensembl_transcript_id', 'go_id'
        filters = "go_parent_term", uniqueRows = TRUE,
        values = GO, mart = ensembl
      )[, 1]
      message("Gene symbols downloaded from ensembl via biomaRt::getBM()")

    } else if (data.base.access == "AnnotationDbi") { message("Running AnnotationDbi::select()")
      genes <- unique(AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db, keys = GO, columns = c("SYMBOL"), keytype = "GOALL")$SYMBOL)
      message("Gene symbols downloaded from AnnotationDbi / org.Hs.eg.db")

    } else {
      stop("Choose either AnnotationDbi or biomaRt for gene list retrieval.")
    }

  } else {
    message("Pre-existing GO_genes results found in obj@misc$gruffi$GO[[make.names(GO)]].")
    genes <- unlist(DOSE::geneInCategory(GO.gene.symbol.slot)[GO])
    message("Gene symbols taken from GO_genes.")
  }
  message("\n", length(genes), " gene symbols used: ", kollapse(head(genes, n = genes.shown)))

  # Intersect the downloaded genes with those expressed in the Seurat object
  genes <- IntersectWithExpressed(obj = obj, genes = genes)

  # Save the result in the Seurat object for later access
  if (is.null(obj@misc$gruffi$"GO")) obj@misc$gruffi$GO <- list()
  obj@misc$gruffi$GO[[make.names(GO)]] <- genes
  message("Genes in ", GO, " are saved under obj@misc$gruffi$GO$", make.names(GO))

  # Open the GO term web page if requ
  if (web.open) system(paste0("open https://www.ebi.ac.uk/QuickGO/search/", GO))
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Show All GO Terms in an Object, and retrieve corresponding names using AnnotationDbi
#'
#' @description Retrieves and optionally annotates all GO term identifiers present in the
#' metadata of a Seurat object. It identifies columns starting with 'Score.GO',
#' extracts the corresponding GO terms, retrieves tge corresponding term-names using AnnotationDbi,
#' and optionally saving them back into the object.
#' @param obj A Seurat object containing metadata with GO term score columns, default: `combined.obj`.
#' @param return.obj Logical indicating whether to return the modified Seurat object.
#'
#' @return A modified Seurat object with GO terms annotated in `misc$gruffi$GO.Lookup` or
#' a vector of GO term names.
#' @seealso
#'  \code{\link[CodeAndRoll2]{grepv}}
#'  \code{\link[AnnotationDbi]{GOID}}, \code{\link[AnnotationDbi]{GOTerms-class}}
#' @export
#' @importFrom CodeAndRoll2 grepv
#' @importFrom AnnotationDbi Term
#' @importFrom gtools mixedsort
#'
#' @export
GetAllGOTermNames <- function(obj = combined.obj, return.obj = TRUE) {
  # Foremerly GetAllGOTerms

  x <- obj@meta.data
  GOz <- CodeAndRoll2::grepv(pattern = "^Score.GO", x = colnames(x), perl = TRUE)
  GOx <- gtools::mixedsort(unique(CodeAndRoll2::grepv(
    pattern = "\\.[0-9]$", x = GOz, perl = TRUE,
    invert = TRUE
  )))
  GO.names <- AnnotationDbi::Term(object = .convert.score.2.GO_term(GOx))

  if (return.obj) {
    obj@misc$gruffi$"GO.Lookup" <- GO.names
    print("GO IDs present in @meta.data are now saved in misc$gruffi$GO.Lookup")
    cat(head(GO.names), "...")
    return(obj)
  } else {
    return(GO.names)
  }
}



# _____________________________________________________________________________________________ ----
# 3. Calculate GO Scores ---------------------------------------------------------------------------


#' @title Assign granule average scores from Gene Ontology (GO) terms.
#'
#' @description  Calculates and assigns granule average scores from Gene Ontology (GO) term in a
#' Seurat object. It can compute new GO term scores, plot gene expressions, and save UMAP visualizations.
#' The function updates the Seurat object with cluster average scores for the specified GO term.
#'
#' @param obj Seurat single cell object, Default: combined.obj
#' @param GO_term The GO term identifier for evaluation, default: 'GO:0034976'.
#' @param new_GO_term_computation Calculate new GO-term score? (or use am existing one?) Default: FALSE
#' @param clustering The clustering identity to use for evaluation, default: `GetGruffiClusteringName(obj)`.
#' @param save.UMAP Logical indicating whether to save UMAP plots, default: `FALSE`.
#' @param plot.each.gene Logical indicating whether to plot each gene's expression on UMAP, default: `FALSE`.
#' @param assay The assay type to use, default: 'RNA'.
#' @param description Description to added to plot title, e.g. GO-terms name, Default: 'NULL'
#' @param mirror Which Ensembl mirror to use in biomaRt::useEnsembl()? If connection to Ensembl
#' fails with default settings, mirror can be specified. Default: 'NULL'
#' @param stat.av How to caluclate the central tendency? Default: c("mean", "median", "normalized.mean", "normalized.median")[3]
#' @param clustering Which clustering to use (from metadata)? Default: GetGruffiClusteringName(obj)
#' e.g. "integrated_snn_res.48.reassigned"
#' @param ... Additional parameters to be passed to the PlotGoTermScores function.
#' @export
#' @importFrom Seurat Idents RenameIdents
#' @importFrom Stringendo iprint

AssignGranuleAverageScoresFromGOterm <- function(obj = combined.obj,
                                                 GO_term = "GO:0034976",
                                                 new_GO_term_computation = FALSE,
                                                 clustering = GetGruffiClusteringName(obj),
                                                 save.UMAP = FALSE,
                                                 plot.each.gene = FALSE,
                                                 assay = "RNA",
                                                 description = NULL,
                                                 mirror = NULL,
                                                 stat.av = c("mean", "median", "normalized.mean", "normalized.median")[3],
                                                 ...) {

  Seurat::Idents(obj) <- obj@meta.data[[clustering]]

  # If new GO term score computation is requested
  if (new_GO_term_computation) {
    obj <- PlotGoTermScores(
      GO = GO_term, save.UMAP = save.UMAP, obj = obj,
      desc = description, plot.each.gene = plot.each.gene, mirror = mirror, ...
    )
  }
  message("Calculating granule average scores")

  cl.av <- CalcClusterAverages_Gruffi(
    obj = obj,
    stat = stat.av,
    col_name = .convert.GO_term.2.score(GO_term),
    split_by = clustering
  )

  # Remove 'cl.' prefix from cluster average score names
  names(cl.av) <- gsub("cl.", "", names(cl.av))
  obj <- Seurat::RenameIdents(obj, cl.av)

  # Store cluster average scores in the metadata under a new column
  mScoreColName <- paste0(clustering, "_cl.av_", GO_term)
  obj@meta.data[mScoreColName] <- Seurat::Idents(obj)

  # Reset the active identities to the original clustering
  Seurat::Idents(obj) <- obj@meta.data[[clustering]]
  message("New meta.data with granule average scores: ", mScoreColName)

  return(obj)
}


# _____________________________________________________________________________________________
#' @title AddGOGeneList.manual
#'
#' @description Add a GO-term gene list under obj@misc$gruffi$GO$xxxx.
#' @param obj Seurat single cell object, Default: combined.obj
#' @param GO GO-term; Default: 'GO:0034976'
#' @param web.open Open weblink, Default: FALSE
#' @param genes Charcter vector of genes, Default: c("A0A140VKG3", "ARX", "CNTN2", "DRD1", "DRD2", "FEZF2", "LHX6")
#' @seealso
#'  \code{\link[Stringendo]{iprint}}
#' @export
#' @importFrom Stringendo iprint

AddGOGeneList.manual <- function(obj = combined.obj, GO = "GO:0034976", web.open = FALSE # Add GO terms via Biomart package.
                                 , genes = c("A0A140VKG3", "ARX", "CNTN2", "DRD1", "DRD2", "FEZF2", "LHX6")) {
  print(head(genes, n = 15))
  genes <- IntersectWithExpressed(obj = obj, genes = genes)

  if (is.null(obj@misc$gruffi$GO)) obj@misc$gruffi$GO <- list()
  obj@misc$gruffi$GO[[make.names(GO)]] <- genes
  Stringendo::iprint("Genes in", GO, "are saved under obj@misc$gruffi$GO$", make.names(GO))
  if (web.open) system(paste0("open https://www.ebi.ac.uk/QuickGO/search/", GO))
  return(obj)
}




# _________________________________________________________________________________________________
#' @title AddGOScore
#' @description Add GO Score to a Seurat object
#' @param obj Seurat single cell object, Default: combined.obj
#' @param GO GO Term. String of form "GO:xxxxxxx", Default: 'GO:0034976'.
#' @param FixName Shall trailing "1" ar the end of the new column name be removed?
#' @seealso
#'  \code{\link[Seurat]{AddModuleScore}}
#' @export
#' @importFrom Seurat AddModuleScore

AddGOScore <- function(obj = combined.obj, GO = "GO:0034976", FixName = TRUE) {
  print("AddGOScore()")
  GO.wDot <- make.names(GO)
  (genes.GO <- list(obj@misc$gruffi$GO[[GO.wDot]]))
  # print(genes.GO)
  ScoreName <- paste0("Score.", make.names(GO))
  if (!is.list(genes.GO)) genes.GO <- list(genes.GO) # idk why this structure is not consistent...
  obj <- Seurat::AddModuleScore(object = obj, features = genes.GO, name = ScoreName)

  if (FixName) obj <- .fix.metad.colname.rm.trailing.1(obj = obj, colname = ScoreName)
  return(obj)
}


# _________________________________________________________________________________________________
#' @title AddCustomScore
#'
#' @description Add a custom gene expression score from a set of genes provided to the
#' Seurat single cell object. Calculates Score for gene set. Fixes name.
#' @param obj Seurat single cell object, Default: combined.obj
#' @param genes Charcter vector of genes, Default: ""
#' @param assay.use Which assay to use for getting the genes? Default: 'RNA'.
#' @param FixName Shall trailing "1" ar the end of the new column name be removed? Default: TRUE
#' @seealso
#'  \code{\link[Seurat]{AddModuleScore}}
#' @export
#' @importFrom Seurat AddModuleScore
#' @importFrom Stringendo ppp

AddCustomScore <- function(obj = combined.obj, genes = "", assay.use = "RNA", FixName = TRUE) {
  ls.genes <- list(genes)
  if (!is.list(ls.genes)) ls.genes <- list(ls.genes) # idk why this structure is not consistent...
  (ScoreName <- Stringendo::ppp("Score", substitute(genes)))
  obj <- Seurat::AddModuleScore(object = obj, features = ls.genes, name = ScoreName, assay = assay.use)

  if (FixName) obj <- .fix.metad.colname.rm.trailing.1(obj = obj, colname = ScoreName)
  return(obj)
}




# _________________________________________________________________________________________________



# _________________________________________________________________________________________________
#' @title CustomScoreEvaluation
#'
#' @description Custom gene-set derived score evaluation for filtering.
#' @param obj Seurat single cell object, Default: combined.obj
#' @param custom.score.name Name of your custom gene set.
#' @param clustering Which clustering to use (from metadata)? Default: GetGruffiClusteringName(obj)
#' e.g. "integrated_snn_res.48.reassigned"
#' @param assay Which assay to use?, Default: 'RNA'
#' @param stat.av How to caluclate the central tendency?
#' Default: c("mean", "median", "normalized.mean", "normalized.median")[3]
#' @seealso
#'  \code{\link[Seurat]{reexports}}
#'  \code{\link[Stringendo]{iprint}}
#' @export
#' @importFrom Seurat Idents RenameIdents
#' @importFrom Stringendo iprint

CustomScoreEvaluation <- function(obj = combined.obj,
                                  custom.score.name = "Score.heGENES",
                                  clustering = GetGruffiClusteringName(obj),
                                  assay = "RNA",
                                  stat.av = c("mean", "median", "normalized.mean", "normalized.median")[3],
                                  ...) {
  Seurat::Idents(obj) <- obj@meta.data[[clustering]]
  all.genes <- rownames(obj@assays[[assay]])

  print("Calculating cl. average score")
  cl.av <- CalcClusterAverages_Gruffi(
    obj = obj,
    stat = stat.av,
    col_name = custom.score.name,
    split_by = clustering
  )
  names(cl.av) <- gsub("cl.", "", names(cl.av))
  obj <- Seurat::RenameIdents(obj, cl.av)
  mScoreColName <- paste0(clustering, "_cl.av_", custom.score.name)
  obj@meta.data[mScoreColName] <- Seurat::Idents(obj)
  Seurat::Idents(obj) <- obj@meta.data[[clustering]]
  print(mScoreColName)
  return(obj)
}


# _____________________________________________________________________________________________ ----
# 4. Filtering ---------------------------------------------------------------------------

#' @title Launch Shiny App for GO Term-based Thresholding
#'
#' @description Launches a Shiny application to interactively apply GO term-based thresholding
#' for identifying stressed cells in a Seurat object. This function calculates the
#' "granule average GO-scores" for individual GO terms and updates the object with a combined stress
#' identification based on these thresholds.#'
#' @param obj Seurat single cell object, with Granule Averages GO-Scores in metadata (ident1-4).
#' @param proposed.method The method for estimating thresholds, default: `"fitted"`.
#' Options are `"fitted"` or `"empirical"`.
#' @param quantile The quantile cutoff for threshold determination, default: `0.99`.
#' Options are `0.99` or `0.9`.
#' @param stress.ident1 Identifier for the first stress-related GO term.
#' @param stress.ident2 Identifier for the second stress-related GO term.
#' @param notstress.ident3 Identifier for the first non-stress-related GO term. Idents 3 and 4 act
#'  as a negative filter for stress identification.
#' @param notstress.ident4 Identifier for the second non-stress-related GO term, optional.
#' @param plot.cluster.shiny The clustering run to be used for plotting in the Shiny app,
#' default: `Seurat.utils::GetClusteringRuns(obj)[1]`.
#'
#' @export
#' @importFrom shiny runApp shinyApp
Shiny.GO.thresh <- function(
    obj = combined.obj,
    proposed.method = c("fitted", "empirical")[1],
    quantile = c(.99, .9)[1],
    stress.ident1,
    stress.ident2,
    notstress.ident3,
    notstress.ident4 = NULL,
    plot.cluster.shiny = Seurat.utils::GetClusteringRuns(obj)[1]) {
  app_env <- new.env()
  meta <- obj@meta.data

  # Ensure that provided GO term identifiers exist in the metadata
  stopifnot(stress.ident1 %in% colnames(meta) | is.null(stress.ident1))
  stopifnot(stress.ident2 %in% colnames(meta) | is.null(stress.ident2))
  stopifnot(notstress.ident3 %in% colnames(meta) | is.null(notstress.ident3))
  stopifnot(notstress.ident4 %in% colnames(meta) | is.null(notstress.ident4))

  # Convert categorical scores to numeric for thresholding
  gr.av.stress.scores1 <- as.numeric(levels(meta[, stress.ident1]))
  gr.av.stress.scores2 <- as.numeric(levels(meta[, stress.ident2]))
  gr.av.notstress.scores3 <- as.numeric(levels(meta[, notstress.ident3]))
  gr.av.notstress.scores4 <- as.numeric(levels(meta[, notstress.ident4]))

  # Compute threshold proposals for each "granule average GO-score" using PlotNormAndSkew function
  app_env$"thresh.stress.ident1" <- PlotNormAndSkew(gr.av.stress.scores1, q = quantile, tresholding = proposed.method, plot.hist = FALSE)
  app_env$"thresh.stress.ident2" <- PlotNormAndSkew(gr.av.stress.scores2, q = quantile, tresholding = proposed.method, plot.hist = FALSE)
  app_env$"thresh.notstress.ident3" <- PlotNormAndSkew(gr.av.notstress.scores3, q = quantile, tresholding = proposed.method, plot.hist = FALSE)
  app_env$"thresh.notstress.ident4" <- PlotNormAndSkew(gr.av.notstress.scores4, q = quantile, tresholding = proposed.method, plot.hist = FALSE)

  # Prepare slider inputs for Shiny app based on min, max, and step values calculated from scores
  # These values define the range and granularity of threshold adjustments within the Shiny app
  min.x.stress.ident1 <- floor(min(gr.av.stress.scores1, app_env$"thresh.stress.ident1"))
  max.x.stress.ident1 <- ceiling(max(gr.av.stress.scores1, app_env$"thresh.stress.ident1"))
  step.stress.ident1 <- 0.001

  min.x.stress.ident2 <- floor(min(gr.av.stress.scores2, app_env$"thresh.stress.ident2"))
  max.x.stress.ident2 <- ceiling(max(gr.av.stress.scores2, app_env$"thresh.stress.ident2"))
  step.stress.ident2 <- 0.001

  min.x.notstress.ident3 <- floor(min(gr.av.notstress.scores3, app_env$"thresh.notstress.ident3"))
  max.x.notstress.ident3 <- ceiling(max(gr.av.notstress.scores3, app_env$"thresh.notstress.ident3"))
  step.notstress.ident3 <- 0.001

  min.x.notstress.ident4 <- floor(min(gr.av.notstress.scores4, app_env$"thresh.notstress.ident4"))
  max.x.notstress.ident4 <- ceiling(max(gr.av.notstress.scores4, app_env$"thresh.notstress.ident4"))
  step.notstress.ident4 <- 0.001

  # Define the environment for the Shiny app as `app_env`.
  app_env$"idents" <- list(
    "stress.ident1" = stress.ident1,
    "stress.ident2" = stress.ident2,
    "notstress.ident3" = notstress.ident3,
    "notstress.ident4" = notstress.ident4
  )

  app_env$"sliders" <- list(
    "min.x.stress.ident1" = min.x.stress.ident1, "max.x.stress.ident1" = max.x.stress.ident1,
    "step.stress.ident1" = step.stress.ident1,
    "min.x.stress.ident2" = min.x.stress.ident2, "max.x.stress.ident2" = max.x.stress.ident2,
    "step.stress.ident2" = step.stress.ident2,
    "min.x.notstress.ident3" = min.x.notstress.ident3, "max.x.notstress.ident3" = max.x.notstress.ident3,
    "step.notstress.ident3" = step.notstress.ident3,
    "min.x.notstress.ident4" = min.x.notstress.ident4, "max.x.notstress.ident4" = max.x.notstress.ident4,
    "step.notstress.ident4" = step.notstress.ident4
  )

  app_env$"average.vec" <- list(
    "gr.av.stress.scores1" = gr.av.stress.scores1,
    "gr.av.stress.scores2" = gr.av.stress.scores2,
    "gr.av.notstress.scores3" = gr.av.notstress.scores3,
    "gr.av.notstress.scores4" = gr.av.notstress.scores4
  )

  app_env$"obj" <- obj
  app_env$"plot.cluster.shiny" <- plot.cluster.shiny

  # Launch the Shiny app
  app_dir <- system.file("shiny", "GO.thresh", package = "gruffi")
  app_ui <- source(file.path(app_dir, "ui.R"),
                   local = new.env(parent = app_env),
                   echo = FALSE, keep.source = TRUE
  )$value

  # Call the server.R file
  app_server <- source(file.path(app_dir, "server.R"),
                       local = new.env(parent = app_env),
                       echo = FALSE, keep.source = TRUE
  )$value

  obj <- shiny::runApp(shiny::shinyApp(app_ui, app_server))
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Automated GO Term-based Thresholding without using a Shiny app.
#'
#' @description Automatically applies granular GO term-based thresholding to identify stressed cells
#' in a Seurat object without launching a Shiny app. This function calculates the
#' "granule average GO-scores" for individual GO terms and updates the object with a combined stress
#' identification based on these thresholds.
#' @param obj Seurat single cell object, with Granule Averages GO-Scores in metadata (ident1-4).
#' @param proposed.method The method for estimating thresholds, default: `"fitted"`.
#' Options are `"fitted"` or `"empirical"`.
#' @param quantile The quantile cutoff for threshold determination, default: `0.99`.
#' Options are `0.99` or `0.9`.
#' @param stress.ident1 Identifier for the first stress-related GO term.
#' @param stress.ident2 Identifier for the second stress-related GO term.
#' @param notstress.ident3 Identifier for the first non-stress-related GO term. Idents 3 and 4 act
#'  as a negative filter for stress identification.
#' @param notstress.ident4 Identifier for the second non-stress-related GO term, optional.
#' @param plot.results Boolean flag to control the plotting of results on UMAPs and histograms,
#' default: `TRUE`.
#'
#' @export
Auto.GO.thresh <- function(
    obj = combined.obj,
    proposed.method = c("fitted", "empirical")[1],
    quantile = c(.99, .9)[1],
    stress.ident1,
    stress.ident2,
    notstress.ident3,
    notstress.ident4 = NULL,
    plot.results = TRUE,
    ...) {
  meta <- obj@meta.data
  # Ensure that provided GO term identifiers exist in the metadata
  stopifnot(stress.ident1 %in% colnames(meta) | is.null(stress.ident1))
  stopifnot(stress.ident2 %in% colnames(meta) | is.null(stress.ident2))
  stopifnot(notstress.ident3 %in% colnames(meta) | is.null(notstress.ident3))
  stopifnot(notstress.ident4 %in% colnames(meta) | is.null(notstress.ident4))

  # Convert categorical scores to numeric for thresholding
  gr.av.stress.scores1 <- as.numeric(levels(meta[, stress.ident1]))
  gr.av.stress.scores2 <- as.numeric(levels(meta[, stress.ident2]))
  gr.av.notstress.scores3 <- as.numeric(levels(meta[, notstress.ident3]))
  gr.av.notstress.scores4 <- as.numeric(levels(meta[, notstress.ident4]))

  # Compute threshold proposals for each "granule average GO-score" using PlotNormAndSkew function
  thresh.stress.ident1 <- PlotNormAndSkew(gr.av.stress.scores1, q = quantile, tresholding = proposed.method, plot.hist = FALSE)
  thresh.stress.ident2 <- PlotNormAndSkew(gr.av.stress.scores2, q = quantile, tresholding = proposed.method, plot.hist = FALSE)
  thresh.notstress.ident3 <- PlotNormAndSkew(gr.av.notstress.scores3, q = quantile, tresholding = proposed.method, plot.hist = FALSE)
  thresh.notstress.ident4 <- PlotNormAndSkew(gr.av.notstress.scores4, q = quantile, tresholding = proposed.method, plot.hist = FALSE)

  # Store computed thresholds in the Seurat object's misc slot for later reference
  if (!is.null(stress.ident1)) obj@misc$gruffi$"thresh.stress.ident1" <- thresh.stress.ident1
  if (!is.null(stress.ident2)) obj@misc$gruffi$"thresh.stress.ident2" <- thresh.stress.ident2
  if (!is.null(notstress.ident3)) obj@misc$gruffi$"thresh.notstress.ident3" <- thresh.notstress.ident3
  if (!is.null(notstress.ident4)) obj@misc$gruffi$"thresh.notstress.ident4" <- thresh.notstress.ident4

  # Stress Filtering & Assignment
  # Cells are considered stressed if their scores are above the threshold for any of the stress identifiers
  # and not above the threshold for any of the non-stress identifiers
  if (!is.null(stress.ident1)) {
    gr.av.scores.1 <- meta[, stress.ident1]
    gr.av.scores.2 <- meta[, stress.ident2]

    # Check if the numeric score is greater than the threshold for stress
    i1.bool <- as.numeric(levels(gr.av.scores.1))[gr.av.scores.1] > thresh.stress.ident1
    if (!is.null(stress.ident2)) {
      i2.bool <- as.numeric(levels(gr.av.scores.2))[gr.av.scores.2] > thresh.stress.ident2
      stress.bool <- i1.bool | i2.bool # Combine both boolean vectors for stress determination
    } else {
      stress.bool <- i1.bool # Use only stress.ident1 for stress determination
    }
  } else {
    # Process only the second stress identifier if the first one is null
    if (!is.null(stress.ident2)) {
      i2.bool <- as.numeric(levels(gr.av.scores.2))[gr.av.scores.2] > thresh.stress.ident2
      stress.bool <- i2.bool
    }
  }

  # Process not stress identifiers similarly
  notstress.bool <- NULL
  gr.av.scores.3 <- meta[, notstress.ident3]
  gr.av.scores.4 <- meta[, notstress.ident4]

  # Determine not stressed cells based on the thresholds for notstress.ident3 and potentially notstress.ident4
  if (!is.null(notstress.ident3)) {
    i3.bool <- as.numeric(levels(gr.av.scores.3))[gr.av.scores.3] > thresh.notstress.ident3
    if (!is.null(notstress.ident4)) {
      i4.bool <- as.numeric(levels(gr.av.scores.4))[gr.av.scores.4] > thresh.notstress.ident4
      notstress.bool <- i3.bool | i4.bool # Combine boolean vectors for notstress.ident3 and notstress.ident4
    } else {
      notstress.bool <- i3.bool # Use only notstress.ident3 for not stress determination
    }
  } else { # aka IF notstress.ident3 is null
    if (!is.null(notstress.ident4)) {
      i4.bool <- as.numeric(levels(gr.av.scores.4))[gr.av.scores.4] > thresh.notstress.ident4
      notstress.bool <- i4.bool
    }
  }

  # Final assignment of stressed cells
  if (!is.null(notstress.bool)) {
    obj$"is.Stressed" <- stress.bool & !notstress.bool
  } else {
    obj$"is.Stressed" <- stress.bool
  }

  if (plot.results) {
    message("plot.results not yet implemented")
  }

  message("Seurat object now contains:")
  message(" --  A new metadata column: `is.Stressed`:")
  message(" --  Threshold values for granule average scores: `obj@misc$gruffi$...`:")
  print(pc_TRUE(obj$is.Stressed, suffix = "stressed cells", NumberAndPC = TRUE))
  print(obj@misc$"gruffi")

  return(obj)
}


# _________________________________________________________________________________________________
#' @title Filter Stressed Cells from a Seurat Object
#'
#' @description Identifies and filters stressed cells based on specified GO terms and a quantile threshold.
#' It supports optional plotting of exclusion results and saving of modified datasets.
#' @param obj Seurat single cell object, Default: combined.obj
#' @param res Clustering resolution, Default: 'integrated_snn_res.30'
#' @param quantile.thr Quantile threshold to cutoff stressed cells, Default: 0.9
#' @param GOterms GO-terms to use for filtering, Default: c(`glycolytic process` = "GO:0006096", `response to endoplasmic reticulum stress` = "GO:0034976")
#' @param direction filtering direction, above or below, Default: 'above'
#' @param saveRDS Logical indicating if the filtered object should be saved as an RDS file, default is TRUE.
#' @param saveRDS.Removed Logical indicating if a separate RDS file for removed cells should be saved, default is FALSE.
#' @param PlotExclusionByEachScore Logical for plotting exclusion results by each score, default is TRUE.
#' @param PlotSingleCellExclusion Logical for plotting single cell exclusion results, default is TRUE.
#' @param GranuleExclusionScatterPlot Logical for plotting a scatter plot of granule exclusion, default is TRUE.
#' @seealso
#'  \code{\link[Stringendo]{iprint}}, \code{\link[Stringendo]{percentage_formatter}}
#'  \code{\link[Seurat.utils]{calc.cluster.averages}}
#'  \code{\link[MarkdownHelpers]{filter_HP}}, \code{\link[MarkdownHelpers]{llprint}}
#'  \code{\link[dplyr]{reexports}}
#' @export
#'
#' @importFrom CodeAndRoll2 matrix.fromNames which_names
#' @importFrom dplyr tibble
#' @importFrom ggExpress qscatter qhistogram
#' @importFrom MarkdownHelpers filter_HP llprint
#' @importFrom Seurat.utils clUMAP calc.cluster.averages isave.RDS
#' @importFrom Stringendo iprint percentage_formatter
FilterStressedCells <- function(
    obj = combined.obj,
    res = "integrated_snn_res.30",
    quantile.thr = 0.9,
    GOterms = c("glycolytic process" = "GO:0006096", "response to endoplasmic reticulum stress" = "GO:0034976"),
    direction = "above",
    saveRDS = TRUE, saveRDS.Removed = FALSE,
    PlotExclusionByEachScore = TRUE,
    PlotSingleCellExclusion = TRUE,
    GranuleExclusionScatterPlot = TRUE) {
  # Check arguments ____________________________________________________________
  Meta <- obj@meta.data
  MetaVars <- colnames(Meta)
  if (!res %in% MetaVars) {
    Stringendo::iprint("res", res, "is not found in the object.")
  }

  ScoreNames <- .convert.GO_term.2.score(GOterms)
  if (!all(ScoreNames %in% MetaVars)) {
    Stringendo::iprint(
      "Some of the GO-term scores were not found in the object:", ScoreNames,
      "Please call first: PlotGoTermScores(), or the actual GetGOTerms()"
    )
  }

  cells.2.granules <- obj[[res]][, 1]

  # Exclude granules ____________________________________________________________
  mScoresFiltPass <- mScores <- CodeAndRoll2::matrix.fromNames(fill = NaN, rowname_vec = sort(unique(Meta[, res])), colname_vec = make.names(GOterms))
  for (i in 1:length(GOterms)) {
    scoreX <- as.character(ScoreNames[i])
    print(scoreX)
    scoreNameX <- names(ScoreNames)[i]
    mScoresFiltPass[, i] <- Seurat.utils::calc.cluster.averages(
      col_name = scoreX, split_by = res, quantile.thr = quantile.thr,
      histogram = TRUE, subtitle = names(ScoreNames)[i],
      filter = direction, obj = obj
    )

    if (GranuleExclusionScatterPlot) {
      mScores[, i] <- Seurat.utils::calc.cluster.averages(
        col_name = scoreX, split_by = res, quantile.thr = quantile.thr,
        plot.UMAP.too = FALSE, plotit = FALSE,
        filter = FALSE
      )
    }

    if (PlotExclusionByEachScore) {
      Seurat.utils::clUMAP(
        ident = res, highlight.clusters = CodeAndRoll2::which_names(mScoresFiltPass[, i]),
        label = FALSE, title = paste0("Stressed cells removed by ", scoreX),
        plotname = paste0("Stressed cells removed by ", scoreX),
        sub = scoreNameX,
        sizes.highlight = .5, raster = FALSE
      )
    }
  }

  # PlotSingleCellExclusion ____________________________________________________________
  if (PlotSingleCellExclusion) {
    mSC_ScoresFiltPass <- CodeAndRoll2::matrix.fromNames(fill = NaN, rowname_vec = rownames(Meta), colname_vec = make.names(GOterms))
    for (i in 1:length(GOterms)) {
      scoreX <- as.character(ScoreNames[i])
      print(scoreX)
      scoreNameX <- names(ScoreNames)[i]
      ScoreValuesSC <- Meta[, scoreX]
      mSC_ScoresFiltPass[, i] <- MarkdownHelpers::filter_HP(
        numeric_vector = ScoreValuesSC,
        threshold = stats::quantile(x = ScoreValuesSC, quantile.thr), breaks = 100
      )
    }
    cells.remove.SC <- rowSums(mSC_ScoresFiltPass) > 0
    table(cells.remove.SC)
    sc.filtered <- CodeAndRoll2::which_names(cells.remove.SC)
    sc.kept <- CodeAndRoll2::which_names(!cells.remove.SC)

    Seurat.utils::clUMAP(
      cells.highlight = sc.filtered, label = FALSE,
      title = "Single-cell filtering is not robust to remove stressed cells",
      plotname = "Single-cell filtering is not robust to remove stressed cells",
      suffix = "Stress.Filtering", sizes.highlight = .5, raster = FALSE
    )


    MeanZ.ER.stress <- c(
      "mean.stressed" = mean(Meta[sc.filtered, scoreX]),
      "mean.kept" = mean(Meta[sc.kept, scoreX])
    )
    qbarplot(MeanZ.ER.stress, ylab = scoreNameX, xlab.angle = 45, xlab = "", col = as.logical(0:1))

    # if (TRUE) {
    #   obj.Removed.by.SingleCells <- subset(x = obj, cells = sc.filtered) # remove stressed
    #   Seurat.utils::isave.RDS(obj.Removed.by.SingleCells, inOutDir = TRUE, suffix = "Removed")
    # }
  } # if (PlotSingleCellExclusion)

  # Plot excluded cells ____________________________________________________________
  PASS <- (rowSums(mScoresFiltPass) == 0)
  granules.excluded <- CodeAndRoll2::which_names(!PASS)
  Seurat.utils::clUMAP(ident = res, highlight.clusters = granules.excluded, label = FALSE, title = "Stressed cells removed", suffix = "Stress.Filtering", sizes.highlight = .5, raster = FALSE)

  if (GranuleExclusionScatterPlot) {
    Av.GO.Scores <- dplyr::tibble(
      "Average ER-stress score" = mScores[, 2],
      "Average Glycolysis score" = mScores[, 1],
      "Stressed" = !PASS,
      "name" = rownames(mScores)
    )

    if (FALSE) colnames(Av.GO.Scores)[1:2] <- names(GOterms)

    ggExpress::qscatter(Av.GO.Scores,
                        cols = "Stressed", w = 6, h = 6,
                        vline = stats::quantile(mScores[, 2], quantile.thr),
                        hline = stats::quantile(mScores[, 1], quantile.thr),
                        title = "Groups of Stressed cells have high scores.",
                        subtitle = "Thresholded at 90th percentile",
                        label = "name",
                        repel = TRUE,
                        label.rectangle = TRUE
    )
  }


  # Exclude cells & subset ____________________________________________________________

  cells.discard <- which(cells.2.granules %in% granules.excluded)
  cells.keep <- which(!(cells.2.granules %in% granules.excluded))
  MarkdownHelpers::llprint(Stringendo::percentage_formatter(length(cells.keep) / length(cells.2.granules)), "cells kept.")

  obj.noStress <- subset(x = obj, cells = cells.keep) # remove stressed
  if (saveRDS) Seurat.utils::isave.RDS(obj.noStress, inOutDir = TRUE, suffix = "Cleaned")
  if (saveRDS.Removed) {
    obj.RemovedCells <- subset(x = obj, cells = cells.discard) # remove stressed
    Seurat.utils::isave.RDS(obj.RemovedCells, inOutDir = TRUE, suffix = "Removed")
  }
  return(obj.noStress)
}






# _____________________________________________________________________________________________ ----
# 5. Visualization ---------------------------------------------------------------------------


# _________________________________________________________________________________________________
#' @title PlotGoTermScores
#'
#' @description Automates retrieving, processing and plotting GO term based gene scores.
#' @param obj Seurat single cell object, Default: combined.obj
#' @param data.base.access Which tool to use to access gene list databases? Options:
#' biomaRt (deault) or AnnotationDbi.
#' @param mirror Ensembl mirror to use, helpful if the default connection fails with
#' default settings, mirror can be specified. Default: 'NULL'
#' @param desc GO-score description on the plot, Default: ""
#' @param only.draw.plot Only show GO term score umap plot, Default: FALSE
#' @param openBrowser Open website for GO term? Default: FALSE
#' @param plot.each.gene Plot each gene's expression on a combined umap? Default: FALSE
#' @param save.UMAP Save umap into a file. Default: TRUE
#' @param verbose Verbose output? Default: TRUE
#' @param GO GO-term; Default: 'GO:0009651'
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @seealso
#'  \code{\link[Seurat]{reexports}}
#'  \code{\link[Stringendo]{iprint}}
#' @export
#' @importFrom Seurat DefaultAssay
#' @importFrom Stringendo iprint ppp
#' @importFrom Seurat.utils clUMAP multiFeaturePlot.A4

PlotGoTermScores <- function(
    obj = combined.obj,
    data.base.access = c("biomaRt", "AnnotationDbi")[1],
    mirror = NULL,
    only.draw.plot = FALSE,
    openBrowser = FALSE,
    plot.each.gene = FALSE,
    desc = "",
    save.UMAP = TRUE,
    verbose = TRUE,
    GO = "GO:0009651",
    ...) {

  backup.assay <- Seurat::DefaultAssay(obj)

  if (grepl(x = GO, pattern = "^GO:", perl = TRUE)) {
    # stopif(condition = grepl(pattern = "Score.", x = GO), message = print("Provide a simple GO-term, like: GO:0009651 - not Score.GO.0006096"))
    stopifnot(
      "Provide a simple GO-term, like: GO:0009651 - not Score.GO.0006096" =
        !grepl(pattern = "Score.", x = GO)
    )

    GO.wDot <- make.names(GO)
    ScoreName <- paste0("Score.", GO.wDot)
    print(ScoreName)

  } else {
    message("Assuming you provided direclty a score name in the 'GO' parameter ", ScoreName)
    ScoreName <- GO
  }


  if (Seurat::DefaultAssay(obj) != "RNA") {
    message("For GO score computation assay is set to RNA. It will be reset to Seurat::DefaultAssay() afterwards.")
    Seurat::DefaultAssay(obj) <- "RNA"
  }

  if (ScoreName %in% colnames(obj@meta.data)) {
    message("GENE SCORE FOUND IN @meta.data")
  } else {
    message(ScoreName, " not found. Need to call GetGOTerms(). set only.draw.plot = FALSE.")
  }

  if (!only.draw.plot) {
    message("\nGENE SCORE WILL NOW BE SET / OVERWRITTEN")
    obj@meta.data[, ScoreName] <- NULL

    obj <- GetGOTerms(obj = obj, GO = GO, web.open = openBrowser, data.base.access = data.base.access, mirror = mirror)
    GO.genes <- obj@misc$gruffi$GO[[GO.wDot]]
    if (verbose) print(head(GO.genes))
    obj <- AddGOScore(obj = obj, GO = GO)
  }

  plot <- FeaturePlotSaveGO(obj = obj, GO.score = ScoreName, save.plot = save.UMAP, name_desc = desc, ...)
  if (plot.each.gene) Seurat.utils::multiFeaturePlot.A4(obj = obj, list.of.genes = GO.genes, foldername = Stringendo::ppp(GO.wDot, "UMAPs"))
  Seurat::DefaultAssay(obj) <- backup.assay

  if (only.draw.plot) {
    return(plot)
  } else {
    return(obj)
  }
}


# _________________________________________________________________________________________________
#' @title FeaturePlotSaveCustomScore
#' @description Plot and save a Seurat FeaturePlot with a custom score in the meta data.
#' @param obj Seurat single cell object, Default: combined.obj
#' @param genes Charcter vector of genes, Default: ""
#' @param name_desc Descriptive name to suffix plot and file name, Default: NULL
#' @param h height, Default: 7
#' @param PNG Save as .png? Default: TRUE
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @seealso
#'  \code{\link[Seurat]{FeaturePlot}}
#'  \code{\link[ggplot2]{labs}}
#'  \code{\link[MarkdownHelpers]{ww.FnP_parser}}
#'  \code{\link[cowplot]{save_plot}}
#' @export
#' @importFrom Seurat FeaturePlot
#' @importFrom ggplot2 labs
#' @importFrom MarkdownHelpers ww.FnP_parser
#' @importFrom cowplot save_plot

FeaturePlotSaveCustomScore <- function(obj = combined.obj, genes = "", name_desc = NULL, h = 7, PNG = TRUE, ...) { # Plot and save a FeaturePlot, e.g. showing gene set scores.
  ScoreName <- paste0("Score.", substitute(genes))

  ggplot.obj <-
    Seurat::FeaturePlot(obj, features = ScoreName, min.cutoff = "q05", max.cutoff = "q95", reduction = "umap", ...) +
    ggplot2::labs(title = paste(ScoreName, name_desc), caption = paste("Score calc. from", length(genes), "expr. genes ."))
  pname <- paste0("FeaturePlot.", ScoreName)
  fname <- MarkdownHelpers::ww.FnP_parser(Stringendo::kpp(pname, name_desc), if (PNG) "png" else "pdf")
  cowplot::save_plot(filename = fname, plot = ggplot.obj, base_height = h)
  ggplot.obj
}

# _________________________________________________________________________________________________
#' @title FeaturePlotSaveGO
#'
#' @description Plot and save a Seurat FeaturePlot from a GO-score.
#' @param obj Seurat single cell object, Default: combined.obj
#' @param GO.score Metadata column name for a GO-term scoreDefault: 'Score.GO.0034976'
#' @param name_desc Descriptive name to suffix plot and file name, Default: NULL
#' @param save.plot Save plot into a file. Default: TRUE
#' @param title_ Plot title, Default: paste(GO.score, name_desc)
#' @param h Plot height, Default: 7
#' @param PNG Save as .png? Default: TRUE
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @seealso
#'  \code{\link[Seurat]{FeaturePlot}}
#'  \code{\link[ggplot2]{labs}}
#'  \code{\link[MarkdownHelpers]{ww.FnP_parser}}
#'  \code{\link[cowplot]{save_plot}}
#' @export
#' @importFrom Seurat FeaturePlot
#' @importFrom ggplot2 labs
#' @importFrom MarkdownHelpers ww.FnP_parser
#' @importFrom cowplot save_plot

FeaturePlotSaveGO <- function(
    obj = combined.obj,
    GO.score = "Score.GO.0034976",
    name_desc = NULL,
    save.plot = TRUE,
    title_ = paste(GO.score, name_desc),
    h = 7, PNG = TRUE, ...) {

  proper.GO <- paste(stringr::str_split_fixed(string = GO.score, pattern = "\\.", n = 3)[2:3], collapse = ":")
  (genes.GO <- obj@misc$gruffi$GO[[make.names(proper.GO)]])

  ggplot.obj <-
    Seurat::FeaturePlot(obj, features = GO.score, min.cutoff = "q05", max.cutoff = "q95", reduction = "umap", ...) +
    ggplot2::labs(title = title_, caption = paste("Score calc. from", length(genes.GO), "expr. genes from BioMart.", paste0("https://www.ebi.ac.uk/QuickGO/search/", proper.GO)))
  pname <- paste0("FeaturePlot.", (GO.score))
  fname <- MarkdownHelpers::ww.FnP_parser(Stringendo::kpp(pname, name_desc), if (PNG) "png" else "pdf")
  if (save.plot) {
    cowplot::save_plot(filename = fname, plot = ggplot.obj, base_height = h)
  }
  ggplot.obj
}


# _________________________________________________________________________________________________
#' @title PlotClustSizeDistr
#'
#' @description Plot cluster size distribution.
#' @param obj Seurat single cell object, Default: combined.obj
#' @param category Clustering or any categorical metadata column name.
#' Default: GetGruffiClusteringName(obj)
#' @param plot Plot results? Default: TRUE
#' @param thr.hist Plot histogram or a barplot? If less than X categories, it draws a barplot. If more, a histogram, Default: 30
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @seealso
#'  \code{\link[ggExpress]{qbarplot}}, \code{\link[ggExpress]{qhistogram}}
#' @export
#' @importFrom ggExpress qbarplot qhistogram

PlotClustSizeDistr <- function(obj = combined.obj,
                               category = GetGruffiClusteringName(obj),
                               plot = TRUE,
                               thr.hist = 30,
                               ...) {
  .Deprecated("Seurat.utils::plotClustSizeDistr")

  clust.size.distr <- table(obj@meta.data[, category])
  print(clust.size.distr)
  resX <- gsub(pattern = ".*res\\.", replacement = "", x = category)

  if (length(clust.size.distr) < thr.hist) {
    ggExpress::qbarplot(clust.size.distr,
                        plotname = Stringendo::ppp("clust.size.distr", (category)),
                        subtitle = paste("Nr.clusters at res.", resX, ":", length(clust.size.distr), " | CV:", percentage_formatter(CodeAndRoll2::cv(clust.size.distr))),
                        ...
    )
  } else {
    ggExpress::qhistogram(
      vec = clust.size.distr, plotname = Stringendo::ppp("clust.size.distr", (category)),
      subtitle = paste("Nr.clusters at res.", resX, ":", length(clust.size.distr), " | CV:", percentage_formatter(CodeAndRoll2::cv(clust.size.distr))),
      ...
    )
  }
}

# _________________________________________________________________________________________________
#' @title PlotNormAndSkew
#'
#' @description Plot normal distribution and skew.
#' @param x Distribution
#' @param q Quantile
#' @param tresholding Tresholding calculation method, Default: c("fitted", "empirical")[1]
#' @param plot.hist Draw a histogram? Default: TRUE
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @export
PlotNormAndSkew <- function(x, q,
                            tresholding = c("fitted", "empirical")[1],
                            plot.hist = TRUE, ...) {
  m <- median(x)
  stde_skewed <- CalcStandDevSkewedDistr(x, mean_x = m)

  thresh <-
    if (tresholding == "fitted") {
      stats::qnorm(q, mean = m, sd = stde_skewed)
    } else if (tresholding == "empirical") {
      stats::quantile(x, q)
    } else {
      print("Unknown tresholding method")
    }

  if (plot.hist) {
    Score.threshold.estimate <- x
    sb <- print(paste(
      "The computed", tresholding, paste0("q", q), "threshold is:",
      signif(thresh, digits = 3), "\nValues above thr:", CodeAndRoll2::pc_TRUE(x > thresh)
    ))
    qhistogram(Score.threshold.estimate,
               vline = thresh, subtitle = sb, xlab = "Score",
               plotname = Stringendo::ppp("Score.threshold.est", substitute(x), tresholding), suffix = Stringendo::ppp("q", q)
    )
  }
  return(thresh)
}





# _________________________________________________________________________________________________
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
GrScoreUMAP <- function(obj = combined.obj,
                        colname = "RNA_snn_res.6.reassigned_cl.av_GO:0042063",
                        miscname = "thresh.stress.ident1",
                        auto = TRUE,
                        ...) {
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

# _________________________________________________________________________________________________
#' @title Histogram Visualization of Granule Scores with Threshold
#'
#' @description Plots a histogram of granule scores from a specified column within a Seurat object,
#' with a vertical line denoting a threshold. The function allows for an optional
#' inversion of filtering logic based on the `miscname` and `auto` parameters.
#'
#' @param obj A Seurat object containing metadata with granule scores and miscellaneous information.
#' Defaults to `combined.obj`.
#' @param colname The name of the column within `obj@meta.data` containing granule scores.
#' Typically something like is `'RNA_snn_res.6.reassigned_cl.av_GO:0042063'`.
#' @param miscname The key within `obj@misc$gruffi` for retrieving the threshold value.
#' Default is `'thresh.stress.ident1'`.
#' @param auto A logical flag indicating whether to automatically invert the filtering logic
#' for a specific `miscname` condition. Default is TRUE.
#' @param ... Additional parameters to be passed to `ggExpress::qhistogram`.
#'
#' @details The function first verifies the structure and content of the provided `obj`, ensuring it
#' contains the necessary components for threshold retrieval and score extraction. It then retrieves
#' the threshold and uses `CodeAndRoll2::as.numeric.wNames.character` to convert the specified granule
#' scores to numeric. The histogram is plotted using `ggExpress::qhistogram`, with customization options
#' available through additional arguments. The `auto` parameter's effect is noted, particularly for the
#' `thresh.notstress.ident3` condition, although it primarily impacts plotting logic interpretation.
#'
#' @return The function directly calls `ggExpress::qhistogram` for plotting, hence does not return a value
#' but generates a histogram plot as a side effect.
#'
#' @examples
#' GrScoreHistogram(
#'   obj = combined.obj,
#'   colname = "RNA_snn_res.6.reassigned_cl.av_GO:0042063",
#'   miscname = "thresh.stress.ident1"
#' )
#' @export
#'
#' @importFrom CodeAndRoll2 as.numeric.wNames.character
#' @importFrom ggExpress qhistogram
GrScoreHistogram <- function(obj = combined.obj,
                             colname = i1,
                             miscname = "thresh.stress.ident1",
                             auto = TRUE,
                             ...) {
  # Argument assertions
  stopifnot(
    inherits(obj, "Seurat"), !is.null(obj@misc$gruffi),
    miscname %in% names(obj@misc$gruffi),
    colname %in% names(obj@meta.data)
  )

  # Retrieve threshold from object
  thr <- obj@misc$gruffi[[miscname]]

  # Convert granule scores to numeric
  granule_scores <- CodeAndRoll2::as.numeric.wNames.character(obj@meta.data[[colname]], verbose = FALSE)


  # Apply threshold to create a logical vector for identification
  obj[["pass"]] <- (granule_scores < thr)

  components <- strsplit(colname, "_cl\\.av_")[[1]]

  nr.clust <- NaN
  try(nr.clust <- length(unique(obj@meta.data[[components[1]]])), silent = TRUE)

  subt <- paste(
    "threshold:", thr, " | Pass: ", pc_TRUE(obj@meta.data[["pass"]]), " | ",
    ncol(obj), "cells |", nr.clust, "clusters."
  )

  colX <- if (miscname == "thresh.notstress.ident3" & auto) 1 else -1

  # Plot histogram with ggExpress::qhistogram
  pobj <- ggExpress::qhistogram(granule_scores,
                                sub = subt, palette_use = "npg",
                                plotname = paste("Granule Thresholding", make.names(components[2])),
                                xlab = "Granule Median Score", ylab = "Nr. of Granules",
                                vline = thr, filtercol = colX, ...
  )
  print(pobj)

  # Optionally, handle the 'auto' and 'miscname' parameters for additional logic
  if (miscname == "thresh.notstress.ident3" & auto) {
    message("Note: 'thresh.notstress.ident3' condition detected, but it affects UMAP plotting logic, not histogram.")
  }
}



# _________________________________________________________________________________________________
#' @title ClusterUMAPthresholding
#' @description clUMAP with threshold value.
#' @param q.meta.col Numeric metadata column name, Default: 'Score.GO.0034976'
#' @param c.meta.col Clustering basis to split the data (Categorical metadata column name), Default: 'integrated_snn_res.30'
#' @param quantile Quantile, Default: 0.95
#' @param absolute.cutoff Absolute cutoff, Default: NULL
#' @param obj Seurat single cell object, Default: combined.obj
#' @param plot.barplot Draw a barplot? Default: FALSE
#' @param subt Plot Subtitle, Default: 'response to ER stress'
#' @param plotUMAP Draw a UMAP? Default: TRUE
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @seealso
#'  \code{\link[Seurat.utils]{calc.cluster.averages}}
#' @export
#' @importFrom Seurat.utils calc.cluster.averages

ClusterUMAPthresholding <- function(
    q.meta.col = "Score.GO.0034976", c.meta.col = "integrated_snn_res.30",
    quantile = .95, absolute.cutoff = NULL,
    obj = combined.obj, plot.barplot = FALSE,
    subt = "response to ER stress", plotUMAP = TRUE,
    ...) {
  clusters.keep <- Seurat.utils::calc.cluster.averages(
    col_name = q.meta.col, split_by = c.meta.col, absolute.thr = absolute.cutoff,
    quantile.thr = quantile, plotit = plot.barplot
  )

  cells.2.granules <- as.character(obj[[c.meta.col]][, 1])

  cluster.nrs.discard <- clusters.keep[which(!clusters.keep)]
  filtered.out <- gsub(pattern = "^cl\\.", replacement = "", x = names(cluster.nrs.discard))
  is.2.highlight <- cells.2.granules %in% filtered.out
  idx.highlight.these <- which(is.2.highlight)
  cells.highlight <- colnames(obj)[idx.highlight.these]

  ttl <- paste("Cells above", "quantile", quantile, "in", q.meta.col)
  if (plotUMAP) {
    Seurat.utils::clUMAP(
      obj = obj, cells.highlight = cells.highlight,
      plotname = Stringendo::ppp("UMAP.thresholding_q", q.meta.col, quantile, c.meta.col),
      title = ttl, sub = c.meta.col, raster = FALSE, label.cex = 5, ...
    ) # , shape.by = c.meta.col
  } else {
    is.2.highlight
  }
}






# _____________________________________________________________________________________________ ----
# 6. Stats and Math ---------------------------------------------------------------------------


# _________________________________________________________________________________________________
#' @title CalcTranscriptomePercentageGO
#' @description Calculate the percentage of transcriptome, for a given GO-term gene set
#' already stored`in obj@misc$gruffi$GO[[GO.score]].
#' @param obj Seurat single cell object, Default: combined.obj
#' @param GO.score GO Term. String of form "GO:xxxxxxx", Default: 'GO.0061621'
#' @seealso
#'  \code{\link[Matrix]{colSums}}
#'  \code{\link[Seurat]{reexports}}
#' @export
#' @importFrom Matrix colSums
#' @importFrom Seurat GetAssayData

CalcTranscriptomePercentageGO <- function(obj = combined.obj, GO.score = "GO.0061621") {
  total_expr <- Matrix::colSums(Seurat::GetAssayData(object = obj))
  Matrix::colSums(obj[obj@misc$gruffi$GO[[GO.score]], ]) / total_expr
}




# _________________________________________________________________________________________________
#' @title CalcTranscriptomePercentage
#' @description Calculate the percentage of transcriptome, for a given gene set.
#' @param obj Seurat single cell object, Default: combined.obj
#' @param genes Charcter vector of genes, Default: genes.GO.0061621.can.glyc
#' @seealso
#'  \code{\link[Matrix]{colSums}}
#'  \code{\link[Seurat]{reexports}}
#' @export
#' @importFrom Matrix colSums
#' @importFrom Seurat GetAssayData

CalcTranscriptomePercentage <- function(obj = combined.obj, genes = genes.GO.0061621.can.glyc) {
  total_expr <- Matrix::colSums(Seurat::GetAssayData(object = obj))
  Matrix::colSums(obj[genes, ]) / total_expr
}



# _________________________________________________________________________________________________
#' @title Calculate Standard Deviation for Skewed Distributions
#'
#' @description Computes a modified standard deviation for skewed distributions by mirroring data points about the mean.
#' This approach aims to reflect the skewed part of the distribution to better understand its spread.
#'
#' @param x A numeric vector representing the distribution.
#' @param mean_x The central tendency of `x`, used to mirror the distribution.
#' Default is the mean of `x`.
#' @param plothist Logical indicating whether to display histograms of the original, mirrored,
#' and combined distributions. Default is `FALSE`.
#' @param breaks The number of bins to use in the histogram if plotted.
#' Default is 50.
#'
#' @return The estimated standard deviation for the skewed distribution.
#'
#' @examples
#' set.seed(123)
#' skewed_data <- rexp(100, rate = 0.1)
#' # Calculate modified standard deviation without plotting
#' modified_sd <- stand_dev_skewed(skewed_data)
#'
#' # Calculate and plot the distribution and its mirrored version
#' stand_dev_skewed(skewed_data, plothist = TRUE)
#'
#' @export

CalcStandDevSkewedDistr <- function(x, mean_x = mean(x),
                                    plothist = FALSE, breaks = 50) {

  x_lower_eq <- x[x <= mean_x]
  x_lower <- x[x < mean_x]
  x_mirror <- c(x_lower_eq, abs(mean_x - x_lower) + mean_x)

  if (plothist) {
    hist(x_lower, xlim = range(x), breaks = breaks)
    hist(x, xlim = range(x), breaks = breaks)
    hist(x_mirror, xlim = range(x), breaks = breaks)
  }

  res <- sum((x_mirror - mean_x)^2) / (length(x_mirror) - 1)
  return(sqrt(res))
}



# _________________________________________________________________________________________________
#' @title CalcClusterAverages_Gruffi
#'
#' @description Calculate granule (cluster) averages scores.
#' @param col_name Numeric metadata column name, Default: 'Score.GO.0006096'
#' @param obj Seurat single cell object, Default: combined.obj
#' @param split_by Clustering to split by (categorical metadata column name), Default: Seurat.utils::GetClusteringRuns(obj)[1]
#' @param stat How to caluclate the central tendency? Default: c("mean", "median", "normalized.mean", "normalized.median")[3]
#' @param scale.zscore Scale as zscore?, Default: FALSE
#' @param quantile.thr Quantile threshold for cutoff, Default: 0.99
#' @param absolute.thr Use an absolute threshold? Default: FALSE
#' @param plotit Draw a plot? Default: FALSE
#' @param save.plot Save plot into a file. Default: FALSE
#' @param return.plot Return the plot as output of the function? Default: FALSE
#' @param max.bin.plot What is considered too many clusters for a histogram to plot? Default: 200
#' @param simplify Simplify output? Default: TRUE
#' @param suffix What suffix to append at the end of the plotname/file name? Default: NULL
#' @param ylab.text Y-axis label. Default: paste("Cluster", stat, "score")
#' @param title Plot title, Default: paste("Cluster", stat, col_name)
#' @param subtitle Plot subtitle, Default: NULL
#' @param width Plot width, Default: 8
#' @param height Plot height, Default: 6
#' @param xlb X-axis label, Default: if (absolute.thr) paste("Threshold at", absolute.thr) else paste("Black lines: ",
#'    Stringendo::kppd(Stringendo::percentage_formatter(quantile.thr)), "quantiles |",
#'    "Cl. >", Stringendo::percentage_formatter(quantile.thr),
#'    "are highlighted. |", split_by)
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @seealso
#'  \code{\link[Stringendo]{percentage_formatter}}, \code{\link[Stringendo]{iprint}}
#'  \code{\link[dplyr]{select_all}}, \code{\link[dplyr]{group_by_all}}, \code{\link[dplyr]{summarise}}, \code{\link[dplyr]{context}}
#'  \code{\link[rlang]{sym}}
#'  \code{\link[CodeAndRoll2]{sem}}, \code{\link[CodeAndRoll2]{sortbyitsnames}}
#' @export
#' @importFrom Stringendo percentage_formatter iprint
#' @importFrom dplyr select_at group_by_at summarize n
#' @importFrom rlang sym
#' @importFrom CodeAndRoll2 sem sortbyitsnames

CalcClusterAverages_Gruffi <- function(
    col_name = "Score.GO.0006096",
    obj = combined.obj,
    split_by = Seurat.utils::GetClusteringRuns(obj)[1],
    stat = c("mean", "median", "normalized.mean", "normalized.median")[3],
    scale.zscore = FALSE,
    quantile.thr = 0.99,
    absolute.thr = FALSE,
    plotit = FALSE,
    save.plot = FALSE,
    return.plot = FALSE,
    max.bin.plot = 200,
    simplify = TRUE,
    suffix = NULL,
    ylab.text = paste("Cluster", stat, "score"),
    title = paste("Cluster", stat, col_name),
    subtitle = NULL,
    width = 8, height = 6,
    xlb = if (absolute.thr) {
      paste("Threshold at", absolute.thr)
    } else {
      paste(
        "Black lines: ", Stringendo::kppd(Stringendo::percentage_formatter(quantile.thr)), "quantiles |",
        "Cl. >", Stringendo::percentage_formatter(quantile.thr), "are highlighted. |", split_by
      )
    },
    ...) {
  Stringendo::iprint(substitute(obj), "split by", split_by)

  df.summary <-
    obj@meta.data %>%
    dplyr::select_at(c(col_name, split_by)) %>%
    dplyr::group_by_at(split_by) %>%
    dplyr::summarize(
      "nr.cells" = dplyr::n(),
      "mean" = mean(!!rlang::sym(col_name), na.rm = TRUE),
      "SEM" = CodeAndRoll2::sem(!!rlang::sym(col_name), na.rm = TRUE),
      "normalized.mean" = mean(!!rlang::sym(col_name), na.rm = TRUE),
      "median" = median(!!rlang::sym(col_name), na.rm = TRUE),
      "SE.median" = 1.2533 * CodeAndRoll2::sem(!!rlang::sym(col_name), na.rm = TRUE),
      "normalized.median" = median(!!rlang::sym(col_name), na.rm = TRUE)
    )

  df.summary$normalized.mean <- sqrt(df.summary$nr.cells) * (df.summary$normalized.mean - median(df.summary$normalized.mean))
  df.summary$normalized.median <- sqrt(df.summary$nr.cells) * (df.summary$normalized.median - median(df.summary$normalized.median))

  if (simplify) {
    av.score <- df.summary[[stat]]
    names(av.score) <- df.summary[[1]]
    av.score <- CodeAndRoll2::sortbyitsnames(av.score)

    if (scale.zscore) av.score <- (scale(av.score)[, 1])

    cutoff <- if (absolute.thr) absolute.thr else stats::quantile(av.score, quantile.thr)

    if (plotit) {
      if (length(av.score) > max.bin.plot) {
        print("Too many clusters for histogram plot.")
      } else {
        p <- qbarplot(
          vec = av.score, save = FALSE,
          hline = cutoff,
          plotname = title,
          suffix = quantile.thr,
          subtitle = subtitle,
          ylab = ylab.text,
          xlab = xlb # Abused
          , xlab.angle = 45
          # , ...
        )
        print(p)
        if (save.plot) {
          title_ <- Stringendo::ppp(title, suffix, Stringendo::flag.nameiftrue(scale.zscore))
          ggExpress::qqSave(ggobj = p, title = title_, fname = Stringendo::ppp(title_, split_by, "png"), w = width, h = height)
        }
        if (return.plot) {
          return(p)
        }
      }
    }
    return(av.score)
  } else {
    return(df.summary)
  }
}




# _____________________________________________________________________________________________ ----
# 7. Helpers  ---------------------------------------------------------------------------


#' @title Retrieve First Clustering Run or First Matching Pattern Run
#'
#' @description Fetches the first clustering run from a Seurat object,
#' optionally filtered by a specific pattern (e.g., ".reassigned").
#'
#' @param obj A Seurat object with clustering results in meta.data.
#' @param pattern A character string pattern to filter clustering runs.
#' Default is ".reassigned".
#'
#' @return The name of the first clustering run that matches the pattern
#' if any; otherwise, the name of the first available clustering run.
#'
#' @examples
#' # Assume `seuratObj` is your Seurat object
#' GetGruffiClusteringName(seuratObj)
#' GetGruffiClusteringName(seuratObj, pattern = "someOtherPattern")
#'
#' @export
#' @importFrom Stringendo iprint
GetGruffiClusteringName <- function(obj, pattern = ".reassigned",
                                    granule.res.slot = "optimal.granule.res") {

  res <- obj@misc$gruffi[[granule.res.slot]]

  if (is.null(res)) {
    message("obj@misc$gruffi$optimal.granule.res not found. Searching for pattern in meta.data: *", pattern)

    # Retrieve all clustering runs from the Seurat object
    clusteringRuns <- Seurat.utils::GetClusteringRuns(obj)

    # Check for clustering runs that match the given pattern
    matchingRuns <- clusteringRuns[grepl(pattern, clusteringRuns)]

    # Return the first matching run if any exist, otherwise return the first clustering run
    if (length(matchingRuns) == 1) {

      if (length(matchingRuns) > 1) {
        warning("Multiple matching clustering runs found. Returning the first one.", immediate. = TRUE)
        Stringendo::iprint(matchingRuns)
      }
      res <- matchingRuns[1]
    } else {
      warning("No matching clustering runs found. Returning the last one.", immediate. = TRUE)
      res <- clusteringRuns[length(clusteringRuns)]
    }
  }
  return(res)
}


# _________________________________________________________________________________________________
#' @title CleanDuplicateScorenames
#' @description Remove duplicate scorenames from obj@mata.data.
#' @param obj Seurat single cell object, Default: obj
#' @seealso
#'  \code{\link[CodeAndRoll2]{grepv}}
#'  \code{\link[Stringendo]{iprint}}
#' @export
#' @importFrom CodeAndRoll2 grepv
#' @importFrom Stringendo iprint
CleanDuplicateScorenames <- function(obj = obj) { # Helper. When AddGOScore(), a '1' is added to the end of the column name. It is hereby removed.
  obj <- combined.obj
  cn <- colnames(obj@meta.data)

  nonGO <- sort(CodeAndRoll2::grepv(x = cn, pattern = paste0("^Score.GO.[0-9]{7}.*"), invert = TRUE))
  clean <- CodeAndRoll2::grepv(x = cn, pattern = paste0("^Score.GO.[0-9]{7}$"))

  appended <- CodeAndRoll2::grepv(x = cn, pattern = paste0("^Score.GO.[0-9]{7}\\.[0-9]$"))
  fixed <- gsub(x = appended, pattern = paste0("\\.[0-9]$"), replacement = "")
  fixed.keep <- which(!(fixed %in% clean))
  uniqueGO <- c(clean, fixed.keep)
  obj@meta.data <- obj@meta.data[, c(nonGO, uniqueGO)]

  Stringendo::iprint(
    length(clean), "GO's are clean, ", length(appended), "GO's are suffixed by .1 etc, of which",
    length(fixed.keep), "GO's had no clean counterpart. All", length(uniqueGO), "scores, are cleaned, fixed and unique now."
  )
  print("Metadata column order re-organized alphabetically, and GO-scores at the end.")
  return(obj)
}




# _________________________________________________________________________________________________
#' @title IntersectWithExpressed
#' @description Intersect a gene set with the list of expressed genes.
#' @param genes Charcter vector of genes
#' @param obj Seurat single cell object, Default: combined.obj
#' @param genes.shown Number of genes shown, Default: 10
#' @seealso
#'  \code{\link[Stringendo]{iprint}}
#' @export
#' @importFrom Stringendo iprint

IntersectWithExpressed <- function(genes, obj = combined.obj, genes.shown = 10) { # Intersect a set of genes with genes in the Seurat object.
  print("IntersectWithExpressed()")
  # print(utils::head(genes, n=15))
  diff <- setdiff(genes, rownames(obj))
  Stringendo::iprint(length(diff), "genes (of", length(genes), ") are MISSING from the Seurat object:", utils::head(diff, genes.shown))
  return(intersect(rownames(obj), genes))
}



# _________________________________________________________________________________________________
#' @title .convert.GO_term.2.score
#' @description Convert a string GO_term-name to Score-name.
#' @param GO_term GO-term; Default: "GO:0006096"
#' @return Score name, e.g.: "Score.GO.0006096"
#'
.convert.GO_term.2.score <- function(GO_term = "GO:0006096") {
  if (grepl(x = GO_term, pattern = "^GO:", perl = TRUE)) {
    gsub(x = GO_term, pattern = ".*?GO:", replacement = "Score.GO.")
  } else {
    print("Input is not a GO-term, e.g.: GO:0006096")
    GO_term
  }
}



# _________________________________________________________________________________________________
#' @title .convert.score.2.GO_term
#' @description Convert a string Score-name to GO_term-name.
#' @param ScoreNames Vector of ScoreNames to convert to GO term, Default: "Score.GO.0006096"
#' @return GO term name, e.g.: "GO:0006096"
#'
.convert.score.2.GO_term <- function(ScoreNames = "Score.GO.0006096") {
  gsub(x = ScoreNames, pattern = "Score.GO.", replacement = "GO:")
}

# _________________________________________________________________________________________________
#' @title .fix.metad.colname.rm.trailing.1
#' @description Fix malformed (suffxed) column names in obj@meta.data. When running AddGOScore(),
#' a '1' is added to the end of the column name. It is hereby removed.
#' @param obj Seurat single cell object, Default: obj
#' @param colname Metadata column name to fix, Default: ScoreName
#' @seealso
#'  \code{\link[Stringendo]{iprint}}
#'
#' @importFrom Stringendo iprint
.fix.metad.colname.rm.trailing.1 <- function(obj = obj, colname = ScoreName) {
  colnames(obj@meta.data) <-
    gsub(
      x = colnames(obj@meta.data),
      pattern = paste0(colname, 1),
      replacement = colname
    )
  Stringendo::iprint("Trailing '1' in metadata column name is removed. Column name:", colname)
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Parse GO Term from Granule-score Name
#'
#' @description Splits a Granule-score column name into its clustering resolution and
#' GO term, and returns the 2nd part, the GO Term.
#' @param ident The column name to be split, e.g., 'RNA_snn_res.6.reassigned_cl.av_GO:0006096'.
#' @param pattern The pattern to split by, default '_cl\\.av_'.
#' @param ... Additional arguments passed to `strsplit`.
#'
#' @return GO Term as string, e.g. "GO:0006096"
#' @examples
#' .parse.GO(ident = "RNA_snn_res.6.reassigned_cl.av_GO:0006096")
.parse.GO <- function(ident = "RNA_snn_res.6.reassigned_cl.av_GO:0006096",
                      pattern = "_cl\\.av_", ...) {
  strsplit(x = ident, split = pattern, ...)[[1]][2]
}

# _____________________________________________________________________________________________ ----
# 8. Deprecated  ---------------------------------------------------------------------------

# These functions may have been used in the publication, but are not needed for the pipeline.


aut.res.clustering <- function() .Deprecated("gruffi::AutoFindGranuleResolution()")
reassign.small.clusters <- function() .Deprecated("gruffi::ReassignSmallClusters()")
GOscoreEvaluation <- function() .Deprecated("gruffi::AssignGranuleAverageScoresFromGOterm()")
GO_score_evaluation <- function() .Deprecated("gruffi::AssignGranuleAverageScoresFromGOterm()")

stand_dev_skewed <- function() .Deprecated("gruffi::CalcStandDevSkewedDistr()")
GetAllGOTerms <- function() .Deprecated("gruffi::GetAllGOTerms()")

calc.cluster.averages.gruff <- function() .Deprecated("gruffi::CalcClusterAverages_Gruffi()")
ww.convert.GO_term.2.score.Rd <- function() .Deprecated("gruffi:::.convert.GO_term.2.score.Rd()")
ww.convert.score.2.GO_term.Rd <- function() .Deprecated("gruffi:::.convert.score.2.GO_term.Rd()")
fix.metad.colname.rm.trailing.1 <- function() .Deprecated("gruffi:::.fix.metad.colname.rm.trailing.1()")
plot.clust.size.distr <- function() .Deprecated("gruffi::PlotClustSizeDistr()")
plot_norm_and_skew <- function() .Deprecated("gruffi::PlotNormAndSkew()")

