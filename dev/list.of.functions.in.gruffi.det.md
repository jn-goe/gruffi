## List of Functions in gruffi.R (32) 

Updated: 2024/02/26 11:15

- #### 1 `ReassignSmallClusters()`

  ReassignSmallClusters. Reassign granules (clusters) smaller than X to nearby granules based on 3D UMAP coordinates.

- #### 2 `CalculateMedianClusterSize()`

  Calculate Median Cluster Size. This function calculates and returns the median size of clusters for a given assay  and resolution in a Seurat object. It provides an option to display the median size via message. 

- #### 3 `GetGOTerms()`

  Retrieve Gene Ontology (GO) Terms and Associated Genes. Fetches genes associated with a specified Gene Ontology (GO) term and optionally  opens the GO term page. The genes are retrieved using either the Ensembl database via biomaRt or   from precomputed GO_genes results stored in the Seurat object. 

- #### 4 `GetAllGOTermNames()`

  Show All GO Terms in an Object, and retrieve corresponding names using AnnotationDbi. Retrieves and optionally annotates all GO term identifiers present in the  metadata of a Seurat object. It identifies columns starting with 'Score.GO',  extracts the corresponding GO terms, retrieves tge corresponding term-names using AnnotationDbi,  and optionally saving them back into the object.

- #### 5 `CalculateAndPlotGoTermScores()`

  Calculate and Plot GO Term Scores. Automates the process of obtaining and calculating of Gene Ontology (GO) term-based    gene set scores and it's visualization through UMAP plots. 

- #### 6 `AssignGranuleAverageScoresFromGOterm()`

  Assign granule average scores from Gene Ontology (GO) terms.. Calculates and assigns granule average scores from Gene Ontology (GO) term in a  Seurat object. It can compute new GO term scores, plot gene expressions, and save UMAP visualizations.  The function updates the Seurat object with cluster average scores for the specified GO term. 

- #### 7 `AddGOGeneList.manual()`

  AddGOGeneList.manual. Add a GO-term gene list under obj@misc$gruffi$GO$xxxx.

- #### 8 `AddGOScore()`

  AddGOScore. Add GO Score to a Seurat object

- #### 9 `AddCustomScore()`

  AddCustomScore. Add a custom gene expression score from a set of genes provided to the  Seurat single cell object. Calculates Score for gene set. Fixes name.

- #### 10 `CustomScoreEvaluation()`

  CustomScoreEvaluation. Custom gene-set derived score evaluation for filtering.

- #### 11 `FindThresholdsShiny()`

  Launch Shiny App for GO Term-based Thresholding. Launches a Shiny application to interactively apply GO term-based thresholding  for identifying stressed cells in a Seurat object. This function calculates the  "granule average GO-scores" for individual GO terms and updates the object with a combined stress  identification based on these thresholds.#'

- #### 12 `FindThresholdsAuto()`

  Automated GO Term-based Thresholding without using a Shiny app.. Automatically applies granular GO term-based thresholding to identify stressed cells  in a Seurat object without launching a Shiny app. This function calculates the  "granule average GO-scores" for individual GO terms and updates the object with a combined stress  identification based on these thresholds.

- #### 13 `FeaturePlotSaveCustomScore()`

  FeaturePlotSaveCustomScore. Plot and save a Seurat FeaturePlot with a custom score in the meta data.

- #### 14 `FeaturePlotSaveGO()`

  FeaturePlotSaveGO. Plot and save a Seurat FeaturePlot from a GO-score.

- #### 15 `PlotClustSizeDistr()`

  PlotClustSizeDistr. Plot cluster size distribution.

- #### 16 `PlotNormAndSkew()`

  PlotNormAndSkew. Plot normal distribution and skew.

- #### 17 `GrScoreUMAP()`

  Wrapper for clUMAP with Thresholding and Plotting. This function applies thresholding to granule scores within a specified object and column,  inverts the thresholding if a specific miscellaneous name is provided, and then calls clUMAP  for plotting based on the thresholded values. 

- #### 18 `GrScoreHistogram()`

  Histogram Visualization of Granule Scores with Threshold. Plots a histogram of granule scores from a specified column within a Seurat object,  with a vertical line denoting a threshold. The function allows for an optional  inversion of filtering logic based on the `miscname` and `auto` parameters. 

- #### 19 `ClusterUMAPthresholding()`

  ClusterUMAPthresholding. clUMAP with threshold value.

- #### 20 `CalcTranscriptomePercentageGO()`

  CalcTranscriptomePercentageGO. Calculate the percentage of transcriptome, for a given GO-term gene set  already stored`in obj@misc$gruffi$GO[[GO.score]].

- #### 21 `CalcTranscriptomePercentage()`

  CalcTranscriptomePercentage. Calculate the percentage of transcriptome, for a given gene set.

- #### 22 `CalcStandDevSkewedDistr()`

  Calculate Standard Deviation for Skewed Distributions. Computes a modified standard deviation for skewed distributions by mirroring data points about the mean.  This approach aims to reflect the skewed part of the distribution to better understand its spread. 

- #### 23 `CalcClusterAverages_Gruffi()`

  Calculate Cluster Averages for Granules in a Seurat Object. Calculates the central tendency (mean or median) of specified scores for clusters (granules)  within a Seurat object and optionally plots these averages. 

- #### 24 `GetGruffiClusteringName()`

  Retrieve First Clustering Run or First Matching Pattern Run. Fetches the first clustering run from a Seurat object,  optionally filtered by a specific pattern (e.g., ".reassigned"). 

- #### 25 `ParseGruffiGranuleScoreName()`

  Generate Gruffi Granule Score Name. This function creates a standardized score name using granule resolution and a  Gene Ontology (GO) identifier.

- #### 26 `CleanDuplicateScorenames()`

  CleanDuplicateScorenames. Remove duplicate scorenames from obj@mata.data.

- #### 27 `IntersectWithExpressed()`

  IntersectWithExpressed. Intersect a gene set with the list of expressed genes.

- #### 28 `ClearGruffi()`

  Remove gruffi results from the Seurat object. Removes corresponding metadata columns and clear @misc$gruffi slot 

- #### 29 `.convert.GO_term.2.score()`

  .convert.GO_term.2.score. Convert a string GO_term-name to Score-name.

- #### 30 `.convert.score.2.GO_term()`

  .convert.score.2.GO_term. Convert a string Score-name to GO_term-name.

- #### 31 `.fix.metad.colname.rm.trailing.1()`

  .fix.metad.colname.rm.trailing.1. Fix malformed (suffxed) column names in obj@meta.data. When running AddGOScore(),  a '1' is added to the end of the column name. It is hereby removed.

- #### 32 `.parse.GO()`

  Parse GO Term from Granule-score Name. Splits a Granule-score column name into its clustering resolution and  GO term, and returns the 2nd part, the GO Term.

