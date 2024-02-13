## List of Functions in gruffi.R (29) 
Updated: 2024/02/13 18:53
- #### 1 `reassign.small.clusters()`
reassign.small.clusters. Reassign granules (clusters) smaller than X to nearby granules based on 3D UMAP coordinates.

- #### 2 `calculateMedianClusterSize()`
Calculate Median Cluster Size. This function calculates and returns the median size of clusters for a given assay  and resolution in a Seurat object. It provides an option to display the median size via message. 

- #### 3 `GetAllGOTerms()`
GetAllGOTerms. Return all GO-terms present in a Seurat single cell object.

- #### 4 `GetGOTerms()`
GetGOTerms. Get GO Terms

- #### 5 `AddGOGeneList.manual()`
AddGOGeneList.manual. Add a GO-term gene list under obj@misc$GO$xxxx.

- #### 6 `AddGOScore()`
AddGOScore. Add GO Score to a Seurat object

- #### 7 `AddCustomScore()`
AddCustomScore. Add a custom gene expression score from a set of genes provided to the  Seurat single cell object. Calculates Score for gene set. Fixes name.

- #### 8 `GO_score_evaluation()`
GO_score_evaluation. GO-score evaluation for filtering.

- #### 9 `CustomScoreEvaluation()`
CustomScoreEvaluation. Custom gene-set derived score evaluation for filtering.

- #### 10 `Shiny.GO.thresh()`
Shiny.GO.thresh. GO-thresholding for the Shiny app.

- #### 11 `FilterStressedCells()`
Filter Stressed Cells from a Seurat Object. Identifies and filters stressed cells based on specified GO terms and a quantile threshold.  It supports optional plotting of exclusion results and saving of modified datasets.

- #### 12 `PlotGoTermScores()`
PlotGoTermScores. Plot GO-term scores.

- #### 13 `FeaturePlotSaveCustomScore()`
FeaturePlotSaveCustomScore. Plot and save a Seurat FeaturePlot with a custom score in the meta data.

- #### 14 `FeaturePlotSaveGO()`
FeaturePlotSaveGO. Plot and save a Seurat FeaturePlot from a GO-score.

- #### 15 `plot.clust.size.distr()`
plot.clust.size.distr. Plot cluster size distribution.

- #### 16 `PlotNormAndSkew()`
PlotNormAndSkew. Plot normal distribution and skew.

- #### 17 `UMAP.3d.cubes()`
UMAP.3d.cubes. UMAP with 3D cubes.

- #### 18 `grScoreUMAP()`
Wrapper for clUMAP with Thresholding and Plotting. This function applies thresholding to granule scores within a specified object and column,  inverts the thresholding if a specific miscellaneous name is provided, and then calls clUMAP  for plotting based on the thresholded values. 

- #### 19 `grScoreHistogram()`
Histogram Visualization of Granule Scores with Threshold. Plots a histogram of granule scores from a specified column within a Seurat object,  with a vertical line denoting a threshold. The function allows for an optional  inversion of filtering logic based on the `miscname` and `auto` parameters. 

- #### 20 `clUMAP.thresholding()`
clUMAP.thresholding. clUMAP with threshold value.

- #### 21 `CalcTranscriptomePercentageGO()`
CalcTranscriptomePercentageGO. Calculate the percentage of transcriptome, for a given GO-term gene set already stored`in obj@misc$GO[[GO.score]].

- #### 22 `CalcTranscriptomePercentage()`
CalcTranscriptomePercentage. Calculate the percentage of transcriptome, for a given gene set.

- #### 23 `stand_dev_skewed()`
stand_dev_skewed. Standard deviation skew.

- #### 24 `calc.cluster.averages.gruffi()`
calc.cluster.averages.gruffi. Calculate granule (cluster) averages scores.

- #### 25 `CleanDuplicateScorenames()`
CleanDuplicateScorenames. Remove duplicate scorenames from obj@mata.data.

- #### 26 `IntersectWithExpressed()`
IntersectWithExpressed. Intersect a gene set with the list of expressed genes.

- #### 27 `PasteUniqueGeneList()`
PasteUniqueGeneList. Paste unique gene list  @seealso   \code{\link[clipr]{read_clip}}

- #### 28 `ww.convert.GO_term.2.score()`
ww.convert.GO_term.2.score. Convert a string GO_term-name to Score-name.

- #### 29 `ww.convert.score.2.GO_term()`
ww.convert.score.2.GO_term. Convert a string Score-name to GO_term-name.

