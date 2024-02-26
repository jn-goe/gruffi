# GRUFFI

#### _Published in [Gruffi: an algorithm for computational removal of stressed cells from brain organoid transcriptomic datasets](https://www.embopress.org/doi/full/10.15252/embj.2022111118)_ 

![status: active](https://raw.githubusercontent.com/vertesy/TheCorvinas/master/GitHub/Badges/active.svg)

The Gruffi R package helps you (1) to identify stressed cells in single-cell RNA-seq datasets using  *granular funcitonal filtering*, or (2) you can use it to analyze any gene set's pathway activity (GO-term or custom defined), and define a cell population based on their combined activity.

`Gruffi` integrates into any `Seurat` analysis pipelione & it comes with a graphical user interface.




## News

- `v1.5.*` is released, that fixes critical problems which arose with changes in dependencies and made the installation / pipeline break.  Additionally it contains 
- Now it is more explicit that GO-annotation can either be obtained via `BioMart` or `AnnotationDbi`, which is helpful, as BioMart connection is not always working. For the same reason, it now allows the usage of alternative mirrors for BioMart (thanks to `@zuzkamat`).
- Major consistency update in function names, see below.
- Consistency update in variable names (in @misc and @meta.data), for the latter, this update is not backward compatible, meaning earlier results needs manual adjustments in those names to run. Recommended is to rerun the new version.
- Major cleanup and of the codebase. Granule average scores are now stored as numeric rather than factor.

<u>*Added  functionality*</u>

- `FindThresholdsAuto()` that is a shiny free version  `FindThresholdsShiny()`, which allows an automated workflow without the need for graphical output (thanks to `@heyfl`).  We still strongly recommend the shiny based workflow that enforces users to look at the results and intermediate steps instead of blidnly taking a TRUE / FALSE column.
- New visualization functions, incuding:
  ```r
  GrScoreUMAP(obj = combined.obj, colname = 'RNA_snn_res.6_cl.av_GO.0042063', miscname = 'thresh.stress.ident1')
  GrScoreHistogram(obj = combined.obj, colname = 'RNA_snn_res.6_cl.av_GO.0042063', miscname = 'thresh.stress.ident1')
  ```
- New helper functions that allow a simpler workflow and shorter codebase.




## Installation

You can install dependencies from **CRAN, Bioconductor and GitHub** via **devtools**:

```R

# Install CRAN & Bioconductor dependencies
install.packages('BiocManager')
BiocManager::install("DOSE")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("sparseMatrixStats")
BiocManager::install("biomaRt")


# Install custom dependencies
install.packages('devtools')
devtools::install_github(repo = "vertesy/Stringendo", upgrade = F)
devtools::install_github(repo = "vertesy/CodeAndRoll2", upgrade = F)
devtools::install_github(repo = "vertesy/ReadWriter", upgrade = F)
devtools::install_github(repo = "vertesy/MarkdownHelpers", upgrade = F)
devtools::install_github(repo = "vertesy/Markdownreports", upgrade = F)
devtools::install_github(repo = "vertesy/ggExpress", upgrade = F)
devtools::install_github(repo = "vertesy/Seurat.utils", upgrade = F)


# Install gruffi
devtools::install_github(repo = "jn-goe/gruffi", upgrade = F)

```

*NOTE: If you type 'all' when R asks to update dependencies, you may get into installation errors / infinite loops. If updating fails, type 'no' when prompted.*


## Usage
<img width="1000" alt="image" src="https://user-images.githubusercontent.com/5101911/170080428-7cacdc2b-c079-4383-a2ec-c9bb3e6271df.png">


### 1. <u>Setup</u>

Load your Seurat single-cell object of interest (in the following `combined.obj`), perform basic analysis steps (Integration, PCA, Clustering) and calculate 2D and 3D umaps.

```R
library(gruffi)

combined.obj <- Seurat.utils::SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = 50, dimensions=3:2, reduction="umap")
# Note that this will recalculate 2D and 3D umaps, and back them up in combined.obj@misc$reductions.backup. 
# See FAQ.md, if you don't want to overwrite the 2D UMAP.
```

If you want to store multiple UMAP's, you have to keep them in backup slot (in combined.obj@misc). Use `combined.obj <- RecallReduction(obj = combined.obj, dim = 2, reduction="umap")` to access them.

Prepare GO-terms and gene-sets

```r
go1 <- "GO:0006096" # Glycolysis
go2 <- "GO:0034976" # ER-stress
go3 <- "GO:0042063" # Gliogenesis, negative filtering
```

### 2. <u>Granule partitioning</u>

Gruffi works best if you partition cells into groups of 100-200 cells. Find the corresponding clustering resolution by:

```  R
combined.obj <- AutoFindGranuleResolution(obj = combined.obj)
```

The optimal resolution found is stored in:

```r
(granule.res.4.gruffi <- GetGruffiClusteringName(combined.obj)) # Recalled from @misc$gruffi$'optimal.granule.res.
```

Some granules have too few cells, therfore their scoring is not robust statistically. Use `ReassignSmallClusters` to assign these cells to the nearest large-enough granule:

```R
combined.obj <- ReassignSmallClusters(combined.obj, ident = granule.res.4.gruffi) # Will be stored in meta data column as "seurat_clusters.reassigned".
````

Above, granules with <30 cells are cell-by-cell re-assigned to a neighboring granule (by default based on Euclidean distance between the mean of cell groups in 3dim UMAP space). The reassigned granules are suffixed as :

```r
(granule.res.4.gruffi <- GetGruffiClusteringName(combined.obj))
```



### 3. <u>Pathway scoring</u>

After finding the right granule resolution, first GO scores per cells, then average GO-scores per granule (normalized for respective cell numbers) are computed. Within the GO score computation an ensembl database will be loaded (that you created above as `ensembl`).

```R
# Glycolytic process	GO:0006096
combined.obj <- AssignGranuleAverageScoresFromGOterm(obj = combined.obj, GO_term = go1, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)

# ER stress 	GO:0034976
combined.obj <- AssignGranuleAverageScoresFromGOterm(obj = combined.obj, GO_term = go2, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)

# Gliogenesis		GO:0042063
combined.obj <- AssignGranuleAverageScoresFromGOterm(obj = combined.obj, GO_term = go3, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)
```

These functions store the resulting scores in `combined.obj@meta.data`.



### 4. <u>Stress filtering</u>

We will now call a Shiny Interface to auto-estimate and/or manually adjust the stress assignment threeholds via `FindThresholdsShiny()`.

- The first 2 scores are for features to filter out (e.g.: cells with high ER stress and glycolysis), 
- and the last 2 scores are for features we want to keep  (e.g.: cells with high gliogenesis). 

Example code for filtering cells high in glycolytic process and ER stress but low in gliogenesis: 
```R
# Create score names:
(i1 <- ParseGruffiGranuleScoreName(goID = go1))
(i2 <- ParseGruffiGranuleScoreName(goID = go3))
(i3 <- ParseGruffiGranuleScoreName(goID = go3))

# Call Shiny app
combined.obj <- FindThresholdsShiny(obj = combined.obj,
                                stress.ident1 = i1,
                                stress.ident2 = i2,
                                notstress.ident3 = i3,
                                plot.cluster.shiny = "orig.ident")

"Dont forget to click the button in the app: Save New Thresholds"

```
The interface allows automatic estimation and manual adjustment of thresholds for stress scores: 

<img width="396" alt="ShinyInterface_Vel Org7 d90 ImpV" src="https://user-images.githubusercontent.com/72695751/156938645-d80fd987-36b1-46dd-b350-4fc26f62035e.png">



After pushing the **<u>Save new thresholds</u>** button in the Shiny graphical user interface, thresholds are saved in `combined.obj@misc$gruffi` and the stress assignment is stored as a new `@meta.data` column, `is.Stressed`. Check results with:

### 5. <u>Visualize results</u>
```r
StressUMAP(combined.obj)
StressBarplotPerCluster(combined.obj,  group.by = GetClusteringRuns()[1])
GrScoreHistogram(combined.obj, colname = i1, miscname = "thresh.stress.ident1")
# and more
```

### 6. <u>Remove stressed cells</u>

```R
cellIDs.keep <- which_names(!combined.obj$'is.Stressed')
subset.obj <- subset(x = combined.obj, cells = cellIDs.keep)  

StressUMAP(combined.obj)
```



## Troubleshooting

For Errors, Problems & Questions: **please [see the FAQ](https://github.com/jn-goe/gruffi/blob/main/FAQ.md)**
**or [raise an issue in the github](https://github.com/jn-goe/gruffi/issues/new) repository**. We try to answer.

