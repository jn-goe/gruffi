# GRUFFI

<img width="100" alt="gruffi" src="https://user-images.githubusercontent.com/72695751/131385126-76993ea7-3746-4681-881f-56997ba28055.png"> 
© Walt Disney

The Gruffi R package helps you (1) to identify stressed cells in single-cell RNA-seq datasets using  *granular funcitonal filtering*, or (2) you can use it to analyze any gene set's pathway activity (GO-term or custom defined) , and define a cell population based on their combined activity.

`Gruffi` integrates into any `Seurat` analysis pipelione & it comes with a graphical user interface.



## Installation

Once the paper is published, we will provide a standalone version, for now, you can install dependencies from **CRAN, Bioconductor and GitHub** via **devtools**:

```R


# Install CRAN & Bioconductor dependencies
install.packages('pdist')
install.packages('raster')
install.packages('rgl')

install.packages('BiocManager')
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

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

### Currently you have to load the package, not simply `require()`

```r
devtools::load_all("~/example/path/to/gruffi")
```

Otherwise the Shiny scripts will not be found. *You also can skip installing **Gruffi**,  but then function help won't work.*



## Usage

### 1. <u>Setup</u>

Load your Seurat single-cell object of interest (in the following `combined.obj`), perform basic analysis steps (Integration, PCA, Clustering) and calculate 2D and 3D umaps.

```R
combined.obj <- Seurat.utils::SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = 50, dimensions=3:2, reduction="umap")
# Note that this will recalculate 2D and 3D umaps, and back them up in combined.obj@misc$reductions.backup. 
# See FAQ.md, if you don't want to overwrite the 2D UMAP.
```

If you want to store multiple UMAP's, you have to keep them in backup slot (in combined.obj@misc). Use `combined.obj <- RecallReduction(obj = combined.obj, dim = 2, reduction="umap")` to access them.

Prepare GO-terms and gene-sets

```r
ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

go1 <- "GO:0006096" # Glycolysis
go2 <- "GO:0034976" # ER-stress
go3 <- "GO:0042063" # Gliogenesis, negative filtering

```

### 2. <u>Granule partitioning</u>

Gruffi works best if you partition cells into groups of 100-200 cells. Find the corresponding clustering resolution by:

```  R
combined.obj <- aut.res.clustering(obj = combined.obj)
```

The optimal resolution found is stored in:

```r
granule.res.4.gruffi <- combined.obj@misc$gruffi$'optimal.granule.res'	
```

Some granules have too few cells, therfore their scoring is not robust statistically. Use `reassign.small.clusters` to assign these cells to the nearest large-enough granule:

```R
combined.obj <- reassign.small.clusters(combined.obj, ident = granule.res.4.gruffi) # will be stored in meta data column as "seurat_clusters.reassigned"
````

Above, granules with <30 cells are cell-by-cell re-assigned to a neighboring granule (by default based on Euclidean distance between the mean of cell groups in 3dim UMAP space). The reassigned granules are suffixed as :

```r
granule.res.4.gruffi <- paste0(granule.res.4.gruffi, '.reassigned')
```



### 3. <u>Pathway scoring</u>

After finding the right granule resolution, first GO scores per cells, then average GO-scores per granule (normalized for respective cell numbers) are computed. Within the GO score computation an ensembl database will be loaded (that you created above as `ensembl`).

```R
# Glycolytic process	GO:0006096
combined.obj <- GO_score_evaluation(obj = combined.obj, GO_term = go1, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi)

# ER stress 	GO:0034976
combined.obj <- GO_score_evaluation(obj = combined.obj, GO_term = go2, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi)

# Gliogenesis		GO:0042063
combined.obj <- GO_score_evaluation(obj = combined.obj, GO_term = go3, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi)
```

These functions store the resulting scores in `combined.obj@meta.data`.



### 4. <u>Stress filtering</u>

We will now call a Shiny Interface to auto-estimate and/or manually adjust the stress assignment threeholds via `Shiny.GO.thresh()`.

- The first 2 scores are for features to filter out (e.g.: cells with high ER stress and glycolysis), 
- and the last 2 scores are for features we want to keep  (e.g.: cells with high gliogenesis). 

Example code for filtering cells high in glycolytic process and ER stress but low in gliogenesis: 
```R
# Create score names:
(i1 <- Stringendo::kppu(granule.res.4.gruffi, go1))
(i2 <- Stringendo::kppu(granule.res.4.gruffi, go2))
(i3 <- Stringendo::kppu(granule.res.4.gruffi, go3))

# Call Shiny app
combined.obj <- Shiny.GO.thresh(obj = combined.obj,
                                stress.ident1 = i1,
                                stress.ident2 = i2,
                                notstress.ident3 = i3,
                                plot.cluster.shiny = "orig.ident")

"Dont forget to click the button in the app: Save New Thresholds"

```
The interfacte allows automatic estimation and manual adjustment of thresholds for stress scores: 

<img width="396" alt="ShinyInterface_Vel Org7 d90 ImpV" src="https://user-images.githubusercontent.com/72695751/156938645-d80fd987-36b1-46dd-b350-4fc26f62035e.png">



After pushing the **<u>Save new thresholds</u>** button in the Shiny graphical user interface, thresholds are saved in `combined.obj@misc$gruffi` and the stress assignment is stored as a new meta data column `is.Stressed`. Check results as:

```r
Seurat.utils::clUMAP('is.Stressed', label =F)
```

### 5. <u>Remove stressed cells</u>

```R
cellIDs.keep <- which_names(!combined.obj$'is.Stressed')
subset.obj <- subset(x = combined.obj, cells = cellIDs.keep)  

Seurat.utils::clUMAP('is.Stressed', label = F, obj = subset.obj)
```





## Troubleshooting

For Errors, Problems & Questions: **please [raise an issue in the github](https://github.com/jn-goe/gruffi/issues/new) repository**. We try to answer.



