# FAQ -  Frequently Asked Questions



## 1. The 2D UMAP looks different after running Gruffi. How can I keep the original UMAP?

The `SetupReductionsNtoKdimensions()` function simply calls the `RunUMAP()` function, and backs it up into `obj@misc$reductions.backup` , so if you provide the exact same parameters, it should reproduce the UMAP you generated earlier.

An alternatively you can 

```R
# 1. backup your old 2D umap
combined.obj@misc$reductions.backup$umap2d <- combined.obj@reductions$umap


# 2. calculate and backup the new 3D umap | "dimensions=3"
combined.obj <- Seurat.utils::SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = 50, dimensions=3, reduction="umap")

# 3. Recall the old 2D umap from @misc
combined.obj@reductions$umap <- combined.obj@misc$reductions.backup$umap2d
# The 3D is in ...reductions.backup$umap3d


# 4. Gruffi will find the 3D umap in @misc - it is needed for the reclassification step.
```




## 2. How to use Gruffi with different integrations?

Gruffi assumes that you used Seurat, thus granule clustering runs on `assays@RNA` or  `assays@integrated`. 

Typically, integration methods simply provide an alternative reduction to `PCA` such as  `iNMF` for LIGER batch correction. Simply provide this reduction when calculateing the 3D UMAP (see 2.).

If you integrated with a method that created a different assay (e.g.: `@BlaBla`), provide the assay name in `AutoFindGranuleResolution( , assay = "BlaBla")`.




## 3. How to use Gruffi when providing a custom set of genes, instead of a GO-term?

Gruffi can work with any set of genes to classify cells. In the current implementation you can provide 2 positive and 2 negative filtering terms (in the default implementation we use 2 positive and 1 negative terms).
```r
# You add scores of your gene set and calculate granule averages:
heGENES <- c("TMSB4X", "NRXN3", "SNTG1","SOX4", "TUBA1A", "NRXN1", "TMSB10", "ACTG1", "ROBO2","ACTB")
combined <- AddCustomScore(obj = combined.obj, genes = heGENES)
# Scorename: 'Score.heGENES'

# Then calculate granule averages
combined.obj <- CustomScoreEvaluation(obj = combined.obj, custom.score.name = 'Score.heGENES')
```

Finally continue with the Gruffi pipeline as normal, by calling the Shiny app. This extension of the Gruffi workflow is less tested, so if you find bugs / encounter an error, [please let us know](https://github.com/jn-goe/gruffi/issues).


## 4. Does it work with Seurat v5+?

Yes, it should. Please  raise an `issue` if you experience problems.

## 5. Does it install and run work with R v4.4.X?

Yes, it should. You may need to manually install `terra` and `rgl` dependencies.
Please  raise an `issue` if you experience problems.
