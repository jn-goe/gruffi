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
