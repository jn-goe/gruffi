# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)


#### 1. Setup ------------------------------------------------------------
require(gruffi)
devtools::load_all(path = "~/GitHub/Packages/gruffi/")

go1 <- "GO:0006096" # Glyco
go2 <- "GO:0034976" # ER
go3 <- "GO:0042063" # Glia

ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

MarkdownReports::create_set_Original_OutDir()


combined.obj <- readRDS("obj.cluster.RDS")


#### 2. Granule partitioning ------------------------------------------------------------

combined.obj <- AutoFindGranuleResolution(combined.obj) # in the following "integrated_snn_res.25.reassigned"
granule.res.4.gruffi <- combined.obj@misc$gruffi$'optimal.granule.res'


combined.obj <- ReassignSmallClusters(combined.obj, ident = granule.res.4.gruffi)
granule.res.4.gruffi <- paste0(granule.res.4.gruffi, '.reassigned')



#### 3. Pathway scoring ------------------------------------------------------------

# Glycolytic process	GO:0006096
combined.obj <- GOscoreEvaluation(obj = combined.obj, GO_term = go1, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi)

# ER stress 	GO:0034976
combined.obj <- GOscoreEvaluation(obj = combined.obj, GO_term = go2, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi)

# Gliogenesis		GO:0042063
combined.obj <- GOscoreEvaluation(obj = combined.obj, GO_term = go3, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi)



#### 4. Stress filtering<  ------------------------------------------------------------

i1 <- kppu(granule.res.4.gruffi, 'cl.av', go1)
i2 <- kppu(granule.res.4.gruffi, 'cl.av', go2)
i3 <- kppu(granule.res.4.gruffi, 'cl.av', go3)

combined.obj <- FindThresholdsShiny(obj = combined.obj,
                                stress.ident1 = i1,
                                stress.ident2 = i2,
                                notstress.ident3 = i3,
                                plot.cluster.shiny = "orig.ident")

Seurat.utils::clUMAP('is.Stressed', label =F)


#### 5. Remove stressed cells  ------------------------------------------------------------

cellIDs.keep <- which_names(!combined.obj$'is.Stressed')
subset.obj <- subset(x = combined.obj, cells = cellIDs.keep)

Seurat.utils::clUMAP('is.Stressed', label = F, obj = subset.obj)


# ######## PLOTS ########
# # stress assignment barplot cell numbers
# ggplot2::ggplot(combined.obj@meta.data, ggplot2::aes(x=combined.obj$orig.ident, fill=is.Stressed, stat = "count")) +
#   ggplot2::scale_fill_manual(values=c(gg_color_hue(2)[2], gg_color_hue(2)[1])) + ggplot2::theme_minimal() +
#   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5, hjust=1)) +
#   ggplot2::geom_bar() + ggplot2::xlab("") +
#   ggplot2::ggsave(paste0("Barplot_stressed.cells.png"), width = 3, height = 7)
#
# # stress assignment UMAP
# clUMAP(obj = combined.obj, ident = "is.Stressed", save.plot = F, label = F, legend = F)+#, splitby = "orig.ident") +
#   ggplot2::theme(text = ggplot2::element_text(size=20)) + Seurat::NoAxes() + ggplot2::ggtitle(ggplot2::element_blank())+ ggplot2::scale_color_manual(values=c(gg_color_hue(2)[2], gg_color_hue(2)[1])) +
#   ggplot2::ggsave(paste0("UMAP_Stress_cluster.png"), width = 10, height = 7)
