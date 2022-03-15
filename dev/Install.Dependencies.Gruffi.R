# install.dependencies.gruffi.R


# Install standard dependencies ------------------------------------------------
install.packages('pdist')
install.packages('raster')
install.packages('rgl')

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

# Install standard dependencies ------------------------------------------------
devtools::install_github(repo = "vertesy/Stringendo", upgrade = F)
devtools::install_github(repo = "vertesy/ReadWriter", upgrade = F)
devtools::install_github(repo = "vertesy/CodeAndRoll2", upgrade = F)
devtools::install_github(repo = "vertesy/MarkdownHelpers", upgrade = F)
devtools::install_github(repo = "vertesy/Markdownreports", upgrade = F)
devtools::install_github(repo = "vertesy/ggExpress", upgrade = F)
devtools::install_github(repo = "vertesy/Seurat.utils", upgrade = F)


"Now try "
devtools::load_all(path = "~/GitHub/Packages/gruffiDev/")

require(gruffiDev)
