######################################################################################################
# Create_the_gruffi_Package.R
######################################################################################################
# source("~/GitHub/Packages/gruffiDev/dev/Create_the_gruffi_Package.R")
# rm(list = ls(all.names = TRUE));
try(dev.off(), silent = TRUE)

# Functions ------------------------
require("devtools")


# Setup ------------------------
PackageName <- "gruffi"
package.version <- "1.0.2"
setwd("~/GitHub/Packages/")

RepositoryDir <- kollapse("~/GitHub/Packages/", PackageName, "/")
fname <- 	kollapse(PackageName, ".R")
Package_FnP <- 	kollapse(RepositoryDir, "R/", fname)

BackupDir <- paste0("~/GitHub/Packages/",PackageName,"/dev/")
dir.create(BackupDir)


# devtools::use_package("vioplot")
DESCRIPTION <- list("Title" = "Gruffi identifies and removes stressed cells from brain organoid single-cell datasets"
                    , "Author" = c(
                      person(given = "Julia", family = "Naas", email = "julia.naas@meduniwien.ac.at", role =  c("aut", "cre") ),
                      person(given = "Abel", family = "Vertesy", email = "abel.vertesy@imba.oeaw.ac.at", role =  c("aut", "cre") )
                    )
                    , "Description" = "The Gruffi R package helps you (1) to identify stressed cells in single-cell
    RNA-seq datasets using  *granular funcitonal filtering*, and (2) you can use it to calculate any
    GO-term defined gene set's pathway activity. Gruffi integrates into single-cell analysis with
    Seurat and comes with a graphical user interface."
                    , "License" = "GPL-3 + file LICENSE"
                    , "Version" = package.version
                    , "Packaged" =  Sys.time()
                    , "Depends" =  "Seurat, magrittr, Stringendo, MarkdownReports"
                    , "Imports" = "cowplot, dplyr, ggplot2, raster, DOSE
                    , MarkdownHelpers (>= 1.0.1), CodeAndRoll2, Seurat.utils, ggExpress, stringr, sm, AnnotationDbi
                    , IRanges, Matrix, biomaRt, clipr, htmlwidgets, org.Hs.eg.db, rgl, rlang, shiny, tictoc, viridis"
                    # , "Suggests" = ""
                    , "BugReports"= "https://github.com/jn-goe/gruffi/issues/"
)



setwd(RepositoryDir)
if ( !dir.exists(RepositoryDir) ) { devtools::create(path = RepositoryDir, description = DESCRIPTION, rstudio = rstudioapi::isAvailable())
} else {
    getwd()
    try(file.remove(c("DESCRIPTION","NAMESPACE", "gruffi.Rproj")))
    usethis::create_package(path = RepositoryDir, fields = DESCRIPTION, open = F, check_name = F)
}



# Compile a package ------------------------------------------------
setwd(RepositoryDir)
getwd()
devtools::document()

# {
#   "update cff version"
#   citpath <- paste0(RepositoryDir, 'CITATION.cff')
#   xfun::gsub_file(file = citpath, perl = T
#                   , "^version: v.+", paste0("version: v", package.version))
# }

# Install your package ------------------------------------------------
devtools::install(RepositoryDir, upgrade = F)
# require("gruffi")
# remove.packages("gruffi")


# Test if you can install from github ------------------------------------------------
# devtools::install_github(repo = "jn-goe/gruffi")
# require("gruffi")


# Clean up if not needed anymore ------------------------------------------------
# View(installed.packages())
# remove.packages("gruffi")

check(RepositoryDir, cran = TRUE)
# as.package(RepositoryDir)
#
#
# # source("https://install-github.me/r-lib/desc")
# # library(desc)
# # desc$set("gruffi", "foo")
# # desc$get(gruffi)
#
#
# system("cd ~/GitHub/gruffi/; ls -a; open .Rbuildignore")
#
# Check package dependencies ------------------------------------------------
depFile <- paste0(RepositoryDir, 'dev/Dependencies.R')

(f.deps <- NCmisc::list.functions.in.file(filename = Package_FnP))
# clipr::write_clip(f.deps)

sink(file = depFile); print(f.deps); sink()
p.deps <- gsub(x = names(f.deps), pattern = 'package:', replacement = '')
write(x = p.deps, file = depFile, append = T)
p.dep.declared <- trimws(unlist(strsplit(DESCRIPTION$Imports, ",")))
p.dep.new <- sort(union( p.deps, p.dep.declared))
# clipr::write_clip(p.dep.new)


