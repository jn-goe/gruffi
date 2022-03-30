######################################################################################################
# Create_the_gruffi_Package.R
######################################################################################################
# source("~/GitHub/Packages/gruffiDev/dev/Create_the_gruffi_Package.R")
# rm(list = ls(all.names = TRUE));
try(dev.off(), silent = TRUE)

# Functions ------------------------
# install_version("devtools", version = "2.0.2", repos = "http://cran.at.r-project.org") # install.packages("devtools")
require("devtools")
require("roxygen2")
require("stringr")

# devtools::install_github(repo = "vertesy/CodeAndRoll2")
require('CodeAndRoll2')
require('Stringendo')


# Setup ------------------------
PackageName = "gruffi"
package.version = "0.7.0"
setwd("~/GitHub/Packages/")

RepositoryDir = kollapse("~/GitHub/Packages/", PackageName, "/")
fname = 	kollapse(PackageName, ".R")
Package_FnP = 	kollapse(RepositoryDir, "R/", fname)

BackupDir = paste0("~/GitHub/Packages/",PackageName,"/dev/")
dir.create(BackupDir)


# devtools::use_package("vioplot")
DESCRIPTION <- list("Title" = "gruffi is the fastest way to create, annotate and export plots in R"
    , "Author" = c(
      person(given = "Julia", family = "Naas", email = "julia.naas@meduniwien.ac.at", role =  c("aut", "cre") ),
      person(given = "Abel", family = "Vertesy", email = "abel.vertesy@imba.oeaw.ac.at", role =  c("aut", "cre") )
      )
    # , "Authors@R" = 'person(given = "Julia", family = "Naas", email = "julia.naas@meduniwien.ac.at", role =  c("aut", "cre") ), person(given = "Abel", family = "Vertesy", email = "a.vertesy@imba.oeaw.ac.at", role =  c("aut", "cre") )'
    , "Description" = "The Gruffi R package helps you (1) to identify stressed cells in single-cell
    RNA-seq datasets using  *granular funcitonal filtering*, and (2) you can use it to calculate any
    GO-term defined gene set's pathway activity. Gruffi integrates into single-cell analysis with
    Seurat and comes with a graphical user interface."
    , "License" = "GPL-3 + file LICENSE"
    , "Version" = package.version
    , "Packaged" =  Sys.time()
    , "Depends" =  "Stringendo, magrittr"
    , "Imports" = "tidyverse, cowplot, ggpubr, graphics, grDevices, raster
      , MarkdownHelpers, MarkdownReports, CodeAndRoll2, Seurat.utils, ggExpress
      , methods, RColorBrewer, sessioninfo, Seurat, sm, stats, AnnotationDbi, IRanges
      , Matrix,biomaRt, clipr, htmlwidgets, org.Hs.eg.db, rgl, rlang, shiny, tictoc, viridis"
    # , 'dplyr'
    # , "Suggests" = ""
    , "BugReports"= "https://github.com/jn-goe/gruffi/issues/"
)


setwd(RepositoryDir)
if ( !dir.exists(RepositoryDir) ) { create(path = RepositoryDir, description = DESCRIPTION, rstudio = TRUE)
} else {
    getwd()
    try(file.remove(c("DESCRIPTION","NAMESPACE", "gruffi.Rproj")))
    create_package(path = RepositoryDir, fields = DESCRIPTION, open = F, check_name = F)
}



# file.copy(from = AnnotatedFile, to = Package_FnP, overwrite = TRUE)


# Compile a package ------------------------------------------------
setwd(RepositoryDir)
getwd()
document()

# Install your package ------------------------------------------------
# # setwd(RepositoryDir)
install(RepositoryDir, upgrade = F)

# require("gruffi")
# # remove.packages("gruffi")
# # Test your package ------------------------------------------------
# help("wplot")
# cat("\014")
# devtools::run_examples()

{
  "update cff version"
  citpath <- paste0(RepositoryDir, 'CITATION.cff')
  xfun::gsub_file(file = citpath, perl = T
                  , "^version: v.+", paste0("version: v", package.version))
}

# Test if you can install from github ------------------------------------------------
# devtools::install_github(repo = "vertesy/gruffi")

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
depFile = paste0(RepositoryDir, 'dev/Dependencies.R')

(f.deps <- NCmisc::list.functions.in.file(filename = Package_FnP))
# clipr::write_clip(f.deps)

sink(file = depFile); print(f.deps); sink()
p.deps <- gsub(x = names(f.deps), pattern = 'package:', replacement = '')
write(x = p.deps, file = depFile, append = T)
p.dep.declared <- trimws(unlist(strsplit(DESCRIPTION$Imports, ",")))
p.dep.new <- sort(union( p.deps, p.dep.declared))
# clipr::write_clip(p.dep.new)


