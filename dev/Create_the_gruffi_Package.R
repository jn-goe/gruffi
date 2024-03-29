######################################################################################################
# Create_the_gruffi_Package.R
######################################################################################################
# source("~/GitHub/Packages/gruffi/dev/Create_the_gruffi_Package.R")
rm(list = ls(all.names = TRUE));
try(dev.off(), silent = TRUE)

# Functions ------------------------
# require(PackageTools)
devtools::load_all("~/GitHub/Packages/PackageTools")
# devtools::load_all("~/GitHub/Packages/gruffi/")

# Setup ------------------------
repository.dir <- "~/GitHub/Packages/gruffi/"
config.path <- file.path(repository.dir, "dev/config.R")

"TAKE A LOOK AT"
file.edit(config.path)
source(config.path)
package.name <- DESCRIPTION$'package.name'

PackageTools::document_and_create_package(repository.dir, config_file = 'config.R', dev = "dev")
'git add commit push to remote'

# Install your package ------------------------------------------------
"disable rprofile by"
rprofile()
devtools::install_local(repository.dir, upgrade = F, force = T)


# Test if you can install from github ------------------------------------------------
remote.path <- file.path(DESCRIPTION$'github.user', package.name)
# pak::pkg_install(remote.path)

devtools::install_github(repo = "jn-goe/gruffi", upgrade = F)

# unload(package.name)
# require(package.name, character.only = TRUE)
# # remove.packages(package.name)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.18")
# BiocManager::install("MatrixGenerics")

# CMD CHECK ------------------------------------------------
checkres <- devtools::check(repository.dir, cran = FALSE)



# Automated Codebase linting to tidyverse style ------------------------------------------------
styler::style_pkg(repository.dir)
styler::style_file(paste0(repository.dir, "inst/shiny/GO.thresh/ui.R"))
styler::style_file(paste0(repository.dir, "inst/shiny/GO.thresh/server.R"))

# Extract package dependencies ------------------------------------------------
PackageTools::extract_package_dependencies(repository.dir)


# Visualize function dependencies within the package------------------------------------------------
{
  warning("works only on the installed version of the package!")
  pkgnet_result <- pkgnet::CreatePackageReport(package.name)
  fun_graph     <- pkgnet_result$FunctionReporter$pkg_graph$'igraph'
  PackageTools::convert_igraph_to_mermaid(graph = fun_graph, openMermaid = T, copy_to_clipboard = T)
}


# Try to find and add missing @importFrom statements------------------------------------------------
devtools::load_all("~/GitHub/Packages/PackageTools/")
(ls.scripts.full.path <- list.files(file.path(repository.dir, "R"), full.names = T))
if (F) {
  (excluded.packages <- unlist(strsplit(DESCRIPTION$'Depends', split = ", ")))
  for (scriptX in ls.scripts.full.path) {
    PackageTools::add_importFrom_statements(scriptX, exclude_packages = excluded.packages)
  }
}


# Generate the list of functions ------------------------------------------------
for (scriptX in ls.scripts.full.path) {
  PackageTools::list_of_funs_to_markdown(scriptX)
}

PackageTools::copy_github_badge("active") # Add badge to readme via clipboard


# Replaces T with TRUE and F with FALSE ------------------------------------------------
for (scriptX in ls.scripts.full.path) {
  PackageTools::replace_short_calls(scriptX, strict_mode = F)
  PackageTools::replace_tf_with_true_false(scriptX, strict_mode = F)
}


