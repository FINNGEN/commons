#!/usr/bin/env Rscript

# R 4.2.3 bundles Matrix 1.5-3 but several packages require >= 1.6.0.
# Current CRAN Matrix requires R > 4.2.3, so install last compatible version from archive.
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz",
  repos = NULL, type = "source"
)

cran_packages <- c(
"data.table",
"tidyverse",
"rjson",
"RColorBrewer",
"ggpubr",
"optparse",
"qqman",
"glue",
"fs",
"R.utils",
"devtools",
"janitor",
"tidylog",
"RNOmni",
"viridis",
"CMplot",
"patchwork",
"duckdb",
"duckdbfs"
)

bioc_packages <- c(
"rtracklayer",
"OmicCircos"
)

install_cran <- function(pack) {
    if (!require(pack, character.only = TRUE, quietly = TRUE)) {
        install.packages(pack, repos = "http://cran.us.r-project.org")
        if (!require(pack, character.only = TRUE, quietly = TRUE)) {
            stop(paste("Failed to install CRAN package:", pack))
        }
    }
}

install_bioc <- function(pack) {
    if (!require(pack, character.only = TRUE, quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", repos = "http://cran.us.r-project.org")
        }
        BiocManager::install(pack, ask = FALSE)
        if (!require(pack, character.only = TRUE, quietly = TRUE)) {
            stop(paste("Failed to install Bioconductor package:", pack))
        }
    }
}

for (pack in cran_packages) {
    install_cran(pack)
}

for (pack in bioc_packages) {
    install_bioc(pack)
}
