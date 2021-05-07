#!/usr/bin/env Rscript

req_packages <- c(
"data.table",
"tidyverse",
"RColorBrewer",
"ggplot2",
"ggpubr",
"optparse",
"qqman",
"purrr",
"glue",
"fs",
"R.utils",
"devtools",
"janitor",
"tidylog",
"RNOmni",
"rtracklayer",
"viridis",
"CMplot",
"patchwork",
"OmicCircos")

for (pack in req_packages) {
    if(!require(pack,character.only = TRUE)) {
        install.packages(pack, repos = "http://cran.us.r-project.org")
    }
}
