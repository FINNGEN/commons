#!/usr/bin/env Rscript

<<<<<<< HEAD
req_packages <- c("data.table", "tidyverse", "RColorBrewer", "ggplot2", "ggpubr", "optparse", "qqman", "purrr", "glue", "fs", "R.utils", "devtools")
=======
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

>>>>>>> master
for (pack in req_packages) {
    if(!require(pack,character.only = TRUE)) {
        install.packages(pack, repos = "http://cran.us.r-project.org")
    }
}
