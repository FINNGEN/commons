# BIOINFORMATICS Docker image

A Swiss knife for "bioinformatics" purposes, Ubuntu 18 base containing:

- python 3.6+
- python packages from [requirements.txt](requirements.txt)
- R built with openblas, shared library
- R packages from [install_packages.R](install_packages.R)
- gawk
- datamash
- plink 1.9
- plink 2
- bcftools
- qctool
- bgenix
- tabix
- bgzip
- UCSC liftOver (build 37/38 chains: /liftover/)

## Building

```
git clone https://github.com/FINNGEN/commons
cd commons/docker/bioinf
docker build -t gcr.io/finngen-refinery-dev/bioinformatics:VERSION -f Dockerfile .
```

## Current image in GCR

The current image is `gcr.io/finngen-refinery-dev/bioinformatics:0.6`

Program versions and libraries in the current image (printed in the end of image build):

```
Linux 854bcdba8d0c 4.15.0-1032-gcp #34-Ubuntu SMP Wed May 8 13:02:46 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux
Python 3.6.9
R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.4 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.0.2
PLINK v1.90b6.18 64-bit (16 Jun 2020)
PLINK v2.00a3LM AVX2 Intel (1 Jul 2020)
[91m
Welcome to qctool[0m[91m
(version: 2.0.8, revision [0m[91m)
[0m[91m
(C) 2009-2017 University of Oxford

[0mUsage: qctool <options>
[91m
Welcome to bgenix
(version: 1.1.7, revision )

(C) 2009-2017 University of Oxford

[0mUsage: bgenix <options>
[91m
Welcome to cat-bgen

(C) 2009-2017 University of Oxford

[0mUsage: cat-bgen <options>
bcftools 1.10.2
Using htslib 1.10.2
Copyright (C) 2019 Genome Research Ltd.
License Expat: The MIT/Expat license
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.
tabix (htslib) 1.10.2
Copyright (C) 2019 Genome Research Ltd.
bgzip (htslib) 1.10.2
Copyright (C) 2019 Genome Research Ltd.
GNU Awk 4.1.4, API: 1.1 (GNU MPFR 4.0.1, GNU MP 6.1.2)
Copyright (C) 1989, 1991-2016 Free Software Foundation.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.
datamash (GNU datamash) 1.2
Copyright (C) 2017 Assaf Gordon
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.

Written by Assaf Gordon.
asn1crypto==0.24.0
cryptography==2.1.4
idna==2.6
joblib==0.16.0
keyring==10.6.0
keyrings.alt==3.0
numpy==1.19.0
pandas==1.0.5
pip==9.0.1
pycrypto==2.6.1
pygobject==3.26.1
python-dateutil==2.8.1
pytz==2020.1
pyxdg==0.25
scikit-learn==0.23.0
scipy==1.5.0
SecretStorage==2.3.1
setuptools==39.0.1
six==1.11.0
threadpoolctl==2.1.0
wheel==0.30.0
                  Package   Version
abind               abind     1.4-5
askpass           askpass       1.1
assertthat     assertthat     0.2.1
backports       backports     1.1.8
base64enc       base64enc     0.1-3
BH                     BH  1.72.0-3
blob                 blob     1.2.1
brew                 brew     1.0-6
broom               broom     0.5.6
calibrate       calibrate     1.7.7
callr               callr     3.4.3
car                   car     3.0-8
carData           carData     3.0-4
cellranger     cellranger     1.1.0
cli                   cli     2.0.2
clipr               clipr     0.7.0
colorspace     colorspace     1.4-1
commonmark     commonmark       1.7
corrplot         corrplot      0.84
covr                 covr     3.5.0
cowplot           cowplot     1.0.0
crayon             crayon     1.3.4
crosstalk       crosstalk   1.1.0.1
curl                 curl       4.3
data.table     data.table    1.12.8
DBI                   DBI     1.1.0
dbplyr             dbplyr     1.4.4
desc                 desc     1.2.0
devtools         devtools     2.3.0
digest             digest    0.6.25
dplyr               dplyr     1.0.0
DT                     DT      0.14
ellipsis         ellipsis     0.3.1
evaluate         evaluate      0.14
fansi               fansi     0.4.1
farver             farver     2.0.3
forcats           forcats     0.5.0
fs                     fs     1.4.2
generics         generics     0.0.2
getopt             getopt    1.20.3
ggplot2           ggplot2     3.3.2
ggpubr             ggpubr     0.4.0
ggrepel           ggrepel     0.8.2
ggsci               ggsci       2.9
ggsignif         ggsignif     0.6.0
gh                     gh     1.1.0
git2r               git2r    0.27.1
glue                 glue     1.4.1
gridExtra       gridExtra       2.3
gtable             gtable     0.3.0
haven               haven     2.3.1
highr               highr       0.8
hms                   hms     0.5.3
htmltools       htmltools     0.5.0
htmlwidgets   htmlwidgets     1.5.1
httr                 httr     1.4.1
ini                   ini     0.3.1
isoband           isoband     0.2.2
jsonlite         jsonlite     1.7.0
knitr               knitr      1.29
labeling         labeling       0.3
later               later   1.1.0.1
lazyeval         lazyeval     0.2.2
lifecycle       lifecycle     0.2.0
lme4                 lme4    1.1-23
lubridate       lubridate     1.7.9
magrittr         magrittr       1.5
maptools         maptools     1.0-1
markdown         markdown       1.1
MatrixModels MatrixModels     0.4-1
memoise           memoise     1.1.0
mime                 mime       0.9
minqa               minqa     1.2.4
modelr             modelr     0.1.8
munsell           munsell     0.5.0
nloptr             nloptr   1.2.2.2
openssl           openssl     1.4.2
openxlsx         openxlsx     4.1.5
optparse         optparse     1.6.6
pbkrtest         pbkrtest   0.4-8.6
pillar             pillar     1.4.4
pkgbuild         pkgbuild     1.0.8
pkgconfig       pkgconfig     2.0.3
pkgload           pkgload     1.1.0
plyr                 plyr     1.8.6
polynom           polynom     1.4-0
praise             praise     1.0.0
prettyunits   prettyunits     1.1.1
processx         processx     3.4.2
progress         progress     1.2.2
promises         promises     1.1.1
ps                     ps     1.3.3
purrr               purrr     0.3.4
qqman               qqman     0.1.4
quantreg         quantreg      5.55
R.methodsS3   R.methodsS3     1.8.0
R.oo                 R.oo    1.23.0
R.utils           R.utils     2.9.2
R6                     R6     2.4.1
rcmdcheck       rcmdcheck     1.3.3
RColorBrewer RColorBrewer     1.1-2
Rcpp                 Rcpp   1.0.4.6
RcppEigen       RcppEigen 0.3.3.7.0
readr               readr     1.3.1
readxl             readxl     1.3.1
rematch           rematch     1.0.1
rematch2         rematch2     2.1.2
remotes           remotes     2.1.1
reprex             reprex     0.3.0
reshape2         reshape2     1.4.4
rex                   rex     1.2.0
rio                   rio    0.5.16
rlang               rlang     0.4.6
rmarkdown       rmarkdown       2.3
roxygen2         roxygen2     7.1.1
rprojroot       rprojroot     1.3-2
rstatix           rstatix     0.6.0
rstudioapi     rstudioapi      0.11
rversions       rversions     2.0.2
rvest               rvest     0.3.5
scales             scales     1.1.1
selectr           selectr     0.4-2
sessioninfo   sessioninfo     1.1.1
sp                     sp     1.4-2
SparseM           SparseM      1.78
statmod           statmod    1.4.34
stringi           stringi     1.4.6
stringr           stringr     1.4.0
sys                   sys       3.3
testthat         testthat     2.3.2
tibble             tibble     3.0.1
tidyr               tidyr     1.1.0
tidyselect     tidyselect     1.1.0
tidyverse       tidyverse     1.3.0
tinytex           tinytex      0.24
usethis           usethis     1.6.1
utf8                 utf8     1.1.4
vctrs               vctrs     0.3.1
viridisLite   viridisLite     0.3.0
whisker           whisker       0.4
withr               withr     2.2.0
xfun                 xfun      0.15
xml2                 xml2     1.3.2
xopen               xopen     1.0.0
yaml                 yaml     2.2.1
zip                   zip     2.0.4
```