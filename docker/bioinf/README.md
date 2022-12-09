# BIOINFORMATICS Docker image

A Swiss knife for "bioinformatics" purposes, Ubuntu 22.04 base containing:

- python 3.10+
- python packages from [requirements.txt](requirements.txt)
- R 4 built with openblas, shared library
- R packages from [install_packages.R](install_packages.R)
- gawk
- datamash
- vim
- emacs
- jq
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
docker build -t eu.gcr.io/finngen-refinery-dev/bioinformatics:VERSION -f Dockerfile .
```

## Current image in GCR

The current image is `eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8`

Program versions and libraries in the current image (printed in the end of image build):

```
Linux bc503bdbaa1e 5.15.0-1022-gcp #29-Ubuntu SMP Mon Oct 24 12:50:24 UTC 2022 x86_64 x86_64 x86_64 GNU/Linux
Python 3.10.6
R version 4.2.2 (2022-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.1 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.2.2
PLINK v1.90b6.26 64-bit (2 Apr 2022)
PLINK v2.00a3.7LM AVX2 Intel (24 Oct 2022)

Welcome to qctool
(version: 2.2.0, revision: unknown)

(C) 2009-2020 University of Oxford

Usage: qctool <options>

Welcome to bgenix
(version: 1.1.7, revision )

(C) 2009-2017 University of Oxford

Usage: bgenix <options>

Welcome to cat-bgen

(C) 2009-2017 University of Oxford

Usage: cat-bgen <options>
bcftools 1.16
Using htslib 1.16
Copyright (C) 2022 Genome Research Ltd.
License Expat: The MIT/Expat license
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.
tabix (htslib) 1.16
Copyright (C) 2022 Genome Research Ltd.
bgzip (htslib) 1.16
Copyright (C) 2022 Genome Research Ltd.
GNU Awk 5.1.0, API: 3.0 (GNU MPFR 4.1.0, GNU MP 6.2.1)
Copyright (C) 1989, 1991-2020 Free Software Foundation.

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
datamash (GNU datamash) 1.7
Copyright (C) 2020 Assaf Gordon
License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>.
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.

Written by Assaf Gordon.
GNU Emacs 27.1
Copyright (C) 2020 Free Software Foundation, Inc.
GNU Emacs comes with ABSOLUTELY NO WARRANTY.
You may redistribute copies of GNU Emacs
under the terms of the GNU General Public License.
For more information about these matters, see the file named COPYING.
VIM - Vi IMproved 8.2 (2019 Dec 12, compiled Sep 13 2022 09:35:02)
Included patches: 1-3995, 4563, 4646, 4774, 4895, 4899, 4901, 4919
Modified by team+vim@tracker.debian.org
Compiled by team+vim@tracker.debian.org
Huge version without GUI.  Features included (+) or not (-):
+acl               +file_in_path      +mouse_urxvt       -tag_any_white
+arabic            +find_in_path      +mouse_xterm       -tcl
+autocmd           +float             +multi_byte        +termguicolors
+autochdir         +folding           +multi_lang        +terminal
-autoservername    -footer            -mzscheme          +terminfo
-balloon_eval      +fork()            +netbeans_intg     +termresponse
+balloon_eval_term +gettext           +num64             +textobjects
-browse            -hangul_input      +packages          +textprop
++builtin_terms    +iconv             +path_extra        +timers
+byte_offset       +insert_expand     -perl              +title
+channel           +ipv6              +persistent_undo   -toolbar
+cindent           +job               +popupwin          +user_commands
-clientserver      +jumplist          +postscript        +vartabs
-clipboard         +keymap            +printer           +vertsplit
+cmdline_compl     +lambda            +profile           +vim9script
+cmdline_hist      +langmap           -python            +viminfo
+cmdline_info      +libcall           +python3           +virtualedit
+comments          +linebreak         +quickfix          +visual
+conceal           +lispindent        +reltime           +visualextra
+cryptv            +listcmds          +rightleft         +vreplace
+cscope            +localmap          -ruby              +wildignore
+cursorbind        -lua               +scrollbind        +wildmenu
+cursorshape       +menu              +signs             +windows
+dialog_con        +mksession         +smartindent       +writebackup
+diff              +modify_fname      +sodium            -X11
+digraphs          +mouse             -sound             -xfontset
-dnd               -mouseshape        +spell             -xim
-ebcdic            +mouse_dec         +startuptime       -xpm
+emacs_tags        +mouse_gpm         +statusline        -xsmp
+eval              -mouse_jsbterm     -sun_workshop      -xterm_clipboard
+ex_extra          +mouse_netterm     +syntax            -xterm_save
+extra_search      +mouse_sgr         +tag_binary        
-farsi             -mouse_sysmouse    -tag_old_static    
   system vimrc file: "$VIM/vimrc"
     user vimrc file: "$HOME/.vimrc"
 2nd user vimrc file: "~/.vim/vimrc"
      user exrc file: "$HOME/.exrc"
       defaults file: "$VIMRUNTIME/defaults.vim"
  fall-back for $VIM: "/usr/share/vim"
Compilation: gcc -c -I. -Iproto -DHAVE_CONFIG_H -Wdate-time -g -O2 -ffile-prefix-map=/build/vim-NA7QBf/vim-8.2.3995=. -flto=auto -ffat-lto-objects -flto=auto -ffat-lto-objects -fstack-protector-strong -Wformat -Werror=format-security -D_REENTRANT -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=1 
Linking: gcc -Wl,-Bsymbolic-functions -flto=auto -ffat-lto-objects -flto=auto -Wl,-z,relro -Wl,-z,now -Wl,--as-needed -o vim -lm -ltinfo -lselinux -lsodium -lacl -lattr -lgpm -L/usr/lib/python3.10/config-3.10-x86_64-linux-gnu -lpython3.10 -lcrypt -ldl -lm -lm 
jq-1.6
Google Cloud SDK 411.0.0
bq 2.0.83
bundled-python3-unix 3.9.12
core 2022.12.05
gcloud-crc32c 1.0.0
gsutil 5.17
contourpy==1.0.6
cycler==0.11.0
dbus-python==1.2.18
fonttools==4.38.0
joblib==1.2.0
kiwisolver==1.4.4
matplotlib==3.6.2
numpy==1.23.5
packaging==22.0
pandas==1.5.2
Pillow==9.3.0
pip==22.0.2
PyGObject==3.42.1
pyparsing==3.0.9
python-dateutil==2.8.2
pytz==2022.6
scikit-learn==1.1.3
scipy==1.9.3
setuptools==59.6.0
six==1.16.0
threadpoolctl==3.1.0
wheel==0.37.1
                    Package    Version
abind                 abind      1.4-5
askpass             askpass        1.1
assertthat       assertthat      0.2.1
backports         backports      1.4.1
base64enc         base64enc      0.1-3
bit                     bit      4.0.5
bit64                 bit64      4.0.5
blob                   blob      1.2.3
brew                   brew      1.0-8
brio                   brio      1.1.3
broom                 broom      1.0.1
bslib                 bslib      0.4.1
cachem               cachem      1.0.6
calibrate         calibrate      1.7.7
callr                 callr      3.7.3
car                     car      3.1-1
carData             carData      3.0-5
cellranger       cellranger      1.1.0
cli                     cli      3.4.1
clipr                 clipr      0.8.0
clisymbols       clisymbols      1.2.0
CMplot               CMplot      4.2.0
colorspace       colorspace      2.0-3
commonmark       commonmark      1.8.1
corrplot           corrplot       0.92
cowplot             cowplot      1.1.1
cpp11                 cpp11      0.4.3
crayon               crayon      1.5.2
credentials     credentials      1.3.2
curl                   curl      4.3.3
data.table       data.table     1.14.6
DBI                     DBI      1.1.3
dbplyr               dbplyr      2.2.1
desc                   desc      1.4.2
diffobj             diffobj      0.3.5
digest               digest     0.6.30
downlit             downlit      0.4.2
dplyr                 dplyr     1.0.10
dtplyr               dtplyr      1.2.2
ellipsis           ellipsis      0.3.2
evaluate           evaluate       0.18
fansi                 fansi      1.0.3
farver               farver      2.1.1
fastmap             fastmap      1.1.0
fontawesome     fontawesome      0.4.0
forcats             forcats      0.5.2
fs                       fs      1.5.2
gargle               gargle      1.2.1
generics           generics      0.1.3
gert                   gert      1.9.2
getopt               getopt     1.20.3
ggplot2             ggplot2      3.4.0
ggpubr               ggpubr      0.5.0
ggrepel             ggrepel      0.9.2
ggsci                 ggsci        2.9
ggsignif           ggsignif      0.6.4
gh                       gh      1.3.1
gitcreds           gitcreds      0.1.2
glue                   glue      1.6.2
googledrive     googledrive      2.0.0
googlesheets4 googlesheets4      1.0.1
gridExtra         gridExtra        2.3
gtable               gtable      0.3.1
haven                 haven      2.5.1
highr                 highr        0.9
hms                     hms      1.1.2
htmltools         htmltools      0.5.4
htmlwidgets     htmlwidgets      1.5.4
httpuv               httpuv      1.6.6
httr                   httr      1.4.4
ids                     ids      1.0.1
ini                     ini      0.3.1
isoband             isoband      0.2.6
janitor             janitor      2.1.0
jquerylib         jquerylib      0.1.4
jsonlite           jsonlite      1.8.4
knitr                 knitr       1.41
labeling           labeling      0.4.2
later                 later      1.3.0
lifecycle         lifecycle      1.0.3
lme4                   lme4     1.1-31
lubridate         lubridate      1.9.0
magrittr           magrittr      2.0.3
MatrixModels   MatrixModels      0.5-1
memoise             memoise      2.0.1
mime                   mime       0.12
miniUI               miniUI    0.1.1.1
minqa                 minqa      1.2.5
modelr               modelr     0.1.10
munsell             munsell      0.5.0
nloptr               nloptr      2.0.3
numDeriv           numDeriv 2016.8-1.1
openssl             openssl      2.0.5
optparse           optparse      1.7.3
patchwork         patchwork      1.1.2
pbkrtest           pbkrtest      0.5.1
pillar               pillar      1.8.1
pkgbuild           pkgbuild      1.4.0
pkgconfig         pkgconfig      2.0.3
pkgload             pkgload      1.3.2
plyr                   plyr      1.8.8
polynom             polynom      1.4-1
praise               praise      1.0.0
prettyunits     prettyunits      1.1.1
processx           processx      3.8.0
profvis             profvis      0.3.7
progress           progress      1.2.2
promises           promises    1.2.0.1
ps                       ps      1.7.2
purrr                 purrr      0.3.5
qqman                 qqman      0.1.8
quantreg           quantreg       5.94
R.methodsS3     R.methodsS3      1.8.2
R.oo                   R.oo     1.25.0
R.utils             R.utils     2.12.2
R6                       R6      2.5.1
rappdirs           rappdirs      0.3.3
rcmdcheck         rcmdcheck      1.4.0
RColorBrewer   RColorBrewer      1.1-3
Rcpp                   Rcpp      1.0.9
RcppArmadillo RcppArmadillo 0.11.4.2.1
RcppEigen         RcppEigen  0.3.3.9.3
readr                 readr      2.1.3
readxl               readxl      1.4.1
rematch             rematch      1.0.1
rematch2           rematch2      2.1.2
remotes             remotes      2.4.2
reprex               reprex      2.0.2
rlang                 rlang      1.0.6
rmarkdown         rmarkdown       2.18
RNOmni               RNOmni      1.0.1
roxygen2           roxygen2      7.2.2
rprojroot         rprojroot      2.0.3
rstatix             rstatix      0.7.1
rstudioapi       rstudioapi       0.14
rversions         rversions      2.1.2
rvest                 rvest      1.0.3
sass                   sass      0.4.4
scales               scales      1.2.1
selectr             selectr      0.4-2
sessioninfo     sessioninfo      1.2.2
shiny                 shiny      1.7.3
snakecase         snakecase     0.11.0
sourcetools     sourcetools      0.1.7
SparseM             SparseM       1.81
stringi             stringi      1.7.8
stringr             stringr      1.5.0
sys                     sys      3.4.1
systemfonts     systemfonts      1.0.4
testthat           testthat      3.1.5
tibble               tibble      3.1.8
tidylog             tidylog      1.0.2
tidyr                 tidyr      1.2.1
tidyselect       tidyselect      1.2.0
tidyverse         tidyverse      1.3.2
timechange       timechange      0.1.1
tinytex             tinytex       0.42
tzdb                   tzdb      0.3.0
urlchecker       urlchecker      1.0.1
usethis             usethis      2.1.6
utf8                   utf8      1.2.2
uuid                   uuid      1.1-0
vctrs                 vctrs      0.5.1
viridis             viridis      0.6.2
viridisLite     viridisLite      0.4.1
vroom                 vroom      1.6.0
waldo                 waldo      0.4.0
whisker             whisker      0.4.1
withr                 withr      2.5.0
xfun                   xfun       0.35
xml2                   xml2      1.3.3
xopen                 xopen      1.0.0
xtable               xtable      1.8-4
yaml                   yaml      2.3.6
zip                     zip      2.2.2
```
