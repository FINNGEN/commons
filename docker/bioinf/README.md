# BIOINFORMATICS Docker image

A Swiss knife for "bioinformatics" purposes, Ubuntu 18.04 base containing:

- python 3.6+
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

The current image is `eu.gcr.io/finngen-refinery-dev/bioinformatics:0.7`

Program versions and libraries in the current image (printed in the end of image build):

```
Linux 21930c8a723a 4.15.0-1032-gcp #34-Ubuntu SMP Wed May 8 13:02:46 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux

Python 3.6.9                                                                                          
R version 4.0.5 (2021-03-31)                                           
Platform: x86_64-pc-linux-gnu (64-bit)                             
Running under: Ubuntu 18.04.5 LTS                                    
                                                                
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
[1] compiler_4.0.5
PLINK v1.90b6.21 64-bit (19 Oct 2020)
PLINK v2.00a3LM AVX2 Intel (13 Apr 2021)

Welcome to qctool
(version: 2.0.8, revision )

(C) 2009-2017 University of Oxford

Usage: qctool <options>

Welcome to bgenix
(version: 1.1.7, revision )

(C) 2009-2017 University of Oxford

Usage: bgenix <options>

Welcome to cat-bgen

(C) 2009-2017 University of Oxford

Usage: cat-bgen <options>
bcftools 1.12
Using htslib 1.12
Copyright (C) 2021 Genome Research Ltd.
License Expat: The MIT/Expat license
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.
tabix (htslib) 1.12
Copyright (C) 2021 Genome Research Ltd.
bgzip (htslib) 1.12
Copyright (C) 2021 Genome Research Ltd.
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
GNU Emacs 25.2.2
Copyright (C) 2017 Free Software Foundation, Inc.
GNU Emacs comes with ABSOLUTELY NO WARRANTY.
You may redistribute copies of GNU Emacs
under the terms of the GNU General Public License.
For more information about these matters, see the file named COPYING.
VIM - Vi IMproved 8.0 (2016 Sep 12, compiled Oct 13 2020 15:49:09)
Included patches: 1-1453
Modified by pkg-vim-maintainers@lists.alioth.debian.org
Compiled by pkg-vim-maintainers@lists.alioth.debian.org
Huge version without GUI.  Features included (+) or not (-):
+acl               +farsi             +mouse_sgr         -tag_any_white
+arabic            +file_in_path      -mouse_sysmouse    -tcl
+autocmd           +find_in_path      +mouse_urxvt       +termguicolors
-autoservername    +float             +mouse_xterm       +terminal
-balloon_eval      +folding           +multi_byte        +terminfo
+balloon_eval_term -footer            +multi_lang        +termresponse
-browse            +fork()            -mzscheme          +textobjects
++builtin_terms    +gettext           +netbeans_intg     +timers
+byte_offset       -hangul_input      +num64             +title
+channel           +iconv             +packages          -toolbar
+cindent           +insert_expand     +path_extra        +user_commands
-clientserver      +job               -perl              +vertsplit
-clipboard         +jumplist          +persistent_undo   +virtualedit
+cmdline_compl     +keymap            +postscript        +visual
+cmdline_hist      +lambda            +printer           +visualextra
+cmdline_info      +langmap           +profile           +viminfo
+comments          +libcall           -python            +vreplace
+conceal           +linebreak         +python3           +wildignore
+cryptv            +lispindent        +quickfix          +wildmenu
+cscope            +listcmds          +reltime           +windows
+cursorbind        +localmap          +rightleft         +writebackup
+cursorshape       -lua               -ruby              -X11
+dialog_con        +menu              +scrollbind        -xfontset
+diff              +mksession         +signs             -xim
+digraphs          +modify_fname      +smartindent       -xpm
-dnd               +mouse             +startuptime       -xsmp
-ebcdic            -mouseshape        +statusline        -xterm_clipboard
+emacs_tags        +mouse_dec         -sun_workshop      -xterm_save
+eval              +mouse_gpm         +syntax
+ex_extra          -mouse_jsbterm     +tag_binary
+extra_search      +mouse_netterm     +tag_old_static
   system vimrc file: "$VIM/vimrc"
     user vimrc file: "$HOME/.vimrc"
 2nd user vimrc file: "~/.vim/vimrc"
      user exrc file: "$HOME/.exrc"
       defaults file: "$VIMRUNTIME/defaults.vim"
  fall-back for $VIM: "/usr/share/vim"
Compilation: gcc -c -I. -Iproto -DHAVE_CONFIG_H   -Wdate-time  -g -O2 -fdebug-prefix-map=/build/vim-EfP9JP/vim-8.0.1453=. -fstack-protector-strong -Wformat -Werror=format-security -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=1
Linking: gcc   -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -Wl,--as-needed -o vim        -lm -ltinfo -lnsl  -lselinux  -lacl -lattr -lgpm -ldl     -L/usr/lib/python3.6/config-3.6m-x86_64-linux-gnu -lpython3.6m -lpthread -ldl -lutil -lm
jq-1.5-1-a5b5cbe
asn1crypto==0.24.0
cryptography==2.1.4
cycler==0.10.0
idna==2.6
joblib==1.0.1
keyring==10.6.0
keyrings.alt==3.0
kiwisolver==1.3.1
matplotlib==3.3.4
numpy==1.19.5
pandas==1.1.5
Pillow==8.2.0
pip==9.0.1
pycrypto==2.6.1
pygobject==3.26.1
pyparsing==2.4.7
python-dateutil==2.8.1
pytz==2021.1
pyxdg==0.25
scikit-learn==0.24.1
scipy==1.5.4
SecretStorage==2.3.1
setuptools==39.0.1
six==1.11.0
threadpoolctl==2.1.0
wheel==0.30.0
                    Package    Version
abind                 abind      1.4-5
askpass             askpass        1.1
assertthat       assertthat      0.2.1
backports         backports      1.2.1
base64enc         base64enc      0.1-3
BH                       BH   1.75.0-0
blob                   blob      1.2.1
brew                   brew      1.0-6
brio                   brio      1.1.1
broom                 broom      0.7.6
cachem               cachem      1.0.4
calibrate         calibrate      1.7.7
callr                 callr      3.6.0
car                     car     3.0-10
carData             carData      3.0-4
cellranger       cellranger      1.1.0
cli                     cli      2.4.0
clipr                 clipr      0.7.1
clisymbols       clisymbols      1.2.0
CMplot               CMplot      3.6.2
colorspace       colorspace      2.0-0
commonmark       commonmark        1.7
conquer             conquer      1.0.2
corrplot           corrplot       0.84
cowplot             cowplot      1.1.1
cpp11                 cpp11      0.2.7
crayon               crayon      1.4.1
credentials     credentials      1.3.0
curl                   curl        4.3
data.table       data.table     1.14.0
DBI                     DBI      1.1.1
dbplyr               dbplyr      2.1.1
desc                   desc      1.3.0
devtools           devtools      2.4.0
diffobj             diffobj      0.3.4
digest               digest     0.6.27
dplyr                 dplyr      1.0.5
ellipsis           ellipsis      0.3.1
evaluate           evaluate       0.14
fansi                 fansi      0.4.2
farver               farver      2.1.0
fastmap             fastmap      1.1.0
forcats             forcats      0.5.1
fs                       fs      1.5.0
generics           generics      0.1.0
gert                   gert      1.3.0
getopt               getopt     1.20.3
ggplot2             ggplot2      3.3.3
ggpubr               ggpubr      0.4.0
ggrepel             ggrepel      0.9.1
ggsci                 ggsci        2.9
ggsignif           ggsignif      0.6.1
gh                       gh      1.2.1
gitcreds           gitcreds      0.1.1
glue                   glue      1.4.2
gridExtra         gridExtra        2.3
gtable               gtable      0.3.0
haven                 haven      2.3.1
highr                 highr        0.8
hms                     hms      1.0.0
htmltools         htmltools    0.5.1.1
httr                   httr      1.4.2
ini                     ini      0.3.1
isoband             isoband      0.2.4
janitor             janitor      2.1.0
jsonlite           jsonlite      1.7.2
knitr                 knitr       1.31
labeling           labeling      0.4.2
lifecycle         lifecycle      1.0.0
lme4                   lme4     1.1-26
lubridate         lubridate     1.7.10
magrittr           magrittr      2.0.1
maptools           maptools      1.1-1
markdown           markdown        1.1
MatrixModels   MatrixModels      0.5-0
matrixStats     matrixStats     0.58.0
memoise             memoise      2.0.0
mime                   mime       0.10
minqa                 minqa      1.2.4
modelr               modelr      0.1.8
munsell             munsell      0.5.0
nloptr               nloptr    1.2.2.2
numDeriv           numDeriv 2016.8-1.1
openssl             openssl      1.4.3
openxlsx           openxlsx      4.2.3
optparse           optparse      1.6.6
patchwork         patchwork      1.1.1
pbkrtest           pbkrtest      0.5.1
pillar               pillar      1.6.0
pkgbuild           pkgbuild      1.2.0
pkgconfig         pkgconfig      2.0.3
pkgload             pkgload      1.2.1
plyr                   plyr      1.8.6
polynom             polynom      1.4-0
praise               praise      1.0.0
prettyunits     prettyunits      1.1.1
processx           processx      3.5.1
progress           progress      1.2.2
ps                       ps      1.6.0
purrr                 purrr      0.3.4
qqman                 qqman      0.1.4
quantreg           quantreg       5.85
R.methodsS3     R.methodsS3      1.8.1
R.oo                   R.oo     1.24.0
R.utils             R.utils     2.10.1
R6                       R6      2.5.0
rappdirs           rappdirs      0.3.3
rcmdcheck         rcmdcheck      1.3.3
RColorBrewer   RColorBrewer      1.1-2
Rcpp                   Rcpp      1.0.6
RcppArmadillo RcppArmadillo 0.10.4.0.0
RcppEigen         RcppEigen  0.3.3.9.1
readr                 readr      1.4.0
readxl               readxl      1.3.1
rematch             rematch      1.0.1
rematch2           rematch2      2.1.2
remotes             remotes      2.3.0
reprex               reprex      2.0.0
rio                     rio     0.5.26
rlang                 rlang     0.4.10
rmarkdown         rmarkdown        2.7
RNOmni               RNOmni      1.0.0
roxygen2           roxygen2      7.1.1
rprojroot         rprojroot      2.0.2
rstatix             rstatix      0.7.0
rstudioapi       rstudioapi       0.13
rversions         rversions      2.0.2
rvest                 rvest      1.0.0
scales               scales      1.1.1
selectr             selectr      0.4-2
sessioninfo     sessioninfo      1.1.1
snakecase         snakecase     0.11.0
sp                       sp      1.4-5
SparseM             SparseM       1.81
statmod             statmod     1.4.35
stringi             stringi      1.5.3
stringr             stringr      1.4.0
sys                     sys        3.4
testthat           testthat      3.0.2
tibble               tibble      3.1.0
tidylog             tidylog      1.0.2
tidyr                 tidyr      1.1.3
tidyselect       tidyselect      1.1.0
tidyverse         tidyverse      1.3.0
tinytex             tinytex       0.31
usethis             usethis      2.0.1
utf8                   utf8      1.2.1
vctrs                 vctrs      0.3.7
viridis             viridis      0.5.1
viridisLite     viridisLite      0.4.0
waldo                 waldo      0.2.5
whisker             whisker        0.4
withr                 withr      2.4.1
xfun                   xfun       0.22
xml2                   xml2      1.3.2
xopen                 xopen      1.0.0
yaml                   yaml      2.2.1
zip                     zip      2.1.1
```
