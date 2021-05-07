#!/usr/bin/env bash

set -e

uname -a
python3 --version
Rscript -e "sessionInfo()"
plink --version
plink2 --version
qctool -help | head -1
bgenix -help | head -1
cat-bgen -help | head -1
bcftools --version
tabix --version
bgzip --version
gawk --version
datamash --version
<<<<<<< HEAD
=======
emacs --version
vim --version
jq --version
>>>>>>> master
pip3 list --format=freeze
print_packages.R
