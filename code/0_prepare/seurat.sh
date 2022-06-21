#!/bin/bash

sudo apt-get update && sudo apt-get install libgeos-dev -y
echo -e "install.packages(c('Seurat'))" > ~/seurat.R
sudo Rscript ~/seurat.R
rm ~/seurat.R
