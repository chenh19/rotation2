#!/bin/bash

echo -e "install.packages('devtools') \ndevtools::install_github('katsevich-lab/sceptre')" > ~/sceptre.R
sudo Rscript ~/sceptre.R
rm ~/sceptre.R
