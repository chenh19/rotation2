#!/bin/bash

# install anaconda
# Side note: conda install kallisto (not necessary for kb)
[ ! -d ./shscript/ ] && mkdir ./shscript/
wget -O ./shscript/Anaconda-latest-Linux-x86_64.sh https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
bash ./shscript/Anaconda-latest-Linux-x86_64.sh && sleep 3
conda update anaconda -y && conda update --all -y
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install kallisto -y
conda config --set auto_activate_base false # disable auto activate base in terminal
#conda activate # activate base when needed
#rm -rf ~/anaconda3/ # uninstall anaconda
rm -rf ./shscript/
