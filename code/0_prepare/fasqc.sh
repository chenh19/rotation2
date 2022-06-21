#!/bin/bash

# install fastqc
sudo apt-get update && sudo apt-get install fastqc -y

# run fastqc in parallel
cd ~/rotation2/
[ ! -d ../Morris_2021/fastqc_results/ ] && mkdir ../Morris_2021/fastqc_results/
fastqc --threads 16 ../Morris_2021/STINGseq_Morris_2021_raw/*.fastq.gz --outdir=../Morris_2021/fastqc_results/
