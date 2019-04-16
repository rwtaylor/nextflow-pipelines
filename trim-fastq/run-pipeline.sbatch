#!/bin/bash
#SBATCH --job-name=p-tfq-ben2
#SBATCH --time=2-0:0:0
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --partition=dpetrov,hns
#################

nextflow run main.nf -resume -c nextflow.config --samples samples.tsv -with-timeline -with-trace -with-dag pipeline.png
