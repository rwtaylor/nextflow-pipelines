#!/bin/bash
#SBATCH --job-name=map-felCat9
#SBATCH --time=12-0:0:0
#SBATCH --partition=control
#################

nextflow run main.nf -resume -c nextflow.config --samples samples.tsv -with-timeline -with-trace -with-dag pipeline.png
