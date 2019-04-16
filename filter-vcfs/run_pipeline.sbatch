#!/bin/bash
#SBATCH --job-name=filter-vcfs
#SBATCH --time=2-0:0:0
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --partition=control
#################

nextflow run main.nf -resume -c nextflow.config -with-timeline -with-trace -with-dag dag.png
