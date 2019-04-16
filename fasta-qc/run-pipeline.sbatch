#!/bin/bash
#SBATCH --job-name=fastq-qc-pipeline
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --partition=control
#################

nextflow run main.nf -resume -c nextflow.config --samples samples.tsv -with-timeline -with-trace -with-dag pipeline.png
