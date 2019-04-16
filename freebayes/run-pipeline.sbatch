#!/bin/bash
#SBATCH --job-name=felCat9-180709-s66
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=7-0:0:0
#SBATCH --partition=dpetrov
#################

module load devel java

nextflow run main.nf -resume -c nextflow.config -with-timeline -with-trace -with-dag dag.png
