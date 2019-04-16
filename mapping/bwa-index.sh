#!/bin/bash
#SBATCH --job-name=bwa-index
#SBATCH --time=2-0:0:0
#SBATCH --ntasks=1 --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=control,normal
#################

singularity exec /zstor/containers/singularity/mapping.img bwa index /zstor/pcg/2016-tiger-wgs/mapping/felCat9/GCF_000181335.3_Felis_catus_9.0_genomic.fna.gz

