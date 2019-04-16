# nextflow pipelines

Some handy-to-me bioinformatic pipelines. These use the amazing [nextflow DSL](https://nextflow.io), and mostly rely on [singularity](https://www.sylabs.io/docs/) containers for easy deployment locally or on a cluster. I add pipelines that I may reuse here. These pipelines are ment to be starting points for new analyses.

* **alignment** This pipeline uses [LAST](http://last.cbrc.jp) to align reference genomes.
* **bam-stats** A pipeline to get basic stats from bamfiles. This was developed for a multiplex pcr workflow.
* **fasta-qc** Very basic fastqc pipeline
* **filter-vcfs** Filtering VCFs from whole genome sequencing.
* **freebayes** A pipeline that calls SNPs from mapped and de-duplicated bams with freebayes. Included is a python script to generate regions from the reference genome for multiprocessing across many nodes.
* **mapping** A pipeline that maps fastq files to a reference using bwa mem. Then duplicate reads in the bams are marked with picardtools markduplicates.
* **mpcrseq** A pipeline for procesing multiplex-pcr sequence data.
* **trim-fastq** Split fastq files for parallel mapping, then trim fastq with trim_galore
