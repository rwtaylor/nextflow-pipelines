#!/usr/bin/env nextflow

/* Multiplex pcr calling workflow
  Started September 2017

  @Authors
  Ryan Taylor <ryan@ryantaylor.net>

*/

// Add fastq files to the fastqFiles channel
fastqFiles = Channel.fromFilePairs("${params.fastq_files_path}")
//fastqFiles = Channel.fromPath("${params.fastq_files_path}")
fastqFiles = fastqFiles.view()

// Extract FileIDs from file names. This uses a regex defined in the nextflow.config file
fastqFiles = fastqFiles.map{ sample, reads ->
  fileID = sample =~ params.regex_fastq_ID
  [fileID[0][1], reads[0], reads[1]]
}

// Update FileID with SampleID if lookup file was given in config
if(params.sample_name_file != null){
  // Get sample info data
  sample_file = file("${params.sample_name_file}")
  sample_info = Channel.from(sample_file.readLines())
  sample_info = sample_info.map { line ->
    lsplit = line.split('\t')
    [lsplit[0], lsplit[1]]
  }

  // Update FileIDs with SampleIDs
  fastqFiles = fastqFiles.phase(sample_info).map {fq, si ->
    [si[1], fq[1], fq[2]]
  }
}


fastqFiles = fastqFiles.view()

////////////////////////////////////////////////////////////////////////////////////////////////////
/* Trimming */

// Trim fastq files if the config parameter 'trim_fastq_files' is true
if (params.trim_fastq_files) {

  process Trim_galore {
    publishDir path:"${params.publish_directory}/trimmed-fastqs", mode: "copy", overwrite: true

    tag "${params.output_prefix}-${sampleID}"

    cpus 1
    memory { task.cpus * 4.GB }
    
    input:
    set sampleID, file(r1), file(r2) from fastqFiles

    output:
    set sampleID, file("*_val_1.fq.gz"), file("*_val_2.fq.gz") into trimmedFastqs

    """
    /usr/local/bin/trim_galore --paired --length 5 --gzip ${r1} ${r2}
    """
  }

} else {

    trimmedFastqs = fastqFiles

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/* Mapping */

if (params.bam_files_path == null) {
  process Mapping_bwa {
    publishDir path:"${params.publish_directory}/bams", mode: "copy", overwrite: true
    tag "${params.output_prefix}-${sampleID}"

    cpus 2
    memory { task.cpus * 4.GB }

    input:
    set sampleID, file(fq1), file(fq2) from trimmedFastqs

    output:
    file("*.bam") into bwaMappedBams
    file("*.bam.bai") into bamIndexes

    script:
    readGroupString="\"@RG\\tID:${sampleID}\\tSM:${sampleID}\\tLB:${sampleID}\\tPL:illumina\""

    """
    set -eo pipefail
    /usr/local/bin/bwa mem -M -R ${readGroupString} -B 3 -t ${task.cpus} ${params.reference} ${fq1} ${fq2} | \
    /usr/local/bin/samtools view -hu - | /usr/local/bin/samtools sort --threads ${task.cpus} -O bam - > ${sampleID}.bam
    /usr/local/bin/samtools index ${sampleID}.bam
    """
  }
} else {
    bwaMappedBams = Channel.fromPath("${params.bam_files_path}/*.bam")
    bamIndexes = Channel.fromPath("${params.bam_files_path}/*.bam.bai")
}




// Pileup and calling

targetTasks = Channel
  .from(file(params.snp_locations_tab).readLines())
  .map { line ->
    list     = line.split("\t")
    name     = list[0]
    name     = name.replaceFirst(/:/, "_")
    scaffold = list[1]
    position = list[2]
    region = scaffold + ":" + position
    [ name, region ]
  }

process Pileup_call_target {
  container = '/zstor/containers/singularity/biobase.img'
  publishDir path:"${params.publish_directory}/vcfs", mode: "copy", overwrite: true
  tag "${params.output_prefix}-${target_name}"

  cpus 1
  memory { 8.GB }
  time { 6.h }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set target_name, target_region from targetTasks
  file(bams) from bwaMappedBams.toList()
  file(bais) from bamIndexes.toList()

  output:
  file("${params.output_prefix}.${target_name}.target.vcf.gz") into target_vcfs
  file("${params.output_prefix}.${target_name}.target.vcf.gz.tbi") into target_vcf_indexes

  """
  set -e -o pipefail
  mkdir -p temp
  /usr/local/bin/bcftools mpileup -r ${target_region} -a INFO/AD,FORMAT/AD,FORMAT/DP -Ou --max-depth 100000 -f ${params.reference} ${bams} |\
   bcftools call -Ou -m | bcftools sort --temp-dir temp -Oz -o ${params.output_prefix}.${target_name}.target.vcf.gz
   tabix -p vcf ${params.output_prefix}.${target_name}.target.vcf.gz
  """
}

target_vcfs = target_vcfs.view()

process ConcatenateTargetVCFs {
  publishDir path:"${params.publish_directory}", mode: "copy", overwrite: true

  cpus {  task.attempt == 1 ? 8: 16  }
  memory { task.attempt == 1 ? 96.GB: 192.GB }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 2
  maxErrors '-1'

  input:
  file(vcfs) from target_vcfs.toList()
  file(vcfindexes) from target_vcf_indexes.toList()

  output:
  set file("${params.output_prefix}.target.vcf.gz" ), file("${params.output_prefix}.target.vcf.gz.tbi") into target_vcf

  script:
  input_vcfs = vcfs.collect{"$it"}.join(' ')

  """
  set -e -o pipefail
  mkdir -p temp
  ls *.target.vcf.gz | awk -F '.' '{print \$1"."\$2}'  | xargs -t -I '{}' cp {}.target.vcf.gz {}.cptarget.vcf.gz
  ls *.target.vcf.gz.tbi | awk -F '.' '{print \$1"."\$2}'  | xargs -t -I '{}' cp {}.target.vcf.gz.tbi {}.cptarget.vcf.gz.tbi
  vcf-concat *.cptarget.vcf.gz | vcf-sort --temporary-directory temp | bgzip -@ ${task.cpus} > ${params.output_prefix}.target.vcf.gz
  tabix -p vcf ${params.output_prefix}.target.vcf.gz
  """
}

vcfs = target_vcf

process ConvertToTSV {
  publishDir path:"${params.publish_directory}", mode: "copy", overwrite: true
  container = '/zstor/containers/singularity/biobase.img'
  publishDir "${params.publish_directory}/tsvs", mode: 'copy'

  cpus 1

  input:
  set file(vcf), file(index) from vcfs

  output:
  file("*.tsv") into tsv_files

  script:
  output_prefix = vcf.baseName - ~/\.vcf*/

  """
  bcftools query -H -f '%CHROM\t%POS\t%INDEL\t%QUAL\t%REF\t%ALT{0}\t%DP[\t%PL:%GT:%AC]\n' -o ${output_prefix}.tsv $vcf
  """
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Pipeline duration: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
