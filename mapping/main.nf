#!/usr/bin/env nextflow

/* Mapping
  Started October 2016

  @Authors
  Ryan Taylor <rwtaylor@stanford.edu>

*/

genome_files = Channel.fromPath("${params.genome_fasta}*").toSortedList()
genome_files = genome_files.view()
genome_files.into{ genome_bwa; genome_markduplicates}

fastqFilesList = []
params.sampleIDs.each {
  sampleDir = file(params.sample_root_dir + "/" + it + "/")
  sampleFiles = sampleDir.eachFileMatch(~/.*\.fq\.gz/) { item ->
    if( item.isFile() ) {
        fastqFilesList.add(item)
    }
  }
}

fastqFiles = Channel.from(fastqFilesList)

fastqInfo = fastqFiles.map { item ->
  fsplit = item.baseName.split('_')
  if(fsplit.length == 6){
    sampleID = fsplit[0]
    libID = fsplit[1]
    laneID = fsplit[2]
    r = fsplit[3]
    rsplit = r.split(/\./)
    splitID = rsplit[1]
  } else if(fsplit.length == 7){
    sampleID = fsplit[0] + '_' + fsplit[1]
    libID = fsplit[2]
    laneID = fsplit[3]
    r = fsplit[4]
    rsplit = r.split(/\./)
    splitID = rsplit[1]
  }
  fq1 = file(params.sample_root_dir + "/" + sampleID + "/" + sampleID + "_" + libID + "_" + laneID + "_R1." + splitID + "_val_1.fq.gz")
  fq2 = file(params.sample_root_dir + "/" + sampleID + "/" + sampleID + "_" + libID + "_" + laneID + "_R2." + splitID + "_val_2.fq.gz")
  [sampleID, libID, laneID, splitID, fq1, fq2]
}

fastqInfo = fastqInfo.unique()
//fastqInfo = fastqInfo.view()

process Mapping_bwa {
  tag "${task.attempt}.${params.pipeline_name}-${sampleID}-${libID}-${laneID}-s.${splitID}"
  publishDir path:"${params.publish_directory}/${params.output_prefix}-bams"

  cpus { 4 }
  memory { task.cpus * 4.GB}
  time { 12.h }
  errorStrategy { 'finish' }

  input:
  file(genome) from genome_bwa
  set sampleID, libID, laneID, splitID, file(fq1), file(fq2) from fastqInfo

  output:
  set sampleID, libID, laneID, file("*.bam") into bwaMappedBams

  script:
  readGroupString="\"@RG\\tID:${sampleID}_${libID}_${laneID}\\tSM:${sampleID}\\tLB:${sampleID}_${libID}\\tPL:illumina\""

  """
  set -eo pipefail
  /usr/local/bin/bwa mem -M -R ${readGroupString} -t ${task.cpus} ${genome[0]} ${fq1} ${fq2} | \
  /usr/local/bin/samtools view -hu - | /usr/local/bin/samtools sort --threads $task.cpus -O bam - > ${sampleID}_${libID}_${laneID}_${splitID}.bam
  """
}

markduplicatesBams_samples = bwaMappedBams.groupTuple(by: 0).map {item -> [item[0], item[3]]}.view()

process MarkDuplicates {
  publishDir path:"${params.publish_directory}/${params.output_prefix}-markduplicates-bams", mode: 'copy', overwrite: true, saveAs: { it == "*.txt" ? "qc/sample-markduplicates-reports/$it" : "marked-duplicates-bams/$it" }
  tag "${task.attempt}.${params.pipeline_name}-${sampleID}"

  cpus { 8 }
  memory { task.cpus * 8.GB}
  time { 7.d }
  errorStrategy { 'finish' }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file(genome) from genome_markduplicates
  set sampleID, file(split_sample_bams) from markduplicatesBams_samples

  output:
  set sampleID, file("*.md.bam"), file("*.md.bai") into sampleBams
  set sampleID, file("*.markduplicates.samples.txt") into markduplicates_results

  script:
  split_sample_bams = split_sample_bams.collect{"I=$it"}.join(' ')

  """
  mkdir -p picard_tmp
  /usr/bin/java -jar /usr/local/opt/picard.jar MarkDuplicates \
    TMP_DIR=picard_tmp \
    ${split_sample_bams} \
    O=${sampleID}.md.bam \
    M=${sampleID}.markduplicates.samples.txt \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=${params.optical_duplicate_pixel_distance} \
    CREATE_INDEX=true
  """
}

sampleBams = sampleBams.view()

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
