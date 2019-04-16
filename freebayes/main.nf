#!/usr/bin/env nextflow
// Create bams channel from samples file

sample_list = file(params.samples_file)
samples = Channel
  .from(sample_list.readLines())
  .map { sampleID ->
    bamfile = file("${params.bams_dir}/${sampleID}.md.bam")
    baifile = file("${params.bams_dir}/${sampleID}.md.bai")
    [ sampleID, bamfile, baifile ]
}

samples = samples.view()
samples.into{ bam_files; bai_files }
bam_files = bam_files.map{sampleID, bam_file, bai_file -> bam_file}
bai_files = bai_files.map{sampleID, bam_file, bai_file -> bai_file}

all_bam_files = bam_files.toList()
all_bai_files = bai_files.toList()

regionTasks = Channel
  .from(file(params.regionTasksFreebayes).readLines())
  .map {line ->
    list       = line.split(",")
    task       = list[0]
    seq_cumsum = list[1]
    flag       = list[2]
    [ task, flag ]
}

process Freebayes1 {
  tag "${params.pipeline_name}-${regionTask}"

  cpus { 1 }
  memory { task.attempt == 1 ? 12.GB: task.attempt == 2 ? 24.GB: 64.GB }
  time { 2.d }
  errorStrategy { 'retry' }
  maxRetries 2
  maxErrors '-1'

  input:
  file(bams) from all_bam_files.first()
  file(bais) from all_bai_files.first()
  set regionTask, regions from regionTasks

  output:
  set regionTask, file("region_${regionTask}.vcf.bgz"), file("*.tbi") into region_VCFs

  script:
  input_bams = bams.collect{"-b $it"}.join(' ')

  """
  set -e -o pipefail
  mkdir -p temp
  freebayes ${params.freebayes_options} \
    -f $params.genome \
    --gvcf \
    $input_bams \
    $regions | vcf-sort --temporary-directory temp | bgzip -@ ${task.cpus} > region_${regionTask}.vcf.bgz
    tabix -p vcf region_${regionTask}.vcf.bgz
  """
}

region_VCFs.into { region_VCFs; region_VCFs2 }
region_vcf_files = region_VCFs.map { id, file, fileindex -> file }
region_vcf_indexes = region_VCFs2.map { id, file, fileindex -> fileindex }
all_VCFs = region_vcf_files.toList()
all_indexes = region_vcf_indexes.toList()

process ConcatenateVCFs {
  tag "${params.pipeline_name}"
  publishDir "outputs"

  cpus {  task.attempt == 1 ? 8: 16  }
  memory { task.attempt == 1 ? 96.GB: 192.GB }
  time { 2.d }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 2
  maxErrors '-1'

  input:
  file(vcfs) from all_VCFs
  file(vcfindexes) from all_indexes

  output:
  set file("*.vcf.gz" ), file("*.tbi") into finalVCF

  script:
  input_vcfs = vcfs.collect{"$it"}.join(' ')

  """
  set -e -o pipefail
  mkdir -p temp
  ls *.vcf.bgz | awk -F '.' '{print \$1}'  | xargs -t -I '{}' cp {}.vcf.bgz {}.vcf.gz
  ls *.vcf.bgz.tbi | awk -F '.' '{print \$1}'  | xargs -t -I '{}' cp {}.vcf.bgz.tbi {}.vcf.gz.tbi
  vcf-concat region*.vcf.gz | vcf-sort --temporary-directory temp | bgzip -@ ${task.cpus} > ${params.output_prefix}.unfiltered.snps.indels.vcf.gz
  tabix -p vcf ${params.output_prefix}.unfiltered.snps.indels.vcf.gz
  """
}

workflow.onComplete {
  println "Pipeline completed at: $workflow.complete"
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
