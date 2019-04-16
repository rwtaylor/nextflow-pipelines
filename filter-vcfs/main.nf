#!/usr/bin/env nextflow

vcfgz = ~/\.vcf.*/

input_vcfs = Channel.from(params.input_vcfs).map{it -> file(it) }.map{ file -> 
  prefix = params.output_prefix
  [prefix, file] }

//print("${params.subsample_ns[1][1]}")

process FilterVCF {
  publishDir "${params.publish_dir}/${prefix}", mode: 'copy'
  tag { prefix + "-snp-" + filterset.filter_name }
  cpus 2
  memory 8.GB
  time 6.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, file(vcf) from input_vcfs
  each filterset from params.filtersets
  
  output:
  set val("${prefix}-snp-${filterset.filter_name}"), file("*.vcf.gz"), file("*.tbi"), file("*.gbi"), prefix into filtered_vcfs

  script:
  removeindividuals = params.excludesamples.collect{"--remove-indv $it"}.join(" ")

  """
  set -e -o pipefail
  mkdir -p temp
  
  # Rename samples so that none have "_" because this borks PLINK
  /usr/local/bin/vcftools --gzvcf ${vcf} --recode --recode-INFO-all --remove-indels ${removeindividuals}\
  --minQ ${filterset.minqual} --minGQ ${filterset.mingq} --hwe ${filterset.hwe} \
  --maf ${filterset.maf} --mac ${filterset.mac} --max-alleles ${filterset.maxa} \
  --stdout | vcf-sort --temporary-directory temp | bgzip -@ ${task.cpus} > ${prefix}-snp-${filterset.filter_name}.vcf.gz
  tabix -p vcf ${prefix}-snp-${filterset.filter_name}.vcf.gz
  grabix index $vcf
  """
}

filtered_vcfs = filtered_vcfs.view()

filtered_vcfs.into { filtered_vcfs; filtered_vcfs_2}

process SampleVCF {
  publishDir "${params.publish_dir}/${output_folder}/subsampled", mode: 'copy'
  tag { prefix + "-ss" + subsamplerate }
  cpus 2
  memory 8.GB
  time 6.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, file(vcf), file(tbindex), file(gbindex), output_folder from filtered_vcfs
  each subn from params.subsample_ns

  output:
  set val("${prefix}-ss${subn[1]}"), file("*.vcf.gz"), file("*.tbi"), file("*.gbi"), output_folder into subsampled_vcfs

  """
  set -e -o pipefail
  mkdir -p temp
  grabix index $vcf
  grabix random $vcf ${subn[0]} | bcftools sort --temp-dir temp -O z -o ${prefix}-ss${subn[1]}.vcf.gz
  tabix -p vcf ${prefix}-ss${subn[1]}.vcf.gz
  grabix index ${prefix}-ss${subn[1]}.vcf.gz
  """
}

vcfs_to_rename = filtered_vcfs_2.mix(subsampled_vcfs)
vcfs_to_rename = vcfs_to_rename.view()
// Plink requires scaffolds to have a character prefix, so if they are integers then rename
// prepending "0" makes the scaff a non-integer (at least as far as plink is concerned), while
// also allowing LSAK (which requires integers) to interpret the scaff ID as an integer... tricky but convenient.

process RenameChromosomes {
  publishDir "${params.publish_dir}/${output_folder}/chrenamed", mode: 'copy'
  tag { prefix }
  cpus 2
  memory 8.GB
  time 6.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, file(vcf), file(tbindex), file(gbindex), output_folder from vcfs_to_rename
  
  output:
  set val("${prefix}-chrenamed"), file("*.vcf.gz"), file("*.tbi"), file("*.gbi"), output_folder into renamed_vcfs

  """
  set -e -o pipefail
  mkdir -p temp
  # 1) Get chromosome IDs, and output chromosome plus renamed chromosome to text file
  #    This strips all non-numeric characters from the chromosome ID. May cause issues if non-numeric characters
  #    are important for name uniquness...
  set -e
  zcat $vcf | grep -oP '^##contig=<ID=.*' | \
    sed -e 's/^##contig=<ID=\\(.*\\),.*/\\1/gm' | \
    awk '{ printf \$1 " "; gsub(/[A-Z_.]/,"", \$1); print 0\$1}' \
    > scaffs.txt
  # 2) Use BCFtools to rename chromosomes in VCF
  bcftools annotate --rename-chrs scaffs.txt $vcf | vcf-sort --chromosomal-order --temporary-directory temp | bgzip -@ ${task.cpus} > ${prefix}-chrenamed.vcf.gz
  tabix -p vcf ${prefix}-chrenamed.vcf.gz
  grabix index ${prefix}-chrenamed.vcf.gz
  """
}

workflow.onComplete {
  println "Pipeline completed at: $workflow.complete"
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
