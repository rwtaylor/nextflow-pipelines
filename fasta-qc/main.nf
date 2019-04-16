
fastqFiles = Channel
  .fromFilePairs("${params.fastq_directory}/*_{1,2}.{fq,fastq}.gz")
  .map { name, reads ->
    namesplit      = name.split("_")
    sampleID    = namesplit[0]
    libID       = namesplit[1]
    laneID      = namesplit[2] + "_" + namesplit[3]
    [ sampleID, libID, laneID, reads ]
}

fastqFiles = fastqFiles.view()

process FastQC {
  tag "${task.attempt}.${params.pipeline_name}-${sampleID}-${libID}-${laneID}"

  cpus 1
  memory 2.GB

  input:
  set sampleID, libID, laneID, file(reads) from fastqFiles

  output:
  file('*_fastqc.{zip,html}') into fastqc_results

  """
  /usr/local/bin/fastqc -q ${reads}
  """
}

fastqc_results = fastqc_results.collect().view()

process SampleMultiQC {

  publishDir "${params.publish_directory}", mode: 'copy', overwrite: true

  cpus 1
  memory 2.GB

  input:
  file ('fastqc/*') from fastqc_results

  output:
  file(qc_report) into multiqc_output

  """
  /usr/local/bin/multiqc -f -o qc-report fastqc/.
  """
}
