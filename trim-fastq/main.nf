
samples = file(params.samples)
fastqFiles = Channel
  .from(samples.readLines())
  .map { line ->
    lsplit      = line.split()
    sampleID    = lsplit[0]
    libID       = lsplit[1]
    laneID      = lsplit[2]
    fastqFile1  = file(lsplit[3]+ "/" + lsplit[4])
    fastqFile2  = file(lsplit[3]+ "/" + lsplit[5])
    [ sampleID, libID, laneID, fastqFile1, fastqFile2 ]
}

process SplitFastq {
  tag "${task.attempt}.${params.pipeline_name}-${sampleID}-${libID}-${laneID}"

  cpus 4
  time { 1.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, file(fq1), file(fq2) from fastqFiles

  output:
  set sampleID, libID, laneID, file("*_R1.*.fq.gz"), file("*_R2.*.fq.gz") into splitFastqs

  script:
  nsplit = 4*params.fastq_chunksize

  """
  zcat ${fq1} | seqtk seq -l0 - | split -d --additional-suffix=.fq -l $nsplit --filter='pigz -p${task.cpus} > \$FILE.gz' - ${sampleID}_${libID}_${laneID}_R1.
  zcat ${fq2} | seqtk seq -l0 - | split -d --additional-suffix=.fq -l $nsplit --filter='pigz -p${task.cpus} > \$FILE.gz' - ${sampleID}_${libID}_${laneID}_R2.
  """
}

splitFastqs = splitFastqs.view()

//[BEN_CI10, mg3, L001, /scratch/PI/dpetrov/2016-tiger-wgs/fastq-data/pipelines/tfq-amu/work/ad/1f09e319b30206700fee36ad9fe9c7/BEN_CI10_mg3_L001_R1.00.fq.gz, /scratch/PI/dpetrov/2016-tiger-wgs/fastq-data/pipelines/tfq-amu/work/ad/1f09e319b30206700fee36ad9fe9c7/BEN_CI10_mg3_L001_R2.00.fq.gz]

splitFastqs = splitFastqs.flatMap{ sampleID, libID, laneID, fq1, fq2 ->
  nfiles = fq1.size() - 1
  out = []
  for (splitID in 0..nfiles ) {
    splitID_fq1 = fq1[splitID].baseName.split(/\./)[1].toInteger()
    splitID_fq2 = fq2[splitID].baseName.split(/\./)[1].toInteger()
    if (splitID_fq1 != splitID_fq2) exit 1, "file splitIDs do not match ${splitID} ${splitID_fq1} ${splitID_fq2}"
    if (splitID_fq1 != splitID) exit 1, "splitID does not match file ID ${splitID} ${splitID_fq1} ${splitID_fq2}"

    out.add([sampleID, libID, laneID, splitID, fq1[splitID], fq2[splitID]])
  }
  out
}

//splitFastqs = splitFastqs.view()

process Trim_galore {
  publishDir path:"${params.publish_directory}/${sampleID}", mode: 'copy', overwrite: true

  tag "${task.attempt}.${params.pipeline_name}-${sampleID}-${libID}-${laneID}-s.${splitID}"

  cpus { task.attempt == 1 ? 2: 4 }
  memory { task.attempt == 1 ? 8.GB: 16.GB }
  time { 1.h * task.attempt }
  errorStrategy { 'retry' }
  maxRetries 5
  maxErrors '-1'

  input:
  set sampleID, libID, laneID, splitID, file(fq1), file(fq2) from splitFastqs

  output:
  set sampleID, libID, laneID, splitID, file("*_val_1.fq.gz"), file("*_val_2.fq.gz") into trimmedFastqs

  """
  /usr/local/bin/trim_galore --paired --length 20 --gzip ${fq1} ${fq2}
  """
}

process FQ_validate {
  publishDir path:"${params.publish_directory}/${sampleID}", mode: 'copy', overwrite: true

  tag "${task.attempt}.${params.pipeline_name}-${sampleID}-${libID}-${laneID}-s.${splitID}"

  cpus 1
  time { 1.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'
  validExitStatus 0,1
  
  input:
  set sampleID, libID, laneID, splitID, file(fq1), file(fq2) from trimmedFastqs

  output:
  set sampleID, libID, laneID, splitID, file("fqvalidate.txt") into fastqvalidation

  """
  /usr/local/bin/fqtools validate ${fq1} ${fq2} &> fqvalidate.txt
  """
}

fqvalidation = fastqvalidation.map{sampleID, libID, laneID, splitID, fqvalfile ->
  fqval = fqvalfile.text.trim()
  [fqval, sampleID, libID, laneID, splitID]
}

fqvalidation.into{ fqvalidation; fqfail }

fqfail = fqfail.filter{ it[4] != "OK"}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    fqvalidation.view()
    fqfail.view()
}
