#!/usr/bin/env nextflow

/* Aligns genome 1 to genome 2 using LAST, calculates statistics. See config file for parameters.
  Started October 2016

  @Authors
  Ryan Taylor <ryan@ryantaylor.net>

  To run:
  $ nextflow run align_genomes.nf -c <nextflow.config>
*/


queries = file(params.query)
query_fasta = Channel
  .from(queries.readLines())
  .map {line ->
    lsplit      = line.split()
    queryID     = lsplit[0]
    fastaFile   = file(lsplit[1])
    [ queryID, fastaFile ]
}

query_fasta = query_fasta.view()

ref_genome = Channel.fromPath(params.ref_genome).map{ file -> [file.baseName, file] }

ref_genome = ref_genome.view()

if(params.minseqlen > 0) {

  process FilterRef {
    tag { prefix }
    publishDir "outputs", mode: 'copy'

    cpus 12
    memory { 48.GB }
    time { 6.h }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 5
    maxErrors '-1'

    input:
    set val(prefix), file(fasta) from ref_genome

    output:
    set prefix, file('*.fasta') into filtered_fasta

    """
    set -e
    /bin/cat $fasta | /usr/local/bin/seqkit seq -m $params.minseqlen > ${prefix}.filtered.fasta 
    """
  }

  filtered_fasta.into { ref_genome2; ref_genome3; ref_genome4 }
}
else {
  ref_genome.into { ref_genome2; ref_genome3; ref_genome4 }
}

process LastDB {
  tag { prefix }
  publishDir "outputs/lastdb"

  cpus 8
  memory { 32.GB }
  time { 6.h }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set val(prefix), file(fasta) from ref_genome2

  output:
  file('refdb*') into database_files

  """
  set -e
  /bin/cat ${fasta} | /usr/local/bin/lastdb -v -P${task.cpus} ${params.lastdb_options} refdb 
  """
}

/*process LastTrain {
  publishDir "outputs/stages/training"

  cpus 32
  memory { 4.GB * task.cpus }
  time { 2.d }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file fasta1 from testgenome1
  file(database_files)

  output:
  file 'genome1.par' into database_params

  """
  /usr/local/bin/last-train -P${task.cpus} ldb ${fasta1} > genome1.par
  """
}
*/

process ShuffleFasta {
  tag { qID }
  publishDir "outputs/stages/shuffled"

   queue 'dpetrov,normal,hns,owners'
  cpus 1
  memory { task.exitStatus == 137 ? (task.attempt > 2 ? 64.GB: 32.GB) : 12.GB}
  time { 2.h }
  errorStrategy { 'retry' }
  maxRetries 5
  maxErrors '-1'

  input:
  set val(qID), file(fasta) from query_fasta

  output:
  set qID, file("*.shuf.fasta") into shuffled_query

  """
  /usr/local/bin/seqkit shuffle ${fasta} > ${qID}.shuf.fasta
  """
}

shuffled_query.splitFasta( by: 50, file: true, elem: 1).set{split_query}

split_query.into{split_query_view; split_query}

//split_query_view.view()

process LastAlign {
  tag { qID }
  publishDir "outputs/stages/alignments"

  cpus 2
  memory { task.exitStatus == 137 ? (task.attempt > 2 ? 64.GB: 32.GB) : 12.GB}
  time { task.attempt == 1 ? 1.h: 24.h }
  errorStrategy { 'retry' }
  maxRetries 20
  maxErrors '-1'
 
  input:
  set val(qID), file(fasta_i) from split_query
  file db_files from database_files.first()
 
  output:
  set qID, file("*.maf") into aligned_mafs
 
  """
  /usr/local/bin/lastal -v -P${task.cpus} ${params.lastal_options} refdb ${fasta_i} | /usr/local/bin/last-split -v -m1 > ${fasta_i}.maf
  """
}

aligned_mafs.groupTuple().into{aligned_mafs_view; aligned_mafs}

aligned_mafs_view.view()

process MergeMAF {
  tag { qID }
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { task.exitStatus == 137 ? (task.attempt > 2 ? 64.GB: 32.GB) : 12.GB}
  time { 1.h}
  errorStrategy { 'retry' }
  maxRetries 5
  maxErrors '-1'

  input:
  set val(qID), file(mafs) from aligned_mafs

  output:
  set qID, file("*.maf") into aligned_maf

  script:
  header_maf = mafs[0]
  input_mafs = mafs.collect{"$it"}.join(' ')

  """
  grep -h "^#.*" $header_maf > ${qID}.maf
  grep -h -v "^#.*" $input_mafs >> ${qID}.maf
  """
}

aligned_maf.into{to_sam; to_blasttab; to_tab}

process MakeSamBam {
  tag { qID }
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { task.exitStatus == 137 ? (task.attempt > 2 ? 64.GB: 32.GB) : 12.GB}
  time { 2.h }
  errorStrategy { 'retry' }
  maxRetries 5
  maxErrors '-1'

  input:
  set val(qID), file(maf_file) from to_sam
  set val(prefix), file(ref_fasta) from ref_genome3.first()

  output:
  set qID, file("${qID}.sam"), file("${qID}.bam"), file("${qID}.bam.bai") into aligned_sam_bam
  file("${qID}.bam") into aligned_bams
  file("${qID}.bam.bai") into aligned_bais

  """
  set -e
  /usr/local/bin/maf-convert -n sam ${maf_file} > temp.sam
  /usr/local/bin/samtools faidx ${ref_fasta}
  /usr/local/bin/samtools view -t ${ref_fasta}.fai temp.sam > ${qID}.sam
  /usr/local/bin/samtools view -bu -t ${ref_fasta}.fai -T ${ref_fasta} ${qID}.sam | samtools addreplacerg -r ID:${qID} -r LB:${qID} -r SM:${qID} -o temp.bam
  /usr/local/bin/samtools sort temp.bam > ${qID}.bam
  /usr/local/bin/samtools index ${qID}.bam
  rm temp.sam temp.bam
  """
}

aligned_bams.into{ aligned_bams, aligned_bams_for_stats }

process MergeBams {
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { task.exitStatus == 137 ? (task.attempt > 2 ? 64.GB: 32.GB) : 12.GB}
  time { 2.h }
  errorStrategy { 'retry' }
  maxRetries 5
  maxErrors '-1'

  input:
  file(bams) from aligned_bams.toList()
  file(bais) from aligned_bais.toList()
  

  output:
  set file("*.bam"), file("*.bam.bai") into merged_bam

  """
  set -e
  /usr/local/bin/samtools merge -r all_pseudohaps.bam $bams
  /usr/local/bin/samtools index all_pseudohaps.bam
  
  """
}

process MakeBlasttab {
  tag { qID }
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { task.exitStatus == 137 ? (task.attempt > 2 ? 64.GB: 32.GB) : 12.GB}
  time { 2.h }
  errorStrategy { 'retry' }
  maxRetries 5
  maxErrors '-1'

  input:
  set val(qID), file(maf_file) from to_blasttab

  output:
  set qID, file("*.blasttab") into aligned_blasttab

  """
  /usr/local/bin/maf-convert -n blasttab ${maf_file} > ${qID}.blasttab
  """
}

process MakeTab {
  tag { qID }
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { task.exitStatus == 137 ? (task.attempt > 2 ? 64.GB: 32.GB) : 12.GB}
  time { 2.h }
  errorStrategy { 'retry' }
  maxRetries 5
  maxErrors '-1'

  input:
  set val(qID), file(maf_file) from to_tab

  output:
  set qID, file("*.tab") into aligned_tab

  """
  /usr/local/bin/maf-convert -n tab ${maf_file} > ${qID}.tab
  """
}

aligned_bams_for_stats = aligned_bams_for_stats.mix(merged_bam)

process Stats {

    tag { qID }
    publishDir "outputs", mode: 'copy'

    cpus 1
    memory { 12.GB}
    time { 2.h }
    errorStrategy { 'retry' }
    maxRetries 5
    maxErrors '-1'

    input:
    file(bam) from aligned_bams_for_stats

    output:
    file("*.stats") into bam_stats

    """
    /usr/bin/bamtools stats -in $bam > $bam.stats
    """

}
