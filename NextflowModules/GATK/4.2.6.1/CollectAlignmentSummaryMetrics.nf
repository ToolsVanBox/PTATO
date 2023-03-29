
process CollectAlignmentSummaryMetrics {
  tag {"GATK CollectAlignmentSummaryMetrics ${sample_id}"}
  label 'GATK_4_2_6_1'
  label 'GATK_4_2_6_1_CollectAlignmentSummaryMetrics'
  clusterOptions = workflow.profile == "sge" ? "-l h_vmem=${params.mem}" : ""
  container = 'docker://broadinstitute/gatk:4.2.6.1'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(bam), path(bai) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.alignment_summary_metrics.txt"), emit: alignment_summary_metrics )

  script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()-4}g -Djava.io.tmpdir=\$TMPDIR" \
    CollectAlignmentSummaryMetrics \
    -I ${bam} \
    -O ${sample_id}.alignment_summary_metrics.txt \
    -R ${params.genome_fasta} \
    ${params.collectalignmentsummarymetrics.optional}
    sed -i 's/picard\\.analysis\\.AlignmentSummaryMetrics/picard\\.analysis\\.CollectAlignmentSummaryMetrics\\\$AlignmentSummaryMetrics/' ${sample_id}.alignment_summary_metrics.txt
    """
}
