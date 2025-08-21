
process CollectWGSMetrics {
  tag {"GATK CollectWGSMetrics ${sample_id}"}
  label 'GATK_4_2_6_1'
  label 'GATK_4_2_6_1_CollectWGSMetrics'
  clusterOptions = workflow.profile == "sge" ? "-l h_vmem=${params.mem}" : ""
  shell = ['/bin/bash', '-euo', 'pipefail']
  //container = 'docker://broadinstitute/gatk:4.2.6.1'
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://broadinstitute/gatk:4.2.6.1':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/gatk4@sha256:f986665e53f97b2726d96a8fb78234afb8a2897582818660b7b84b65a973f5d2' }"
           
  input:
    tuple( val(donor_id), val(sample_id), path(bam), path(bai) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.wgs_metrics.txt"), emit: wgs_metrics )

  script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()-4}g -Djava.io.tmpdir=\$TMPDIR" \
    CollectWgsMetrics \
    -I ${bam} \
    -O ${sample_id}.wgs_metrics.txt \
    -R ${params.genome_fasta} \
    ${params.collectwgsmetrics.optional}
    sed -i 's/picard\\.analysis\\.WgsMetrics/picard\\.analysis\\.CollectWgsMetrics\\\$WgsMetrics/' ${sample_id}.wgs_metrics.txt
    """
}
