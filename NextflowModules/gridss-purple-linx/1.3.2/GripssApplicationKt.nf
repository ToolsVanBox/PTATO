
process GripssApplicationKt {

  tag {"GripssApplicationKt ${normal_sample_id} ${tumor_sample_id}"}
  label 'GripssApplicationKt'
  label 'gripss_purple_linx_1_3_2_GripssApplicationKt'
  shell = ['/bin/bash', '-euo', 'pipefail']
  //container = 'docker://gridss/gridss-purple-linx:1.3.2'
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://gridss/gridss-purple-linx:1.3.2':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/gridss-purple-linx@sha256:ede83e9703a9a6ad667484d8a12b2fcf38393c9e1711fdba4cd17e40c3d617be' }"
        
  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(gridss_unfiltered_vcf), path(gridss_unfiltered_tbi) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.gripss.somatic.vcf.gz"), path("${tumor_sample_id}.gripss.somatic.vcf.gz.tbi"), emit: gripss_somatic_vcf )

  script:
    """
    java -Xmx${task.memory.toGiga()-8}g -cp /opt/hmftools/gripss-1.9.jar com.hartwig.hmftools.gripss.GripssApplicationKt \
    -ref_genome ${params.genome_fasta} \
    -input_vcf ${gridss_unfiltered_vcf} \
    -output_vcf ${tumor_sample_id}.gripss.somatic.vcf.gz \
    -tumor ${tumor_sample_id} \
    -reference ${normal_sample_id} \
    -breakpoint_hotspot ${params.gripss.breakpoint_hotspot} \
    -breakend_pon ${params.gripss.breakend_pod} \
    -breakpoint_pon ${params.gripss.breakpoint_pon} \
    ${params.gripssapplicationkt.optional}
    """
}
