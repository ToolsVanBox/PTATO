
process GripssHardFilterApplicationKt {

  tag {"GripssHardFilterApplicationKt ${normal_sample_id} ${tumor_sample_id}"}
  label 'GripssHardFilterApplicationKt'
  label 'gripss_purple_linx_1_3_2_GripssHardFilterApplicationKt'
  shell = ['/bin/bash', '-euo', 'pipefail']
  //container = 'docker://gridss/gridss-purple-linx:1.3.2'
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://gridss/gridss-purple-linx:1.3.2':
        params.artifact_registry_path + '/gridss-purple-linx:1.3.3' }"

//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/gridss-purple-linx@sha256:803c157e841181c58af91fac6820317843b3fc053ae9ea99aa3a61e26c979d8c' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/gridss-purple-linx@sha256:ede83e9703a9a6ad667484d8a12b2fcf38393c9e1711fdba4cd17e40c3d617be' }"
        
  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(gripss_somatic_vcf), path(gripss_somatic_tbi) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.gripss.somatic.filtered.vcf.gz"), path("${tumor_sample_id}.gripss.somatic.filtered.vcf.gz.tbi"), emit: gripss_somatic_filtered_vcf )

  script:
    """
    java -Xmx${task.memory.toGiga()-4}g -cp /opt/hmftools/gripss-1.9.jar  com.hartwig.hmftools.gripss.GripssHardFilterApplicationKt \
    -input_vcf ${gripss_somatic_vcf} \
    -output_vcf ${tumor_sample_id}.gripss.somatic.filtered.vcf.gz \
    ${params.gripsshardfilterapplicationkt.optional}
    """
}
