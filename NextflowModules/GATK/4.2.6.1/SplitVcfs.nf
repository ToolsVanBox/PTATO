
process SplitVcfs {
  tag {"GATK SplitVcfs ${sample_id}"}
  label 'GATK_4_2_6_1'
  label 'GATK_4_2_6_1_SplitVcfs'
  clusterOptions = workflow.profile == "sge" ? "-l h_vmem=${params.mem}" : ""
  shell = ['/bin/bash', '-euo', 'pipefail']
  //container = 'docker://broadinstitute/gatk:4.2.6.1'
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://broadinstitute/gatk:4.2.6.1':
        params.artifact_registry_path + '/gatk4:4.4.0.0--py36hdfd78af_0' }"

//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/gatk4@sha256:f986665e53f97b2726d96a8fb78234afb8a2897582818660b7b84b65a973f5d2' }"
        
  input:
    tuple( val(donor_id), val(sample_id), path(vcf), path(tbi) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.snvs.vcf.gz"), path("${sample_id}.snvs.vcf.gz.tbi"), path("${sample_id}.indels.vcf.gz"), path("${sample_id}.indels.vcf.gz.tbi"), emit: split_vcfs )

  script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()-4}g -Djava.io.tmpdir=\$TMPDIR" \
    SplitVcfs \
    -I ${vcf} \
    -SNP_OUTPUT ${sample_id}.snvs.vcf.gz \
    -INDEL_OUTPUT ${sample_id}.indels.vcf.gz \
    -STRICT false \
    ${params.splitvcfs.optional}
    """
}
