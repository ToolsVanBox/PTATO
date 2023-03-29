
process SplitVcfs {
  tag {"GATK SplitVcfs ${sample_id}"}
  label 'GATK_4_2_6_1'
  label 'GATK_4_2_6_1_SplitVcfs'
  clusterOptions = workflow.profile == "sge" ? "-l h_vmem=${params.mem}" : ""
  container = 'docker://broadinstitute/gatk:4.2.6.1'
  shell = ['/bin/bash', '-euo', 'pipefail']

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
