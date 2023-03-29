
process GetSampleName {
  tag {"GATK GetSampleName ${sample_id}"}
  label 'GATK_4_2_6_1'
  label 'GATK_4_2_6_1_GetSampleName'
  clusterOptions = workflow.profile == "sge" ? "-l h_vmem=${params.mem}" : ""
  container = 'docker://broadinstitute/gatk:4.2.6.1'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(bam), path(bai) )

  output:
    tuple( val(donor_id), stdout, path(bam), path(bai), emit: samplename_bams )

  script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()-4}g -Djava.io.tmpdir=\$TMPDIR" \
    GetSampleName \
    -I ${bam} \
    -O /dev/stdout
    ${params.getsamplename.optional}
    """
}
