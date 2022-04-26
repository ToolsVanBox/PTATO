
process GetSampleName {
  tag {"GATK GetSampleName ${sample_id}"}
  label 'GATK_4_1_3_0'
  label 'GATK_4_1_3_0_GetSampleName'
  clusterOptions = workflow.profile == "sge" ? "-l h_vmem=${params.mem}" : ""
  container = 'library://sawibo/default/bioinf-tools:gatk4.1.3.0'
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
