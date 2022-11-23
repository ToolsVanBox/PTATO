
process gridss {
  tag {"gridss ${normal_sample_id} ${tumor_sample_id}"}
  label 'gridss'
  label 'gridss_2_13_2'
  container = 'docker://gridss/gridss:2.13.2'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(normal_sample_id), path(normal_bam), path(normal_bai), val(tumor_sample_id), path(tumor_bam), path(tumor_bai) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.gridss.driver.vcf.gz"), path("${tumor_sample_id}.gridss.driver.vcf.gz.tbi"), path("${tumor_sample_id}.gridss.driver.vcf.gz.assembly.bam"), emit: gridss_driver_vcf )

  script:
    """
    gridss \
    --jvmheap ${task.memory.toGiga()-4}g \
    -o ${tumor_sample_id}.gridss.driver.vcf.gz \
    -r ${params.genome_fasta} \
    -t ${task.cpus} \
    --labels ${normal_sample_id},${tumor_sample_id} \
    ${params.gridss.optional} \
    ${normal_bam} \
    ${tumor_bam}
    """
}
