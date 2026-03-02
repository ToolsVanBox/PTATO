
process gridss {
  tag {"gridss ${normal_sample_id} ${tumor_sample_id}"}
  label 'gridss'
  label 'gridss_2_13_2'
  shell = ['/bin/bash', '-euo', 'pipefail']
  //container = 'docker://gridss/gridss:2.13.2'
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://gridss/gridss:2.13.2':
        params.artifact_registry_path + '/gridss:2.13.2' }"

//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/gridss:2.13.2@sha256:3a0b53fb9c37891cf6429b24ca70fd77454fd3103966469b3c5124b55997453d' }"
//        'biocontainers/gridss:2.13.2--h270b39a_0' }"
        
  input:
    tuple( val(donor_id), val(normal_sample_id), path(normal_bam), path(normal_bai), val(tumor_sample_id), path(tumor_bam), path(tumor_bai) )
    tuple( val(fast_meta), path( genome_fasta ) )
    tuple( val(fai_meta), path( genome_fai ) )
    tuple( val(dict_meta), path( genome_dict ) )
  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.gridss.driver.vcf.gz"), path("${tumor_sample_id}.gridss.driver.vcf.gz.tbi"), path("${tumor_sample_id}.gridss.driver.vcf.gz.assembly.bam"), emit: gridss_driver_vcf )

  script:
    """
    gridss \
    --jvmheap ${task.memory.toGiga()-4}g \
    -o ${tumor_sample_id}.gridss.driver.vcf.gz \
    -r ${genome_fasta} \
    -t ${task.cpus} \
    --labels ${normal_sample_id},${tumor_sample_id} \
    ${params.gridss.optional} \
    ${normal_bam} \
    ${tumor_bam}
    """
}
