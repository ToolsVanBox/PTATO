
process AnnotateInsertedSequence {
  tag {"AnnotateInsertedSequence ${normal_sample_id} ${tumor_sample_id}"}
  label 'gridss'
  label 'gridss_2_13_2_AnnotateInsertedSequence'
  shell = ['/bin/bash', '-euo', 'pipefail']
  //container = 'docker://gridss/gridss:2.13.2'
//  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'docker://gridss/gridss:2.13.2':
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/gridss:2.13.2@sha256:3a0b53fb9c37891cf6429b24ca70fd77454fd3103966469b3c5124b55997453d' }"
//        'biocontainers/gridss:2.13.2--h270b39a_0' }"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"


  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(gridss_driver_vcf), path(gridss_driver_tbi) )
    path( viralreference )
    path( viral_fai )
  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.gridss.unfiltered.vcf.gz"), path("${tumor_sample_id}.gridss.unfiltered.vcf.gz.tbi"), emit: gridss_unfiltered_vcf )

  script:
//    """
//    java -Xmx${task.memory.toGiga()-4}g \
//    -Dsamjdk.create_index=true \
//		-Dsamjdk.use_async_io_read_samtools=true \
//		-Dsamjdk.use_async_io_write_samtools=true \
//		-Dsamjdk.use_async_io_write_tribble=true \
//		-Dsamjdk.buffer_size=4194304 \
//		-cp /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar gridss.AnnotateInsertedSequence \
//		REFERENCE_SEQUENCE=${viralreference} \
//		INPUT=${gridss_driver_vcf} \
//		OUTPUT=${tumor_sample_id}.gridss.unfiltered.vcf.gz \
//		ALIGNMENT=APPEND WORKER_THREADS=${task.cpus}
//    """
    """
    AnnotateInsertedSequence \\
        REFERENCE_SEQUENCE=${viralreference} \
		INPUT=${gridss_driver_vcf} \
		OUTPUT=${tumor_sample_id}.gridss.unfiltered.vcf.gz \
        WORKER_THREADS=${task.cpus} \

    """    
//    cat <<-END_VERSIONS > versions.yml
//    "${task.process}":
//        gridss_annotateinsertedsequence: ${VERSION}
//    END_VERSIONS
//    """    

}
