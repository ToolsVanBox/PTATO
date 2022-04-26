
process AnnotateInsertedSequence {
  tag {"AnnotateInsertedSequence ${normal_sample_id} ${tumor_sample_id}"}
  label 'gridss'
  label 'gridss_2_13_2_AnnotateInsertedSequence'
  container = 'docker://gridss/gridss:2.13.2'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(gridss_driver_vcf), path(gridss_driver_tbi) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.gridss.unfiltered.vcf.gz"), path("${tumor_sample_id}.gridss.unfiltered.vcf.gz.tbi"), emit: gridss_unfiltered_vcf )

  script:
    """
    java -Xmx${task.memory.toGiga()-4}g \
    -Dsamjdk.create_index=true \
		-Dsamjdk.use_async_io_read_samtools=true \
		-Dsamjdk.use_async_io_write_samtools=true \
		-Dsamjdk.use_async_io_write_tribble=true \
		-Dsamjdk.buffer_size=4194304 \
		-cp /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar gridss.AnnotateInsertedSequence \
		REFERENCE_SEQUENCE=${params.gridss.viralreference} \
		INPUT=${gridss_driver_vcf} \
		OUTPUT=${tumor_sample_id}.gridss.unfiltered.vcf.gz \
		ALIGNMENT=APPEND WORKER_THREADS=${task.cpus}
    """
}
