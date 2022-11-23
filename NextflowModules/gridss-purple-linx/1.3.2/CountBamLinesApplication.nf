
process CountBamLinesApplication {

  tag {"CountBamLinesApplication ${normal_sample_id} ${tumor_sample_id}"}
  label 'CountBamLinesApplication'
  label 'gripss_purple_linx_1_3_2_CountBamLinesApplication'
  container = 'docker://gridss/gridss-purple-linx:1.3.2'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(normal_sample_id), path(normal_bam), path(normal_bai), val(tumor_sample_id), path(tumor_bam), path(tumor_bai) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}/*"), path("${tumor_sample_id}/${tumor_sample_id}.cobalt.ratio.tsv"), emit: cobalt_files )

  script:
    """
    java -Xmx${task.memory.toGiga()-4}g \
    -Dsamjdk.reference_fasta=${params.genome_fasta} \
    -Dsamjdk.use_async_io_read_samtools=true \
    -Dsamjdk.use_async_io_write_samtools=true \
    -Dsamjdk.use_async_io_write_tribble=true \
    -Dsamjdk.buffer_size=4194304 \
    -Dsamjdk.async_io_read_threads=${task.cpus} \
		-cp /opt/hmftools/cobalt-1.11.jar com.hartwig.hmftools.cobalt.CountBamLinesApplication \
		-threads ${task.cpus} \
		-tumor ${tumor_sample_id} \
		-tumor_bam ${tumor_bam} \
		-ref_genome ${params.genome_fasta} \
		-output_dir ./${tumor_sample_id} \
		-gc_profile ${params.cobalt.gc_profile} \
		-reference ${normal_sample_id} \
    -reference_bam ${normal_bam} \
    ${params.countbamlinesapplication.optional}
		"""
}
