
process CountBamLinesApplication {

  tag {"CountBamLinesApplication ${normal_sample_id} ${tumor_sample_id}"}
  label 'CountBamLinesApplication'
  label 'gripss_purple_linx_1_3_2_CountBamLinesApplication'
  //container = 'docker://gridss/gridss-purple-linx:1.3.2'
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://gridss/gridss-purple-linx:1.3.2':
        params.artifact_registry_path + '/gridss-purple-linx:1.3.3' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/gridss-purple-linx@sha256:803c157e841181c58af91fac6820317843b3fc053ae9ea99aa3a61e26c979d8c' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/gridss-purple-linx@sha256:ede83e9703a9a6ad667484d8a12b2fcf38393c9e1711fdba4cd17e40c3d617be' }"
  shell = ['/bin/bash', '-euo', 'pipefail']
        
  input:
    tuple( val(donor_id), val(normal_sample_id), path(normal_bam), path(normal_bai), val(tumor_sample_id), path(tumor_bam), path(tumor_bai) )
    tuple( val(fast_meta), path( genome_fasta ) )
    tuple( val(fai_meta), path( genome_fai ) )
    tuple( val(dict_meta), path( genome_dict ) )
    path( gc_profile )
  output:
//    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}/*"), path("${tumor_sample_id}/${tumor_sample_id}.cobalt.ratio.tsv"), emit: cobalt_files )
    tuple val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}/*"), emit: cobalt_files 
    tuple val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}/${tumor_sample_id}.cobalt.ratio.tsv"), emit: cobalt_ratio_file

  script:
    """
    java -Xmx${task.memory.toGiga()-4}g \
    -Dsamjdk.reference_fasta=${genome_fasta} \
    -Dsamjdk.use_async_io_read_samtools=true \
    -Dsamjdk.use_async_io_write_samtools=true \
    -Dsamjdk.use_async_io_write_tribble=true \
    -Dsamjdk.buffer_size=4194304 \
    -Dsamjdk.async_io_read_threads=${task.cpus} \
		-cp /opt/hmftools/cobalt-1.11.jar com.hartwig.hmftools.cobalt.CountBamLinesApplication \
		-threads ${task.cpus} \
		-tumor ${tumor_sample_id} \
		-tumor_bam ${tumor_bam} \
		-ref_genome ${genome_fasta} \
		-output_dir ./${tumor_sample_id} \
		-gc_profile ${gc_profile} \
		-reference ${normal_sample_id} \
    -reference_bam ${normal_bam} \
    ${params.countbamlinesapplication.optional}
    """
}
