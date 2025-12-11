process AnnotateInsertedSequence {
  tag {"AnnotateInsertedSequence ${normal_sample_id} ${tumor_sample_id}"}
  label 'gridss'
  label 'gridss_2_13_2_AnnotateInsertedSequence'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"

    input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(gridss_driver_vcf), path(gridss_driver_tbi) )
    tuple( val(meta), path(viralreference) )
    path(repeatmaskerbed)

    output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.gridss.unfiltered.vcf.gz"), path("${tumor_sample_id}.gridss.unfiltered.vcf.gz.tbi"), emit: gridss_unfiltered_vcf )


    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${tumor_sample_id}"
    def args = task.ext.args ?: ''
    def repeatmaskerbed = repeatmaskerbed ? "REPEAT_MASKER_BED=${repeatmaskerbed}" : ""
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """ 
    FASTA=`find -L ./ -name "*.fa"`

    ls ./
    
    AnnotateInsertedSequence \\
        ALIGNMENT=APPEND \\
        REFERENCE_SEQUENCE=\$FASTA \\
	INPUT=${gridss_driver_vcf} \\
	OUTPUT=${prefix}.gridss.unfiltered.vcf.gz \\
        WORKER_THREADS=${task.cpus} \\
        ${repeatmaskerbed}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss_annotateinsertedsequence: ${VERSION}
    END_VERSIONS
    """    

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.gridss.unfiltered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss_annotateinsertedsequence: ${VERSION}
    END_VERSIONS
    """
}
