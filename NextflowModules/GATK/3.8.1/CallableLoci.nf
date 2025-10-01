
process CallableLoci {
  tag {"GATK CallableLoci ${sample_id}"}
  label 'GATK_3_8_1'
  label 'GATK_3_8_1_CallableLoci'
  clusterOptions = workflow.profile == "sge" ? "-l h_vmem=${params.mem}" : ""
  shell = ['/bin/bash', '-euo', 'pipefail']
  //container = 'quay.io/biocontainers/gatk:3.8--hdfd78af_11'
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/gatk:3.8--hdfd78af_11':
        'biocontainers/gatk:3.8--hdfd78af_11' }"

  input:
    tuple val(donor_id), val(sample_id), path(bam), path(bai) 
    tuple val(fasta_meta), path(genome_fasta)
    tuple val(fai_meta), path(genome_fai)
    tuple val(dict_meta), path(genome_dict)
    
  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.callableloci.bed"), path("${sample_id}.callableloci.txt"), emit: callableloci_files )

  script:
    """

    java -Djava.io.tmpdir=\$TMPDIR \
    -Xmx${task.memory.toGiga()-4}g -jar \
    /usr/local/opt/gatk-3.8/GenomeAnalysisTK.jar \
    -T CallableLoci \
    -I ${bam} \
    -R ${genome_fasta} \
    -o ${sample_id}.callableloci.bed \
    --summary ${sample_id}.callableloci.txt \
    ${params.callableloci.optional}
    """

    // """
    // gatk3 --java-options "-Xmx${task.memory.toGiga()-4}g -Djava.io.tmpdir=\$TMPDIR" \
    // CallableLoci \
    // -I ${bam} \
    // -R ${params.genome_fasta} \
    // -o ${sample_id}.callableloci.bed \
    // --summary ${sample_id}.callableloci.txt \
    // ${params.callableloci.optional}
    // """
}
