include { ESVEE } from '../moduleNextflowModuless/esvee/main.nf' params(params)

workflow get_esvee_vcfs {
    take:
        normal_bams
        tumor_bams
        genome_fasta
        genome_fai
        genome_dict

    main:
        input_esvee = normal_bams
            .combine( tumor_bams, by: [0] )
        ESVEE( input_esvee, genome_fasta, genome_fai, genome_dict )
        ESVEE.out.view()
}