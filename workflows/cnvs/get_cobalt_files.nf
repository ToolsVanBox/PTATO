include { CountBamLinesApplication } from '../../NextflowModules/gridss-purple-linx/1.3.2/CountBamLinesApplication.nf' params(params)

workflow get_cobalt_files {
  take:
    normal_bams
    tumor_bams
    genome_fasta
    genome_fai
    genome_dict
  main:
    input_cobalt = normal_bams
      .combine( tumor_bams, by: [0] )

    def gc_profile = file( params.cobalt.gc_profile, checkIfExists: true )
      
    CountBamLinesApplication( input_cobalt, genome_fasta, genome_fai, genome_dict, gc_profile )

    cobalt_ratio_tsvs = CountBamLinesApplication.out.cobalt_ratio_file
    .map{
        donor_id, normal_sample_id, tumor_sample_id, cobalt_ratio_tsv ->
        tsv_name = cobalt_ratio_tsv.getName()
        cobalt_ratio_tsv = cobalt_ratio_tsv.copyTo("${params.out_dir}/intermediate/cnvs/cobalt/${donor_id}/${normal_sample_id}/${tsv_name}")
            [ donor_id, normal_sample_id, tumor_sample_id, cobalt_ratio_tsv ]	
    }

    cobalt_files = CountBamLinesApplication.out.cobalt_files
      .transpose()
      .map{
        donor_id, normal_sample_id, tumor_sample_id, cobalt_file ->
        file_name = cobalt_file.getName()
        cobalt_file.copyTo("${params.out_dir}/intermediate/cnvs/cobalt/${donor_id}/${normal_sample_id}/${tumor_sample_id}/${file_name}")
        [ donor_id, normal_sample_id, tumor_sample_id, cobalt_file ]
      }
    
//    cobalt_ratio_tsvs = CountBamLinesApplication.out
//      .transpose()
//      .map{
//        donor_id, normal_sample_id, tumor_sample_id, cobalt_file, cobalt_ratio_tsv ->
//        file_name = cobalt_file.getName()
//        tsv_name = cobalt_ratio_tsv.getName()
//        cobalt_file.copyTo("${params.out_dir}/intermediate/cnvs/cobalt/${donor_id}/${normal_sample_id}/${tumor_sample_id}/${file_name}")
//        cobalt_ratio_tsv = cobalt_ratio_tsv.copyTo("${params.out_dir}/intermediate/cnvs/cobalt/${donor_id}/${normal_sample_id}/${tsv_name}")
//        [ donor_id, normal_sample_id, tumor_sample_id, cobalt_ratio_tsv ]
//      }
//      .unique()


  emit:
    cobalt_ratio_tsvs
}
