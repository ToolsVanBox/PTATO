include { CountBamLinesApplication } from '../../NextflowModules/gridss-purple-linx/1.3.2/CountBamLinesApplication.nf' params(params)

workflow get_cobalt_files {
  take:
    normal_bams
    tumor_bams
  main:
    input_cobalt = normal_bams
      .combine( tumor_bams, by: [0] )
    CountBamLinesApplication( input_cobalt )

    cobalt_ratio_tsvs = CountBamLinesApplication.out
      .transpose()
      .map{
        donor_id, normal_sample_id, tumor_sample_id, cobalt_file, cobalt_ratio_tsv ->
        file_name = cobalt_file.getName()
        tsv_name = cobalt_ratio_tsv.getName()
        cobalt_file.copyTo("${params.out_dir}/intermediate/cnvs/cobalt/${donor_id}/${normal_sample_id}/${tumor_sample_id}/${file_name}")
        cobalt_ratio_tsv = cobalt_ratio_tsv.copyTo("${params.out_dir}/intermediate/cnvs/cobalt/${donor_id}/${normal_sample_id}/${tsv_name}")
        [ donor_id, normal_sample_id, tumor_sample_id, cobalt_ratio_tsv ]
      }
      .unique()
  emit:
    cobalt_ratio_tsvs
}
