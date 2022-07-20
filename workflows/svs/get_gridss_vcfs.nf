include { gridss } from '../../NextflowModules/gridss/2.13.2/gridss.nf' params(params)
include { AnnotateInsertedSequence } from '../../NextflowModules/gridss/2.13.2/AnnotateInsertedSequence.nf' params(params)
include {
    extractGridssDriverVcfFromDir
} from '../../NextflowModules/Utils/getFilesFromDir.nf' params(params)

workflow get_gridss_vcfs {
  take:
    normal_bams
    tumor_bams
  main:
    input_gridss = normal_bams
      .combine( tumor_bams, by: [0] )
    if ( params.optional.svs.gridss_driver_vcfs_dir ) {
      raw_gridss_driver_vcfs = extractGridssDriverVcfFromDir( params.optional.svs.gridss_driver_vcfs_dir )
      get_gzipped_vcfs( raw_gridss_driver_vcfs )
      gridss_driver_vcfs = get_gzipped_vcfs.out
        .map{
          donor_id, sample_id, vcf, tbi ->
          m = sample_id =~ /(.+);(.+)/
          normal_sample_id = m[0][1]
          tumor_sample_id = m[0][2]
          [ donor_id, normal_sample_id, tumor_sample_id, vcf, tbi]
        }
    } else {
      gridss( input_gridss )
      gridss_driver_vcfs = gridss.out
        .map{
          donor_id, normal_sample_id, tumor_sample_id, gridss_driver_vcf, gridss_driver_tbi, gridss_driver_bam ->
          vcf_name = gridss_driver_vcf.getName()
          tbi_name = gridss_driver_tbi.getName()
          bam_name = gridss_driver_bam.getName()
          gridss_driver_vcf = gridss_driver_vcf.copyTo("${params.out_dir}/intermediate/svs/gridss/${donor_id}/${normal_sample_id}/${vcf_name}")
          gridss_driver_tbi = gridss_driver_tbi.copyTo("${params.out_dir}/intermediate/svs/gridss/${donor_id}/${normal_sample_id}/${tbi_name}")
          gridss_driver_bam.copyTo("${params.out_dir}/intermediate/svs/gridss/${donor_id}/${normal_sample_id}/${bam_name}")
          [ donor_id, normal_sample_id, tumor_sample_id, gridss_driver_vcf, gridss_driver_tbi ]
        }
    }
    AnnotateInsertedSequence( gridss_driver_vcfs )
    gridss_unfiltered_vcfs = AnnotateInsertedSequence.out
      .map{
        donor_id, normal_sample_id, tumor_sample_id, gridss_unfiltered_vcf, gridss_unfiltered_tbi ->
        vcf_name = gridss_unfiltered_vcf.getName()
        tbi_name = gridss_unfiltered_tbi.getName()
        gridss_unfiltered_vcf = gridss_unfiltered_vcf.copyTo("${params.out_dir}/intermediate/svs/gridss/${donor_id}/${normal_sample_id}/${vcf_name}")
        gridss_unfiltered_tbi = gridss_unfiltered_tbi.copyTo("${params.out_dir}/intermediate/svs/gridss/${donor_id}/${normal_sample_id}/${tbi_name}")
        [ donor_id, normal_sample_id, tumor_sample_id, gridss_unfiltered_vcf, gridss_unfiltered_tbi ]
      }
  emit:
    gridss_unfiltered_vcfs
}
