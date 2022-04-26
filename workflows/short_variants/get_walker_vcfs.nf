include { walker } from '../../NextflowModules/walker/2.1.2/walker.nf' params(params)
include { bgzip } from '../../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../../NextflowModules/htslib/1.15/tabix.nf' params(params)


workflow get_walker_vcfs {
  take:
    input_somatic_vcfs
    input_germline_vcfs
    input_bams
  main:
    input_walker = input_germline_vcfs
      .join( input_bams, by: [0] )
      .combine( input_somatic_vcfs, by: [0] )

    walker( input_walker )
    bgzip( walker.out
      .map{ donor_id, sample_id, walker_vcf, walker_bed, walker_txt ->
        bed_name = walker_bed.getName()
        txt_name = walker_txt.getName()
        walker_bed.copyTo("${params.out_dir}/intermediate/short_variants/walker/${donor_id}/${bed_name}")
        walker_txt.copyTo("${params.out_dir}/intermediate/short_variants/walker/${donor_id}/${txt_name}")
        [ donor_id, sample_id, walker_vcf ]
      }
    )
    tabix( bgzip.out )

    walker_vcfs = tabix.out
      .map{ donor_id, sample_id, vcf_gz, vcf_tbi ->
        vcf_name = vcf_gz.getName()
        tbi_name = vcf_tbi.getName()
        vcf_gz = vcf_gz.copyTo("${params.out_dir}/intermediate/short_variants/walker/${donor_id}/${vcf_name}")
        vcf_tbi = vcf_tbi.copyTo("${params.out_dir}/intermediate/short_variants/walker/${donor_id}/${tbi_name}")
        [ donor_id, sample_id, vcf_gz, vcf_tbi ]
      }
  emit:
    walker_vcfs
}
