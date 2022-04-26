include { smurf } from '../../NextflowModules/SMuRF/2.1.5/SMuRF.nf' params(params)
include { bgzip } from '../../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../../NextflowModules/htslib/1.15/tabix.nf' params(params)


workflow get_somatic_vcfs {
  take:
    input_germline_vcfs
    input_bams
  main:
    input_SMuRF = input_germline_vcfs
      .join( input_bams, by: [0] )
      .join( Channel.from( params.bulk_names ).groupTuple( by: [0] ) )

    smurf( input_SMuRF )

    bgzip( smurf.out
      .transpose()
      .map{ donor_id, smurf_vcf, smurf_filtered_vcf, vaf_pdf, somatic_vcf ->
        smurf_name = smurf_vcf.getName()
        smurf_filtered_name = smurf_filtered_vcf.getName()
        vaf_pdf_name = vaf_pdf.getName()
        smurf_vcf.copyTo("${params.out_dir}/intermediate/short_variants/SMuRF/${donor_id}/${smurf_name}")
        smurf_filtered_vcf.copyTo("${params.out_dir}/intermediate/short_variants/SMuRF/${donor_id}/${smurf_filtered_name}")
        vaf_pdf.copyTo("${params.out_dir}/intermediate/short_variants/SMuRF/${donor_id}/${vaf_pdf_name}")
        sample_id = somatic_vcf.getName().toString().replaceAll(/.SMuRF.filtered.vcf$/, '')
        [ donor_id, sample_id, somatic_vcf ]
      }
    )
    tabix( bgzip.out )
    somatic_vcfs = tabix.out
      .map{ donor_id, sample_id, vcf_gz, vcf_tbi ->
        vcf_name = vcf_gz.getName()
        tbi_name = vcf_tbi.getName()
        vcf_gz = vcf_gz.copyTo("${params.out_dir}/intermediate/short_variants/SMuRF/${donor_id}/${vcf_name}")
        vcf_tbi = vcf_tbi.copyTo("${params.out_dir}/intermediate/short_variants/SMuRF/${donor_id}/${tbi_name}")
        [ donor_id, sample_id, vcf_gz, vcf_tbi ]
      }

  emit:
    somatic_vcfs
}
