#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { extractGermlineVcfGzFromDir } from '../NextflowModules/Utils/getFilesFromDir.nf' params(params)
include { SnpSift } from '../NextflowModules/SnpSift/4.3.1t--1/SnpSift.nf' params(params)
include { bgzip } from '../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../NextflowModules/htslib/1.15/tabix.nf' params(params)

workflow get_germline_vcfs {
  take:
    input_vcfs
  main:
    bulk_names = Channel.from( params.bulk_names ).groupTuple( by: [0] )

    if ( params.optional.germline_vcfs_dir ) {
      germline_vcfs = extractGermlineVcfGzFromDir( params.optional.germline_vcfs_dir )
    } else {
      input_snpsift = input_vcfs.join( bulk_names, by: [0] )

      SnpSift( input_snpsift )
      bgzip( SnpSift.out )
      tabix( bgzip.out )
      germline_vcfs = tabix.out
        .map{ donor_id, sample_id, vcf_gz, vcf_tbi ->
          vcf_name = vcf_gz.getName()
          tbi_name = vcf_tbi.getName()
          vcf_gz = vcf_gz.copyTo("${params.out_dir}/intermediate/germline/${donor_id}/${vcf_name}")
          vcf_tbi = vcf_tbi.copyTo("${params.out_dir}/intermediate/germline/${donor_id}/${tbi_name}")
          [ donor_id, sample_id, vcf_gz, vcf_tbi ]
        }
    }
  emit:
    germline_vcfs
}
