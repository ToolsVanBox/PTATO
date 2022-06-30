#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { bgzip } from '../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../NextflowModules/htslib/1.15/tabix.nf' params(params)

workflow get_gzipped_vcfs {
  take:
    input_vcfs
  main:
    input_unzipped_vcfs = input_vcfs.filter{ it[2] =~ /.vcf$/ }
    input_unindexed_vcfs = input_vcfs.filter{ it[2] =~ /.vcf.gz$/ }.filter{ it.size() == 3 }
    gzipped_vcf_files = input_vcfs.filter{ it[2] =~ /.vcf.gz$/ }.filter{ it[3] =~ /.tbi$/ }

    bgzip( input_unzipped_vcfs )
    input_unindexed_vcfs = input_unindexed_vcfs.concat( bgzip.out )

    tabix( input_unindexed_vcfs )
    gzipped_vcf_files = gzipped_vcf_files.concat( tabix.out )

  emit:
    gzipped_vcf_files
}
