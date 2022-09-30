#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { Index } from '../NextflowModules/Sambamba/0.8.2/Index.nf' params(params)
include { GetSampleName } from '../NextflowModules/GATK/4.2.6.1/GetSampleName.nf' params(params)

workflow get_indexed_bams {
  take:
    input_bams
  main:
    input_unindexed_bams = input_bams.filter{ it.size() == 3 }
    indexed_bams = input_bams
      .filter{ it[2] =~ /.bam$/ }
      .filter{ it[3] =~ /.bai$/ }

    Index( input_unindexed_bams )

    indexed_bams = indexed_bams
      .concat( Index.out )

    GetSampleName( indexed_bams )
    indexed_bams = GetSampleName.out

  emit:
    indexed_bams
}
