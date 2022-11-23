#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
  extractGridssUnfilteredVcfFromDir;
  extractGripssSomaticFilteredVcfFromDir;
  extractCobaltRatioTsvFromDir;
} from '../NextflowModules/Utils/getFilesFromDir.nf' params(params)

include { get_gzipped_vcfs } from './get_gzipped_vcfs.nf' params(params)
include { get_cobalt_files } from './cnvs/get_cobalt_files.nf' params(params)
include { filter_cnv_files } from './cnvs/filter_cnv_files.nf' params(params)

workflow cnvs {
  take:
    bams
    germline_vcfs
  main:
    bulk_names = Channel.from( params.bulk_names ).combine(Channel.from("Normal"))

    normal_bams = bulk_names
      .join( bams.transpose(), by: [0,1], remainder: true )
      .filter{ it[2] == 'Normal' }
      .map{
        donor_id, sample_id, label, bam, bai ->
        [donor_id, sample_id, bam, bai]
      }

    tumor_bams = bulk_names
        .join( bams.transpose(), by: [0,1], remainder: true )
        .filter{ it[2] == null }
        .map{
          donor_id, sample_id, label, bam, bai ->
          [donor_id, sample_id, bam, bai]
        }

    if ( params.optional.cnvs.cobalt_ratio_tsv_dir ) {
      cobalt_ratio_tsv_files = extractCobaltRatioTsvFromDir( params.optional.cnvs.cobalt_ratio_tsv_dir )
    } else {
      get_cobalt_files( normal_bams, tumor_bams )
      cobalt_ratio_tsv_files = get_cobalt_files.out
    }
  //
    filter_cnv_files( cobalt_ratio_tsv_files, germline_vcfs, normal_bams, tumor_bams )
    filtered_cnv_files = filter_cnv_files.out

  emit:
    filtered_cnv_files
}
