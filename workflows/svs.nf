#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
  extractGridssUnfilteredVcfGzFromDir;
  extractGripssSomaticFilteredVcfGzFromDir;
  extractCobaltRatioTsvFromDir;
} from '../NextflowModules/Utils/getFilesFromDir.nf' params(params)

include { get_gridss_vcfs } from './svs/get_gridss_vcfs.nf' params(params)
include { get_gripss_vcfs } from './svs/get_gripss_vcfs.nf' params(params)
include { get_cobalt_files } from './svs/get_cobalt_files.nf' params(params)
include { filter_sv_files } from './svs/filter_sv_files.nf' params(params)
// include { filter_cobalt_files } from './svs/filter_cobalt_files.nf' params(params)
// include { filter_baf_files } from './svs/filter_baf_files.nf' params(params)
// include { filter_gripss_vcfs } from './svs/filter_gripss_vcfs.nf' params(params)
// include { integrate_files } from './svs/integrate_files.nf' params(params)

workflow svs {
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

    if ( params.optional.svs.gripss_somatic_filtered_vcfs_dir ) {
      gripss_somatic_filtered_vcfs = extractGripssSomaticFilteredVcfGzFromDir( params.optional.svs.gripss_somatic_filtered_vcfs_dir )
    } else {
      if ( params.optional.svs.gridss_unfiltered_vcfs_dir ) {
        gridss_unfiltered_vcfs = extractGridssUnfilteredVcfGzFromDir( params.optional.svs.gridss_unfiltered_vcfs_dir )
      } else {
        get_gridss_vcfs( normal_bams, tumor_bams )
        gridss_unfiltered_vcfs = get_gridss_vcfs.out
      }
      get_gripss_vcfs( gridss_unfiltered_vcfs )
      gripss_somatic_filtered_vcfs = get_gripss_vcfs.out
    }

    if ( params.optional.svs.cobalt_ratio_tsv_dir ) {
      cobalt_ratio_tsv_files = extractCobaltRatioTsvFromDir( params.optional.svs.cobalt_ratio_tsv_dir )
    } else {
      get_cobalt_files( normal_bams, tumor_bams )
      cobalt_ratio_tsv_files = get_cobalt_files.out
    }

    filter_sv_files( cobalt_ratio_tsv_files, germline_vcfs, normal_bams, tumor_bams, gripss_somatic_filtered_vcfs )

}
