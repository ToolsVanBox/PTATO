#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
  extractGridssUnfilteredVcfFromDir;
  extractGripssSomaticFilteredVcfFromDir;
  extractCobaltRatioTsvFromDir;
} from '../NextflowModules/Utils/getFilesFromDir.nf' params(params)

include { get_gzipped_vcfs } from './get_gzipped_vcfs.nf' params(params)
include { get_gridss_vcfs } from './svs/get_gridss_vcfs.nf' params(params)
include { get_gripss_vcfs } from './svs/get_gripss_vcfs.nf' params(params)
include { combine_cnvs_svs } from './svs/combine_cnvs_svs.nf' params(params)


workflow svs {
  take:
    bams
    germline_vcfs
    filtered_cnv_files
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
      raw_gripss_somatic_filtered_vcfs = extractGripssSomaticFilteredVcfFromDir( params.optional.svs.gripss_somatic_filtered_vcfs_dir )
      get_gzipped_vcfs( raw_gripss_somatic_filtered_vcfs )
      gripss_somatic_filtered_vcfs = get_gzipped_vcfs.out
        .map{
          donor_id, sample_id, vcf, tbi ->
          m = sample_id =~ /(.+);(.+)/
          normal_sample_id = m[0][1]
          tumor_sample_id = m[0][2]
          [ donor_id, normal_sample_id, tumor_sample_id, vcf, tbi]
        }

    } else {
      if ( params.optional.svs.gridss_unfiltered_vcfs_dir ) {
        raw_gridss_unfiltered_vcfs = extractGridssUnfilteredVcfFromDir( params.optional.svs.gridss_unfiltered_vcfs_dir )
        get_gzipped_vcfs( raw_gridss_unfiltered_vcfs )
        gridss_unfiltered_vcfs = get_gzipped_vcfs.out
          .map{
            donor_id, sample_id, vcf, tbi ->
            m = sample_id =~ /(.+);(.+)/
            normal_sample_id = m[0][1]
            tumor_sample_id = m[0][2]
            [ donor_id, normal_sample_id, tumor_sample_id, vcf, tbi]
          }
      } else {
        get_gridss_vcfs( normal_bams, tumor_bams )
        gridss_unfiltered_vcfs = get_gridss_vcfs.out
      }
      get_gripss_vcfs( gridss_unfiltered_vcfs )
      gripss_somatic_filtered_vcfs = get_gripss_vcfs.out
    }

    combine_cnvs_svs( filtered_cnv_files, gripss_somatic_filtered_vcfs )
}
