#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { get_gzipped_vcfs } from './get_gzipped_vcfs.nf' params(params)
include { get_gzipped_vcfs as get_gzipped_vcfs2 } from './get_gzipped_vcfs.nf' params(params)

include {
  extractSomaticVcfFromDir;
  extractWalkerVcfFromDir;
  extractAbTableFromDir;
  extractContextBedFromDir;
  extractFeaturesBedFromDir;
  extractPtaVcfFromDir;
  extractNoptaVcfFromDir;
} from '../NextflowModules/Utils/getFilesFromDir.nf' params(params)

include { get_ab_tables } from './short_variants/get_ab_tables.nf' params(params)
include { get_walker_vcfs } from './short_variants/get_walker_vcfs.nf' params(params)
include { get_somatic_vcfs } from './short_variants/get_somatic_vcfs.nf' params(params)
include { get_context_beds } from './short_variants/get_context_beds.nf' params(params)
include { closest_feature } from './short_variants/get_features_beds.nf' params(params)
include { intersect_feature } from './short_variants/get_features_beds.nf' params(params)
include { merge_features } from './short_variants/get_features_beds.nf' params(params)

include { SplitVcfs } from '../NextflowModules/GATK/4.2.6.1/SplitVcfs.nf' params(params)
include { bgzip } from '../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../NextflowModules/htslib/1.15/tabix.nf' params(params)

include { intersect_ptato_vcfs } from './short_variants/intersect_ptato_vcfs.nf' params(params)
include { merge_ptato_vcfs } from './short_variants/merge_ptato_vcfs.nf' params(params)

include {
  snvs;
  snvs_train;
} from './snvs.nf' params(params)

include {
  indels;
  indels_train;
} from './indels.nf' params(params)


workflow short_variants {
  take:
    input_vcfs
    bams
    germline_vcfs
  main:
    if ( params.optional.short_variants.somatic_vcfs_dir ) {
      raw_somatic_vcfs = extractSomaticVcfFromDir( params.optional.short_variants.somatic_vcfs_dir )
      get_gzipped_vcfs( raw_somatic_vcfs )
      somatic_vcfs = get_gzipped_vcfs.out
    } else {
      get_somatic_vcfs( input_vcfs, bams )
      somatic_vcfs = get_somatic_vcfs.out
    }

    if ( params.optional.short_variants.ab_tables_dir ) {
      ab_tables = extractAbTableFromDir( params.optional.short_variants.ab_tables_dir )
    } else {
      get_ab_tables( germline_vcfs, somatic_vcfs )
      ab_tables = get_ab_tables.out
    }

    if ( params.optional.short_variants.walker_vcfs_dir ) {
      raw_walker_vcfs = extractWalkerVcfFromDir( params.optional.short_variants.walker_vcfs_dir )
      get_gzipped_vcfs2( raw_walker_vcfs )
      walker_vcfs = get_gzipped_vcfs2.out
    } else {
      get_walker_vcfs( somatic_vcfs, germline_vcfs, bams )
      walker_vcfs = get_walker_vcfs.out
    }

    if ( params.optional.short_variants.context_beds_dir ) {
      context_beds = extractContextBedFromDir( params.optional.short_variants.context_beds_dir )
    } else {
      get_context_beds( somatic_vcfs )
      context_beds = get_context_beds.out
    }

    if ( params.optional.short_variants.features_beds_dir ) {
      features_beds = extractFeaturesBedFromDir( params.optional.short_variants.features_beds_dir )
    } else {
      closest_features = Channel.from( params.features.closest )
      intersect_features = Channel.from( params.features.intersect )

      closest_feature( context_beds, closest_features )
      input_feature_beds = closest_feature.out.groupTuple( by: [0,1] )
      intersect_feature( context_beds, intersect_features )
      if ( input_feature_beds ) {
        input_feature_beds = intersect_feature.out.groupTuple( by: [0,1] ).join( input_feature_beds, by: [0,1] )
      } else {
        input_feature_beds = intersect_feature.out.groupTuple( by: [0,1] )
      }
      merge_features( context_beds, input_feature_beds )
      features_beds = merge_features.out
    }

    SplitVcfs( somatic_vcfs )

    somatic_snv_vcfs = SplitVcfs.out
      .map{
        donor_id, sample_id, snv_vcf, snv_tbi, indel_vcf, indel_tbi ->
        [donor_id, sample_id, snv_vcf, snv_tbi]
      }

    somatic_indel_vcfs = SplitVcfs.out
      .map{
        donor_id, sample_id, snv_vcf, snv_tbi, indel_vcf, indel_tbi ->
        [donor_id, sample_id, indel_vcf, indel_tbi]
      }

    if( params.run.snvs ) {
      snvs( ab_tables, features_beds, somatic_snv_vcfs, walker_vcfs )
      // snvs_ptato_vcfs = snvs.out
    } else {
      snvs_ptato_vcfs = Channel.empty()
    }

    // if ( params.run.indels ) {
    //   indels( ab_tables, features_beds, somatic_indel_vcfs, walker_vcfs )
    //   indels_ptato_vcfs = indels.out
    // } else {
    //   indels_ptato_vcfs = Channel.empty()
    // }
    //
    // intersect_ptato_vcfs( input_vcfs, snvs_ptato_vcfs, indels_ptato_vcfs )
    // ptato_intersect_vcfs = intersect_ptato_vcfs.out
    //
    // merge_ptato_vcfs( ptato_intersect_vcfs, snvs_ptato_vcfs, indels_ptato_vcfs )

}

workflow short_variants_train {
  take:
    input_vcfs
    germline_vcfs
  main:
    raw_pta_vcfs = extractPtaVcfFromDir( params.pta_vcfs_dir )
    raw_nopta_vcfs = extractNoptaVcfFromDir( params.nopta_vcfs_dir )
    raw_input_vcfs = raw_pta_vcfs.concat( raw_nopta_vcfs )

    get_gzipped_vcfs( raw_input_vcfs )
    somatic_vcfs = get_gzipped_vcfs.out

    pta_vcfs = raw_pta_vcfs
      .map{ [it[0], it[1]] }
      .combine( somatic_vcfs, by: [0,1] )

    nopta_vcfs = raw_nopta_vcfs
      .map{ [it[0], it[1]] }
      .combine( somatic_vcfs, by: [0,1] )

    if ( params.optional.short_variants.ab_tables_dir ) {
      ab_tables = extractAbTableFromDir( params.optional.short_variants.ab_tables_dir )
    } else {
      get_ab_tables( germline_vcfs, somatic_vcfs )
      ab_tables = get_ab_tables.out
    }

    if ( params.optional.short_variants.context_beds_dir ) {
      context_beds = extractContextBedFromDir( params.optional.short_variants.context_beds_dir )
    } else {
      get_context_beds( somatic_vcfs )
      context_beds = get_context_beds.out
    }

    if ( params.optional.short_variants.features_beds_dir ) {
      features_beds = extractFeaturesBedFromDir( params.optional.short_variants.features_beds_dir )
    } else {
      closest_features = Channel.from( params.features.closest )
      intersect_features = Channel.from( params.features.intersect )

      closest_feature( context_beds, closest_features )
      input_feature_beds = closest_feature.out.groupTuple( by: [0,1] )
      intersect_feature( context_beds, intersect_features )
      if ( input_feature_beds ) {
        input_feature_beds = intersect_feature.out.groupTuple( by: [0,1] ).join( input_feature_beds, by: [0,1] )
      } else {
        input_feature_beds = intersect_feature.out.groupTuple( by: [0,1] )
      }
      merge_features( context_beds, input_feature_beds )
      features_beds = merge_features.out
    }

    if( params.run.snvs || params.run.indels ) {
      label_info = pta_vcfs.map{
        donor_id, sample_id, vcf_files, tbi_files ->
        [ donor_id, sample_id ]
      }
      .combine( Channel.from('PTA') )
      .concat(
        nopta_vcfs.map{
          donor_id, sample_id, vcf_files, tbi_files ->
          [ donor_id, sample_id ]
        }
        .combine( Channel.from('noPTA') )
      )
      if ( params.run.snvs ) {
        snvs_train( ab_tables, features_beds, label_info )
      }
      if ( params.run.indels ) {
        indels_train( ab_tables, features_beds, label_info )
      }
    }
}
