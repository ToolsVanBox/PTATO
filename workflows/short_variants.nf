#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
  extractSomaticVcfGzFromDir;
  extractWalkerVcfGzFromDir;
  extractAbTableFromDir;
  extractContextBedFromDir;
  extractFeaturesBedFromDir;
  extractPtaVcfGzFromDir;
  extractNoptaVcfGzFromDir;
} from '../NextflowModules/Utils/getFilesFromDir.nf' params(params)

include { get_ab_tables } from './short_variants/get_ab_tables.nf' params(params)
include { get_walker_vcfs } from './short_variants/get_walker_vcfs.nf' params(params)
include { get_somatic_vcfs } from './short_variants/get_somatic_vcfs.nf' params(params)
include { get_context_beds } from './short_variants/get_context_beds.nf' params(params)
include { closest_feature } from './short_variants/get_features_beds.nf' params(params)
include { intersect_feature } from './short_variants/get_features_beds.nf' params(params)
include { merge_features } from './short_variants/get_features_beds.nf' params(params)

include { SplitVcfs } from '../NextflowModules/GATK/4.1.3.0/SplitVcfs.nf' params(params)
include { bgzip } from '../NextflowModules/htslib/1.15/bgzip.nf' params(params)
include { tabix } from '../NextflowModules/htslib/1.15/tabix.nf' params(params)

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
      somatic_vcfs = extractSomaticVcfGzFromDir( params.optional.short_variants.somatic_vcfs_dir )
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
      walker_vcfs = extractWalkerVcfGzFromDir( params.optional.short_variants.walker_vcfs_dir )
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
    }

    if ( params.run.indels ) {
      indels( ab_tables, features_beds, somatic_indel_vcfs, walker_vcfs )
    }
}

workflow short_variants_train {
  take:
    input_vcfs
    bams
    germline_vcfs
  main:
    pta_vcfs = extractPtaVcfGzFromDir( params.pta_vcfs_dir )
    nopta_vcfs = extractNoptaVcfGzFromDir( params.nopta_vcfs_dir )

    somatic_vcfs = pta_vcfs.concat( nopta_vcfs )

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
