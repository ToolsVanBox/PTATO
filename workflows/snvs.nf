#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { get_snvs_rf_tables } from './snvs/get_rf_tables.nf' params(params)
include { get_ptato_vcfs } from './snvs/get_ptato_vcfs.nf' params(params)
include { filter_ptato_vcfs } from './snvs/filter_ptato_vcfs.nf' params(params)

include { get_snv_rf_files } from './snvs/get_rf_files.nf' params(params)


workflow snvs {
  take:
    ab_tables
    features_beds
    somatic_snv_vcfs
    walker_vcfs
  main :
    get_snvs_rf_tables( ab_tables, features_beds )
    rf_tables = get_snvs_rf_tables.out

    get_ptato_vcfs( somatic_snv_vcfs, rf_tables )
    ptato_vcfs = get_ptato_vcfs.out

    filter_ptato_vcfs( ptato_vcfs, walker_vcfs )
    ptato_filtered_vcfs = filter_ptato_vcfs.out

    ptato_combined_vcfs = ptato_vcfs.combine( ptato_filtered_vcfs, by: [0,1])

  emit:
    ptato_combined_vcfs
}

workflow snvs_train{
  take:
    ab_tables
    features_beds
    label_info
  main:
    get_snvs_rf_tables( ab_tables, features_beds )
    rf_tables = get_snvs_rf_tables.out

    get_snv_rf_files( rf_tables, label_info )
    rf_files = get_snv_rf_files.out
}
