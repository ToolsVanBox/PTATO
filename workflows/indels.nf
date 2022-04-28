#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { get_indel_rf_tables } from './indels/get_rf_tables.nf' params(params)
include { get_ptato_indel_vcfs } from './indels/get_ptato_vcfs.nf' params(params)
include { filter_ptato_vcfs } from './indels/filter_ptato_vcfs.nf' params(params)

include { get_indel_rf_files } from './indels/get_rf_files.nf' params(params)
include { create_blacklist_bed } from './indels/create_blacklist_bed.nf' params(params)

workflow indels {
  take:
    ab_tables
    features_beds
    somatic_indel_vcfs
    walker_vcfs
  main :
    get_indel_rf_tables( ab_tables, features_beds )
    rf_tables = get_indel_rf_tables.out

    get_ptato_indel_vcfs( somatic_indel_vcfs, rf_tables )
    ptato_vcfs = get_ptato_indel_vcfs.out

    filter_ptato_vcfs( ptato_vcfs, walker_vcfs )
    ptato_filtered_vcfs = filter_ptato_vcfs.out
}

workflow indels_train {
  take:
    ab_tables
    features_beds
    label_info
  main:
    create_blacklist_bed( features_beds, label_info )

    get_indel_rf_tables( ab_tables, features_beds )
    rf_tables = get_indel_rf_tables.out

    get_indel_rf_files( rf_tables, label_info )
    rf_files = get_indel_rf_files.out
}
