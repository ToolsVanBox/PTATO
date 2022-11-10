#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { get_indel_rf_tables } from './indels/get_rf_tables.nf' params(params)
include { get_excludeindellist_filtered_vcfs } from './indels/get_excludeindellist_filtered_vcfs.nf' params(params)
include { get_ptato_vcfs } from './indels/get_ptato_vcfs.nf' params(params)
include { filter_ptato_vcfs } from './indels/filter_ptato_vcfs.nf' params(params)


include { get_indel_rf_files } from './indels/get_rf_files.nf' params(params)
include { create_excludeindellist_bed } from './indels/create_excludeindellist_bed.nf' params(params)

workflow indels {
  take:
    somatic_indel_vcfs
    context_beds

  main :
    get_ptato_vcfs( somatic_indel_vcfs, context_beds )
    ptato_vcfs = get_ptato_vcfs.out

    filter_ptato_vcfs( ptato_vcfs )
    ptato_filtered_vcfs = filter_ptato_vcfs.out

  emit:
    ptato_vcfs
}

workflow indels_train {
  take:
    ab_tables
    features_beds
    label_info
  main:
    create_excludeindellist_bed( features_beds, label_info )

    get_indel_rf_tables( ab_tables, features_beds )
    rf_tables = get_indel_rf_tables.out

    get_indel_rf_files( rf_tables, label_info )
    rf_files = get_indel_rf_files.out
}
