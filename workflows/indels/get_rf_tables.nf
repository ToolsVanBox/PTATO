include { extractIndelRfTableRdsFromDir } from '../../NextflowModules/Utils/getFilesFromDir' params(params)
include { createIndelRfTable } from '../../NextflowModules/Utils/createRfTable.nf' params(params)

workflow get_indel_rf_tables {
  take:
    ab_tables
    features_beds
  main:
    if ( params.optional.indels.rf_tables_dir ) {
      rf_tables = extractIndelRfTableRdsFromDir( params.optional.indels.rf_tables_dir )
    } else {
      createIndelRfTable( ab_tables.combine( features_beds, by: [0,1] ) )
      rf_tables = createIndelRfTable.out
        .map{
          donor_id, sample_id, rf_table_rds ->
          rds_name = rf_table_rds.getName()
          rf_table_rds = rf_table_rds.copyTo("${params.out_dir}/indels/intermediate/rf/${donor_id}/${rds_name}")
          [ donor_id, sample_id, rf_table_rds ]
        }
    }
  emit:
    rf_tables
}
