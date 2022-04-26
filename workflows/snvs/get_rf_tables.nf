include { extractSnvRfTableRdsFromDir } from '../../NextflowModules/Utils/getFilesFromDir' params(params)
include { createSnvRfTable } from '../../NextflowModules/Utils/createRfTable.nf' params(params)

workflow get_snvs_rf_tables {
  take:
    ab_tables
    features_beds
  main:
    if ( params.optional.snvs.rf_tables_dir ) {
      rf_tables = extractSnvRfTableRdsFromDir( params.optional.snvs.rf_tables_dir )
    } else {
      createSnvRfTable( ab_tables.combine( features_beds, by: [0,1] ) )
      rf_tables = createSnvRfTable.out
        .map{
          donor_id, sample_id, rf_table_rds ->
          rds_name = rf_table_rds.getName()
          rf_table_rds = rf_table_rds.copyTo("${params.out_dir}/snvs/intermediate/rf/${donor_id}/${rds_name}")
          [ donor_id, sample_id, rf_table_rds ]
        }
    }
  emit:
    rf_tables
}
