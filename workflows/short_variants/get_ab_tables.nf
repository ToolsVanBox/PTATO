include { extractPhasedVcfGzFromDir } from '../../NextflowModules/Utils/getFilesFromDir.nf' params(params)
include { shapeit } from '../../NextflowModules/shapeit/4.2.2/shapeit.nf' params(params)
include {
  createABtable;
  mergeABtable
} from '../../NextflowModules/Utils/ABtable.nf' params(params)
include { tabix } from '../../NextflowModules/htslib/1.15/tabix.nf' params(params)


workflow get_ab_tables {
  take:
    input_germline_vcfs
    input_somatic_vcfs

  main:
    chroms_file = file(params.shapeit.chroms)
    chroms = chroms_file.readLines()
    input_germline_vcfs_chroms = input_germline_vcfs.combine( chroms )
    bulk_names = Channel.from( params.bulk_names )

    if ( params.optional.short_variants.phased_vcfs_dir ) {
      phased_vcfs = extractPhasedVcfGzFromDir( params.optional.short_variants.phased_vcfs_dir )
    } else {
      shapeit( input_germline_vcfs_chroms )
      tabix( shapeit.out
        .map{ donor_id, sample_id, chrom, phased_vcf ->
          [ donor_id, sample_id, phased_vcf ]
        }
      )

      phased_vcfs = tabix.out
        .map{ donor_id, sample_id, vcf_gz, vcf_tbi ->
          vcf_name = vcf_gz.getName()
          tbi_name = vcf_tbi.getName()
          vcf_gz = vcf_gz.copyTo("${params.out_dir}/intermediate/short_variants/shapeit/${donor_id}/${vcf_name}")
          vcf_tbi = vcf_tbi.copyTo("${params.out_dir}/intermediate/short_variants/shapeit/${donor_id}/${tbi_name}")
          m = vcf_name =~ /.+\.(.+?)\.phased.vcf.gz/
          chrom = m[0][1]
          [ donor_id, sample_id, chrom, vcf_gz, vcf_tbi ]
        }
    }

    createABtable(
      phased_vcfs
        .combine( input_germline_vcfs, by: 0 )
        .combine( input_somatic_vcfs, by: 0 )
        .combine( bulk_names, by: 0 )
        .map{
          donor_id, phased_sample_id, chrom, phased_vcf_gz, phased_tbi, germline_sample_id, germline_vcf_gz, germline_tbi, somatic_sample_id, somatic_vcf_gz, somatic_tbi, bulk_name ->
          [ donor_id, somatic_sample_id, chrom, phased_vcf_gz, phased_tbi, germline_vcf_gz, germline_tbi, somatic_vcf_gz, somatic_tbi, bulk_name ]
        }
    )

    mergeABtable (
      createABtable.out
        .groupTuple( by: [0,1] )
    )

    ab_tables = mergeABtable.out
      .map{
        donor_id, sample_id, ab_table ->
        file_name = ab_table.getName()
        ab_table = ab_table.copyTo("${params.out_dir}/intermediate/short_variants/ab/${donor_id}/${file_name}")
        [ donor_id, sample_id, ab_table ]
      }

  emit:
    ab_tables
}
