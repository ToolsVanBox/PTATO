include { createIndelExcludeIndelList } from '../../NextflowModules/Utils/createExcludeIndelList.nf' params(params)
include { sort } from '../../NextflowModules/bedtools/2.30.0/sort.nf'

workflow create_excludeindellist_bed {
  take:
    features_beds
    label_info
  main:
    input_excludeindellist = features_beds
      .combine( label_info, by: [0,1] )
      .map{
        donor_id, sample_id, feature_bed, label ->
        [ label, feature_bed ]
      }
      .groupTuple( by: [0] )
      .filter{ it[0] == 'PTA' }

    createIndelExcludeIndelList( input_excludeindellist )
    sort( createIndelExcludeIndelList.out )
    excludeindellist_bed = sort.out
      .map{ donor_id, sample_id, bed ->
        bed_name = bed.getName()
        bed = bed.copyTo("${params.out_dir}/indels/excludeindellist/${bed_name}")
        [ donor_id, sample_id, bed ]
      }

  emit:
      excludeindellist_bed
}
