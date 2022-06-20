include { getContext } from '../../NextflowModules/Utils/getContext.nf' params(params)
include { sort } from '../../NextflowModules/bedtools/2.30.0/sort.nf' params(params)

workflow get_context_beds {
  take:
    vcfs
  main:
    getContext( vcfs )
    sort( getContext.out )
    context_beds = sort.out
      .map{ donor_id, sample_id, bed ->
        bed_name = bed.getName().toString().replace(".sorted.bed",".context.sorted.bed")
        bed = bed.copyTo("${params.out_dir}/intermediate/short_variants/context/${donor_id}/${bed_name}")
        [ donor_id, sample_id, bed ]
      }

  emit:
    context_beds
}
