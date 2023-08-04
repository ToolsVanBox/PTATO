def extractNoptaVcfFromDir( nopta_dir ) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    nopta_dir = nopta_dir.tokenize().collect{"$it/*/*.{vcf,vcf.gz}"}
    Channel
      .fromPath(nopta_dir, type:'file')
      .ifEmpty { error "No .vcf(.gz) files found in ${nopta_dir}." }
      .map { nopta_vcf_path ->
          nopta_vcf_file = nopta_vcf_path
          nopta_tbi_file = nopta_vcf_path+".tbi"
          nopta_sample_id = nopta_vcf_path.getName().toString().replaceAll(/.vcf(.gz)*$/, '')
          file(nopta_vcf_path.getBaseName()).getBaseName()
          nopta_donor_id = nopta_vcf_path.getParent().getName()
          if ( file(nopta_tbi_file).exists() ) {
            [nopta_donor_id, nopta_sample_id, nopta_vcf_file, nopta_tbi_file]
          } else {
            [nopta_donor_id, nopta_sample_id, nopta_vcf_file]
          }
      }
}

def extractPtaVcfFromDir( pta_dir ) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    pta_dir = pta_dir.tokenize().collect{"$it/*/*.{vcf,vcf.gz}"}
    Channel
      .fromPath(pta_dir, type:'file')
      .ifEmpty { error "No .vcf(.gz) files found in ${pta_dir}." }
      .map { pta_vcf_path ->
          pta_vcf_file = pta_vcf_path
          pta_tbi_file = pta_vcf_path+".tbi"
          pta_sample_id = pta_vcf_path.getName().toString().replaceAll(/.vcf(.gz)*$/, '')
          pta_donor_id =  pta_vcf_path.getParent().getName()
          if ( file(pta_tbi_file).exists() ) {
            [pta_donor_id, pta_sample_id, pta_vcf_file, pta_tbi_file]
          } else {
            [pta_donor_id, pta_sample_id, pta_vcf_file]
          }

      }
}

def extractInputVcfFromDir( input_dir ) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    input_dir = input_dir.tokenize().collect{"$it/*/*.{vcf,vcf.gz}"}
    Channel
      .fromPath(input_dir, type:'file')
      .ifEmpty { error "No .vcf(.gz) files found in ${input_dir}." }
      .map { input_vcf_path ->
          input_vcf_file = input_vcf_path
          input_tbi_file = input_vcf_path+".tbi"
          input_sample_id = input_vcf_path.getName().toString().replaceAll(/.vcf(.gz)*$/, '')
          input_donor_id =  input_vcf_path.getParent().getName()
          if ( file(input_tbi_file).exists() ) {
            [input_donor_id, input_sample_id, input_vcf_file, input_tbi_file]
          } else {
            [input_donor_id, input_sample_id, input_vcf_file]
          }
      }
}

def extractGermlineVcfFromDir( germline_vcfs_dir ) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    germline_vcfs_dir = germline_vcfs_dir.tokenize().collect{"$it/*/*.{vcf,vcf.gz}"}
    Channel
      .fromPath(germline_vcfs_dir, type:'file')
      .ifEmpty { error "No .vcf(.gz) files found in ${germline_vcfs_dir}." }
      .map { germline_vcf_path ->
          germline_vcf_file = germline_vcf_path
          germline_tbi_file = germline_vcf_path+".tbi"
          germline_sample_id = germline_vcf_path.getName().toString().replaceAll(/(.germline)*.vcf(.gz)*$/, '')
          germline_donor_id =  germline_vcf_path.getParent().getName()
          if ( file(germline_tbi_file).exists() ) {
            [germline_donor_id, germline_sample_id, germline_vcf_file, germline_tbi_file]
          } else {
            [germline_donor_id, germline_sample_id, germline_vcf_file]
          }
      }
}

def extractSomaticVcfFromDir( somatic_dir ) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    somatic_dir = somatic_dir.tokenize().collect{"$it/*/*.{vcf,vcf.gz}"}
    Channel
      .fromPath(somatic_dir, type:'file')
      .ifEmpty { error "No .vcf(.gz) files found in ${somatic_dir}." }
      .map { somatic_vcf_path ->
          somatic_vcf_file = somatic_vcf_path
          somatic_tbi_file = somatic_vcf_path+".tbi"
          somatic_sample_id = somatic_vcf_path.getName().toString().replaceAll(/(.SMuRF.filtered)*.vcf(.gz)*$/, '')
          somatic_donor_id =  somatic_vcf_path.getParent().getName()
          if ( file(somatic_tbi_file).exists() ) {
            [somatic_donor_id, somatic_sample_id, somatic_vcf_file, somatic_tbi_file]
          } else {
            [somatic_donor_id, somatic_sample_id, somatic_vcf_file]
          }
      }
}

def extractBamsFromDir( bams_dir ) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    bams_dir = bams_dir.tokenize().collect{"$it/*/*.bam"}
    Channel
      .fromPath(bams_dir, type:'file')
      .ifEmpty { error "No bam files found in ${bams_dir}." }
      .map { bam_path ->
          bam_file = bam_path
          bai_file = bam_path.toString().replace('.bam','.bai')
          bam_sample_id = bam_path.getName().toString().replaceAll(/(_dedup)*.bam$/, '')
          bam_donor_id =  bam_path.getParent().getName()
          if ( file(bai_file).exists() ) {
            [bam_donor_id, bam_sample_id, bam_file, bai_file]
          } else {
            bai_file = bam_path+".bai"
            if ( file(bai_file).exists() ) {
              [bam_donor_id, bam_sample_id, bam_file, bai_file]
            } else {
              [bam_donor_id, bam_sample_id, bam_file]
            }
          }
      }
}

def extractWalkerVcfFromDir( walker_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  walker_dir = walker_dir.tokenize().collect{"$it/*/*.{vcf,vcf.gz}"}
  Channel
    .fromPath(walker_dir, type:'file')
    .ifEmpty { error "No .vcf(.gz) files found in ${walker_dir}." }
    .map { walker_vcf_path ->
        walker_vcf_file = walker_vcf_path
        walker_tbi_file = walker_vcf_path+".tbi"
        walker_sample_id = walker_vcf_path.getName().toString().replaceAll(/(.walker)*.vcf(.gz)*$/, '')
        walker_donor_id =  walker_vcf_path.getParent().getName()
        if ( file(walker_tbi_file).exists() ) {
          [walker_donor_id, walker_sample_id, walker_vcf_file, walker_tbi_file]
        } else {
          [walker_donor_id, walker_sample_id, walker_vcf_file]
        }
    }
}

def extractAbTableFromDir( ab_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  ab_dir = ab_dir.tokenize().collect{"$it/*/*.abtable.txt"}
  Channel
    .fromPath(ab_dir, type:'file')
    .ifEmpty { error "No .abtable.txt files found in ${ab_dir}." }
    .map { ab_path ->
        ab_file = ab_path
        ab_sample_id = ab_path.getName().toString().replaceAll(/.abtable.txt$/, '')
        ab_donor_id =  ab_path.getParent().getName()
        [ab_donor_id, ab_sample_id, ab_file]
    }
}

def extractContextBedFromDir( context_bed_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  context_bed_dir = context_bed_dir.tokenize().collect{"$it/*/*.bed"}
  Channel
    .fromPath(context_bed_dir, type:'file')
    .ifEmpty { error "No .bed files found in ${context_bed_dir}." }
    .map { context_bed_path ->
        context_bed_file = context_bed_path
        context_sample_id = context_bed_path.getName().toString().replaceAll(/(.context)*(.sorted)*.bed$/, '')
        context_donor_id =  context_bed_path.getParent().getName()
        [context_donor_id, context_sample_id, context_bed_file]
    }
}

def extractFeaturesBedFromDir( features_bed_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  features_bed_dir = features_bed_dir.tokenize().collect{"$it/*/*.bed"}
  Channel
    .fromPath(features_bed_dir, type:'file')
    .ifEmpty { error "No .bed files found in ${features_bed_dir}." }
    .map { features_bed_path ->
        features_bed_file = features_bed_path
        features_sample_id = features_bed_path.getName().toString().replaceAll(/(.features)*(.groupby)*.bed$/, '')
        features_donor_id =  features_bed_path.getParent().getName()
        [features_donor_id, features_sample_id, features_bed_file]
    }
}

def extractPhasedVcfGzFromDir( phased_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  phased_dir = phased_dir.tokenize().collect{"$it/*/*.phased.vcf.gz"}
  Channel
    .fromPath(phased_dir, type:'file')
    .ifEmpty { error "No .phased.vcf.gz files found in ${phased_dir}." }
    .map { phased_vcf_path ->
        phased_vcf_file = phased_vcf_path
        phased_tbi_file = phased_vcf_path+".tbi"
        phased_sample_id = phased_vcf_path.getName().toString().replaceAll(/\.(.+?)\.phased.vcf.gz$/, '')
        phased_donor_id =  phased_vcf_path.getParent().getName()
        m = phased_vcf_path =~ /.+\.(.+?)\.phased.vcf.gz/
        phased_chrom = m[0][1]
        [phased_donor_id, phased_sample_id, phased_chrom, phased_vcf_file, phased_tbi_file]
    }
}

def extractExcludeIndelListFilteredVcfFromDir( excludeindellist_filtered_vcfs_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  excludeindellist_filtered_vcfs_dir = excludeindellist_filtered_vcfs_dir.tokenize().collect{"$it/*/*.{vcf,vcf.gz}"}
  Channel
    .fromPath(excludeindellist_filtered_vcfs_dir, type:'file')
    .ifEmpty { error "No .vcf(.gz) files found in ${excludeindellist_filtered_vcfs_dir}." }
    .map { excludeindellist_filtered_vcf_path ->
        excludeindellist_filtered_vcf_file = excludeindellist_filtered_vcf_path
        excludeindellist_filtered_tbi_file = excludeindellist_filtered_vcf_path+".tbi"
        excludeindellist_filtered_sample_id = excludeindellist_filtered_vcf_path.getName().toString().replaceAll(/(.excludeindellist.filtered)*.vcf(.gz)*$/, '')
        excludeindellist_filtered_donor_id =  excludeindellist_filtered_vcf_path.getParent().getName()
        if ( file(excludeindellist_filtered_tbi_file).exists() ) {
          [excludeindellist_filtered_donor_id, excludeindellist_filtered_sample_id, excludeindellist_filtered_vcf_file, excludeindellist_filtered_tbi_file]
        } else {
          [excludeindellist_filtered_donor_id, excludeindellist_filtered_sample_id, excludeindellist_filtered_vcf_file]
        }
    }
}

def extractSnvRfTableRdsFromDir( snv_rf_tables_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  snv_rf_tables_dir = snv_rf_tables_dir.tokenize().collect{"$it/*/*.rftable.rds"}
  Channel
    .fromPath(snv_rf_tables_dir, type:'file')
    .ifEmpty { error "No .rftable.rds files found in ${snv_rf_tables_dir}." }
    .map { snv_rf_rds_path ->
        snv_rf_rds_file = snv_rf_rds_path
        snv_rf_rds_sample_id = snv_rf_rds_path.getName().toString().replaceAll(/.rftable.rds$/, '')
        snv_rf_rds_donor_id =  snv_rf_rds_path.getParent().getName()
        [snv_rf_rds_donor_id, snv_rf_rds_sample_id, snv_rf_rds_file]
    }
}

def extractIndelRfTableRdsFromDir( indel_rf_tables_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  indel_rf_tables_dir = indel_rf_tables_dir.tokenize().collect{"$it/*/*.rftable.rds"}
  Channel
    .fromPath(indel_rf_tables_dir, type:'file')
    .ifEmpty { error "No .rftable.rds files found in ${indel_rf_tables_dir}." }
    .map { indel_rf_rds_path ->
        indel_rf_rds_file = indel_rf_rds_path
        indel_rf_rds_sample_id = indel_rf_rds_path.getName().toString().replaceAll(/.rftable.rds$/, '')
        indel_rf_rds_donor_id =  indel_rf_rds_path.getParent().getName()
        [indel_rf_rds_donor_id, indel_rf_rds_sample_id, indel_rf_rds_file]
    }
}

def extractPtatoVcfFromDir( ptato_vcfs_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  ptato_vcfs_dir = ptato_vcfs_dir.tokenize().collect{"$it/*/*.{vcf,vcf.gz}"}
  Channel
    .fromPath(ptato_vcfs_dir, type:'file')
    .ifEmpty { error "No .vcf(.gz) files found in ${ptato_vcfs_dir}." }
    .map { ptato_vcf_path ->
        ptato_vcf_file = ptato_vcf_path
        ptato_tbi_file = ptato_vcf_path+".tbi"
        ptato_sample_id = ptato_vcf_path.getName().toString().replaceAll(/(.ptato)*.vcf(.gz)*$/, '')
        ptato_donor_id =  ptato_vcf_path.getParent().getName()
        if ( file(ptato_tbi_file).exists() ) {
          [ptato_donor_id, ptato_sample_id, ptato_vcf_file, ptato_tbi_file]
        } else {
          [ptato_donor_id, ptato_sample_id, ptato_vcf_file]
        }
    }
}

// You are adding here 
def extractPtatoTableFromDir( ptato_vcfs_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  ptato_vcfs_dir = ptato_vcfs_dir.tokenize().collect{"$it/*/*/*.txt"}
  Channel
    .fromPath(ptato_vcfs_dir, type:'file')
    .ifEmpty { error "No .txt files found in ${ptato_vcfs_dir}." }
    .map { ptato_vcf_path ->
        ptato_table_file = ptato_vcf_path
        ptato_sample_id = ptato_vcf_path.getName().toString().replaceAll(/(.ptatotable)*.txt*$/, '')
        ptato_donor_id =  ptato_vcf_path.getParent().getParent().getName()
        [ptato_donor_id, ptato_sample_id, ptato_table_file]
    }
}
// You are adding here 

def extractCombinedPtatoVcfFromDir( ptato_vcfs_dir ){
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  ptato_vcfs_dir = ptato_vcfs_dir.tokenize().collect{"$it/*/*.{vcf,vcf.gz}"}
  Channel
    .fromPath(ptato_vcfs_dir, type:'file')
    .ifEmpty { error "No .vcf(.gz) files found in ${ptato_vcfs_dir}." }
    .map { ptato_vcf_path ->
        ptato_vcf_file = ptato_vcf_path
        ptato_tbi_file = ptato_vcf_path+".tbi"
        sample_id = ptato_vcf_path.getName().toString().replaceAll(/.snvs.ptato.vcf(.gz)*$/, '')
        donor_id =  ptato_vcf_path.getParent().getName()
        ptato_filt_vcf_file = ptato_vcf_file.getParent() / sample_id / sample_id+".snvs.ptato.filtered.vcf.gz"
        ptato_filt_tbi_file = ptato_filt_vcf_file+".tbi"
        [donor_id, sample_id, ptato_vcf_file, ptato_tbi_file, ptato_filt_vcf_file, ptato_filt_tbi_file]
    }
}

def extractWGSMetricsFromDir( wgs_metrics_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  wgs_metrics_dir = wgs_metrics_dir.tokenize().collect{"$it/*/*.wgs_metrics*"}
  Channel
    .fromPath(wgs_metrics_dir, type:'file')
    .ifEmpty { error "No .wgs_metrics* files found in ${wgs_metrics_dir}." }
    .map { wgs_metrics_path ->
        wgs_metrics_file = wgs_metrics_path
        wgs_metrics_sample_id = wgs_metrics_path.getName().toString().replaceAll(/.wgs_metrics.*$/, '')
        wgs_metrics_donor_id =  wgs_metrics_path.getParent().getName()
        [wgs_metrics_donor_id, wgs_metrics_sample_id, wgs_metrics_file]
    }
}

def extractAlignmentSummaryMetricsFromDir( alignment_summary_metrics_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  alignment_summary_metrics_dir = alignment_summary_metrics_dir.tokenize().collect{"$it/*/*.alignment_summary_metrics*"}
  Channel
    .fromPath(alignment_summary_metrics_dir, type:'file')
    .ifEmpty { error "No .alignment_summary_metrics* files found in ${alignment_summary_metrics_dir}." }
    .map { alignment_summary_metrics_path ->
        alignment_summary_metrics_file = alignment_summary_metrics_path
        alignment_summary_metrics_sample_id = alignment_summary_metrics_path.getName().toString().replaceAll(/(.multiple_metrics)*.alignment_summary_metrics.*$/, '')
        alignment_summary_metrics_donor_id =  alignment_summary_metrics_path.getParent().getName()
        [alignment_summary_metrics_donor_id, alignment_summary_metrics_sample_id, alignment_summary_metrics_file]
    }
}

def extractAutosomalCallableLociFromDir( callableloci_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  callableloci_dir = callableloci_dir.tokenize().collect{"$it/*/*.callableloci.autosomal.txt"}
  Channel
    .fromPath(callableloci_dir, type:'file')
    .ifEmpty { error "No .callableloci.bed files found in ${callableloci_dir}." }
    .map { callableloci_path ->
        autosomalcallableloci_file = callableloci_path
        autosomalcallableloci_sample_id = callableloci_path.getName().toString().replaceAll(/.callableloci.autosomal.txt$/, '')
        autosomalcallableloci_donor_id =  callableloci_path.getParent().getName()
        [autosomalcallableloci_donor_id, autosomalcallableloci_sample_id, autosomalcallableloci_file]
    }
}

def extractCallableLociBedFromDir( callableloci_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  callableloci_dir = callableloci_dir.tokenize().collect{"$it/*/*.callableloci.bed"}
  Channel
    .fromPath(callableloci_dir, type:'file')
    .ifEmpty { error "No .callableloci.bed files found in ${callableloci_dir}." }
    .map { callableloci_path ->
        callableloci_bed = callableloci_path
        callableloci_txt = callableloci_bed.toString().replace('.bed','.txt')
        callableloci_sample_id = callableloci_path.getName().toString().replaceAll(/.callableloci.bed$/, '')
        callableloci_donor_id =  callableloci_path.getParent().getName()
        if ( file(callableloci_txt).exists() ) {
          [callableloci_donor_id, callableloci_sample_id, callableloci_bed, callableloci_txt]
        }
    }
}

def extractGridssDriverVcfFromDir( gridss_driver_vcfs_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  gridss_driver_vcfs_dir = gridss_driver_vcfs_dir.tokenize().collect{"$it/*/*.gridss.driver.{vcf,vcf.gz}"}
  Channel
    .fromPath(gridss_driver_vcfs_dir, type:'file')
    .ifEmpty { error "No .gridss.driver.vcf(.gz) files found in ${gridss_driver_vcfs_dir}." }
    .map { gridss_driver_vcf_path ->
        gridss_driver_vcf_file = gridss_driver_vcf_path
        gridss_driver_tbi_file = gridss_driver_vcf_path+".tbi"
        gridss_driver_tumor_sample_id = gridss_driver_vcf_path.getName().toString().replaceAll(/(.gridss.driver)*.vcf(.gz)*$/, '')
        gridss_driver_normal_sample_id =  gridss_driver_vcf_path.getParent().getName()
        gridss_driver_donor_id =  gridss_driver_vcf_path.getParent().getParent().getName()
        if ( file(gridss_driver_tbi_file).exists() ) {
          [gridss_driver_donor_id, gridss_driver_normal_sample_id+";"+gridss_driver_tumor_sample_id, gridss_driver_vcf_file, gridss_driver_tbi_file]
        } else {
          [gridss_driver_donor_id, gridss_driver_normal_sample_id+";"+gridss_driver_tumor_sample_id, gridss_driver_vcf_file]
        }
    }
}


def extractGridssUnfilteredVcfFromDir( gridss_unfiltered_vcfs_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  gridss_unfiltered_vcfs_dir = gridss_unfiltered_vcfs_dir.tokenize().collect{"$it/*/*/*.gridss.unfiltered.{vcf,vcf.gz}"}
  Channel
    .fromPath(gridss_unfiltered_vcfs_dir, type:'file')
    .ifEmpty { error "No .gridss.unfiltered.vcf(.gz) files found in ${gridss_unfiltered_vcfs_dir}." }
    .map { gridss_unfiltered_vcf_path ->
        gridss_unfiltered_vcf_file = gridss_unfiltered_vcf_path
        gridss_unfiltered_tbi_file = gridss_unfiltered_vcf_path+".tbi"
        gridss_unfiltered_tumor_sample_id = gridss_unfiltered_vcf_path.getName().toString().replaceAll(/(.gridss.unfiltered)*.vcf(.gz)*$/, '')
        gridss_unfiltered_normal_sample_id =  gridss_unfiltered_vcf_path.getParent().getName()
        gridss_unfiltered_donor_id =  gridss_unfiltered_vcf_path.getParent().getParent().getName()
        if ( file(gridss_unfiltered_tbi_file).exists() ) {
          [gridss_unfiltered_donor_id, gridss_unfiltered_normal_sample_id+";"+gridss_unfiltered_tumor_sample_id, gridss_unfiltered_vcf_file, gridss_unfiltered_tbi_file]
        } else {
          [gridss_unfiltered_donor_id, gridss_unfiltered_normal_sample_id+";"+gridss_unfiltered_tumor_sample_id, gridss_unfiltered_vcf_file]
        }
    }
}

def extractGripssSomaticFilteredVcfFromDir( gripss_somatic_filtered_vcfs_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  gripss_somatic_filtered_vcfs_dir = gripss_somatic_filtered_vcfs_dir.tokenize().collect{"$it/*/*/*.gripss.somatic.filtered.{vcf,vcf.gz}"}
  Channel
    .fromPath(gripss_somatic_filtered_vcfs_dir, type:'file')
    .ifEmpty { error "No .gripss.somatic.filtered.vcf(.gz) files found in ${gripss_somatic_filtered_vcfs_dir}." }
    .map { gripss_somatic_filtered_vcf_path ->
        gripss_somatic_filtered_vcf_file = gripss_somatic_filtered_vcf_path
        gripss_somatic_filtered_tbi_file = gripss_somatic_filtered_vcf_path+".tbi"
        gripss_somatic_filtered_tumor_sample_id = gripss_somatic_filtered_vcf_path.getName().toString().replaceAll(/(.gripss.somatic.filtered)*.vcf(.gz)*$/, '')
        gripss_somatic_filtered_normal_sample_id =  gripss_somatic_filtered_vcf_path.getParent().getName()
        gripss_somatic_filtered_donor_id =  gripss_somatic_filtered_vcf_path.getParent().getParent().getName()
        if ( file(gripss_somatic_filtered_tbi_file).exists() ) {
          [gripss_somatic_filtered_donor_id, gripss_somatic_filtered_normal_sample_id+";"+gripss_somatic_filtered_tumor_sample_id, gripss_somatic_filtered_vcf_file, gripss_somatic_filtered_tbi_file]
        } else {
          [gripss_somatic_filtered_donor_id, gripss_somatic_filtered_normal_sample_id+";"+gripss_somatic_filtered_tumor_sample_id, gripss_somatic_filtered_vcf_file]
        }
    }
}

def extractCobaltRatioTsvFromDir( cobalt_ratio_tsv_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  cobalt_ratio_tsv_dir = cobalt_ratio_tsv_dir.tokenize().collect{"$it/*/*/*.cobalt.ratio.tsv"}
  Channel
    .fromPath(cobalt_ratio_tsv_dir, type:'file')
    .ifEmpty { error "No .cobalt.ratio.tsv files found in ${cobalt_ratio_tsv_dir}." }
    .map { cobalt_ratio_tsv_path ->
        cobalt_ratio_tsv_file = cobalt_ratio_tsv_path
        cobalt_ratio_tumor_sample_id = cobalt_ratio_tsv_path.getName().toString().replaceAll(/.cobalt.ratio.tsv$/, '')
        cobalt_ratio_normal_sample_id =  cobalt_ratio_tsv_path.getParent().getName()
        cobalt_ratio_donor_id =  cobalt_ratio_tsv_path.getParent().getParent().getName()
        [cobalt_ratio_donor_id, cobalt_ratio_normal_sample_id, cobalt_ratio_tumor_sample_id, cobalt_ratio_tsv_file]
    }
}

def extractCobaltFilteredReadCounts( cobalt_filtered_readcounts_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  cobalt_filtered_readcounts_dir = cobalt_filtered_readcounts_dir.tokenize().collect{"$it/*/*/*.readcounts.*.txt"}
  Channel
    .fromPath(cobalt_filtered_readcounts_dir, type:'file')
    .ifEmpty { error "No .readcounts.*.txt files found in ${cobalt_filtered_readcounts_dir}." }
    .map { cobalt_filtered_readcounts_path ->
        cobalt_filtered_readcounts_file = cobalt_filtered_readcounts_path
        cobalt_filtered_readcounts_tumor_sample_id = cobalt_filtered_readcounts_path.getName().toString().replaceAll(/.readcounts.filtered(.*).txt$/, '')
        cobalt_filtered_readcounts_normal_sample_id =  cobalt_filtered_readcounts_path.getParent().getName()
        cobalt_filtered_readcounts_donor_id =  cobalt_filtered_readcounts_path.getParent().getParent().getName()
        [cobalt_filtered_readcounts_donor_id, cobalt_filtered_readcounts_normal_sample_id, cobalt_filtered_readcounts_tumor_sample_id, cobalt_filtered_readcounts_file]
    }
}

def extractCobaltFilteredReadCountsSegments( cobalt_filtered_readcounts_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  cobalt_filtered_readcounts_dir = cobalt_filtered_readcounts_dir.tokenize().collect{"$it/*/*/*.readcounts.segments.*"}
  Channel
    .fromPath(cobalt_filtered_readcounts_dir, type:'file')
    .ifEmpty { error "No .readcounts.segments.* files found in ${cobalt_filtered_readcounts_dir}." }
    .map { cobalt_filtered_readcounts_bedpe_path ->
        cobalt_filtered_readcounts_bedpe = cobalt_filtered_readcounts_bedpe_path
        cobalt_filtered_readcounts_bedpe_tumor_sample_id = cobalt_filtered_readcounts_bedpe_path.getName().toString().replaceAll(/.readcounts.segments.*$/, '')
        cobalt_filtered_readcounts_bedpe_normal_sample_id =  cobalt_filtered_readcounts_bedpe_path.getParent().getName()
        cobalt_filtered_readcounts_bedpe_donor_id = cobalt_filtered_readcounts_bedpe_path.getParent().getParent().getName()
        [cobalt_filtered_readcounts_bedpe_donor_id, cobalt_filtered_readcounts_bedpe_normal_sample_id, cobalt_filtered_readcounts_bedpe_tumor_sample_id, cobalt_filtered_readcounts_bedpe]
    }
}

def extractBafFilteredFiles( baf_filtered_files_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  baf_filtered_files_dir = baf_filtered_files_dir.tokenize().collect{"$it/*/*/*.baf.filtered.txt"}
  Channel
    .fromPath(baf_filtered_files_dir, type:'file')
    .ifEmpty { error "No .baf.* files found in ${baf_filtered_files_dir}." }
    .map { baf_filtered_file_path ->
        baf_filtered_file = baf_filtered_file_path
        baf_filtered_tumor_sample_id = baf_filtered_file_path.getName().toString().replaceAll(/.baf.filtered.txt$/, '')
        baf_filtered_normal_sample_id =  baf_filtered_file_path.getParent().getName()
        baf_filtered_donor_id =  baf_filtered_file_path.getParent().getParent().getName()
        [baf_filtered_donor_id, baf_filtered_normal_sample_id, baf_filtered_tumor_sample_id, baf_filtered_file]
    }
}

def extractBafBinnedFiles( baf_binned_files_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  baf_binned_files_dir = baf_binned_files_dir.tokenize().collect{"$it/*/*/*.baf.binned*.txt"}
  Channel
    .fromPath(baf_binned_files_dir, type:'file')
    .ifEmpty { error "No .baf.binned*.txt files found in ${baf_binned_files_dir}." }
    .map { baf_binned_file_path ->
        baf_binned_file = baf_binned_file_path
        baf_binned_tumor_sample_id = baf_binned_file_path.getName().toString().replaceAll(/.baf.binned(.*).txt$/, '')
        baf_binned_normal_sample_id =  baf_binned_file_path.getParent().getName()
        baf_binned_donor_id =  baf_binned_file_path.getParent().getParent().getName()
        [baf_binned_donor_id, baf_binned_normal_sample_id, baf_binned_tumor_sample_id, baf_binned_file]
    }
}

def extractBafSegmentsFiles( baf_segments_files_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  baf_segments_files_dir = baf_segments_files_dir.tokenize().collect{"$it/*/*/*.baf.segments.*"}
  Channel
    .fromPath(baf_segments_files_dir, type:'file')
    .ifEmpty { error "No .baf.segments.* files found in ${baf_segments_files_dir}." }
    .map { baf_segments_files_path ->
        baf_segments_file = baf_segments_files_path
        baf_segments_tumor_sample_id = baf_segments_files_path.getName().toString().replaceAll(/.baf.segments.*$/, '')
        baf_segments_normal_sample_id =  baf_segments_files_path.getParent().getName()
        baf_segments_donor_id =  baf_segments_files_path.getParent().getParent().getName()
        [baf_segments_donor_id, baf_segments_normal_sample_id, baf_segments_tumor_sample_id, baf_segments_file]
    }
}

def extractGripssFilteredFiles( gripss_filtered_files_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  gripss_filtered_files_dir = gripss_filtered_files_dir.tokenize().collect{"$it/*/*/*.svs.postfilter.{vcf,vcf.gz}"}
  Channel
    .fromPath(gripss_filtered_files_dir, type:'file')
    .ifEmpty { error "No .svs.postfilter.vcf(.gz) files found in ${gripss_filtered_files_dir}." }
    .map { gripss_filtered_file_path ->
        gripss_filtered_file = gripss_filtered_file_path
        gripss_filtered_tumor_sample_id = gripss_filtered_file_path.getName().toString().replaceAll(/(.svs.postfilter)*.vcf(.gz)*$/, '')
        gripss_filtered_normal_sample_id =  gripss_filtered_file_path.getParent().getName()
        gripss_filtered_donor_id =  gripss_filtered_file_path.getParent().getParent().getName()
        [gripss_filtered_donor_id, gripss_filtered_normal_sample_id+";"+gripss_filtered_tumor_sample_id, gripss_filtered_file]
    }
}

def extractIntegratedSvFiles( integrated_sv_files_dir ) {
  // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
  integrated_sv_files_dir = integrated_sv_files_dir.tokenize().collect{"$it/*/*/*.integrated.*"}
  Channel
    .fromPath(integrated_sv_files_dir, type:'file')
    .ifEmpty { error "No .integrated. files found in ${integrated_sv_files_dir}." }
    .map { integrated_sv_file_path ->
        integrated_sv_file = integrated_sv_file_path
        integrated_sv_tumor_sample_id = integrated_sv_file_path.getName().toString().replaceAll(/.integrated.*/, '')
        integrated_sv_normal_sample_id =  integrated_sv_file_path.getParent().getName()
        integrated_sv_donor_id =  integrated_sv_file_path.getParent().getParent().getName()
        [integrated_sv_donor_id, integrated_sv_normal_sample_id, integrated_sv_tumor_sample_id, integrated_sv_file]
    }
}
