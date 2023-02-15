# PTATO
PTa Analysis TOolkit

Includes QC and SNV, INDEL, SV and CNV filtration.

## Dependencies

- Nextflow
- Singularity
- R/4.1.2

### R libraries
- BSgenome
- copynumber
- cowplot
- ggplot2
- gtools
- MutationalPatterns
- randomForest
- scales
- StructuralVariantAnnotation
- VariantAnnotation

## Installing & Setup

1. [Install Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
2. [Install Singularity](https://sylabs.io/guides/3.5/admin-guide/)
3. [Pull/Clone PTATO](#pull-or-clone)

## Pull or Clone 
```
git clone git@github.com:ToolsVanBox/PTATO.git
```

## Resources
First extract the following resources files:
- resources/hg38/gripss/gridss_pon_breakpoint.tar.gz
- resources/hg38/cobalt/COBALT_PTA_Normalized_Full.tar.gz
- resources/hg38/smurf/Mutational_blacklists/Fetal_15x_raw_variants_hg38.tar.gz
- resources/hg38/smurf/Mutational_blacklists/MSC_healthyBM_raw_variants_hg38.tar.gz

### Reference genome
Please download the reference genome fasta file. Must have the following files:
- *.dict
- *.fasta
- *.fasta.amb
- *.fasta.ann
- *.fasta.bwt
- *.fasta.dict
- *.fasta.fai
- *.fasta.gridsscache
- *.fasta.img
- *.fasta.pac
- *.fasta.sa
- *.full.len
- *.genomesize.txt
- *.len
- *.len.genome

And put them in this folder `resources/hg38`

### SHAPEIT resources
Can you download the SHAPEIT resources files here:
- http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
- https://github.com/odelaneau/shapeit4/blob/master/maps/genetic_maps.b38.tar.gz

And put them in this folder, respectively:
- resources/hg38/shapeit/Phasing_reference/
- resources/hg38/shapeit/shapeit_maps/


## Running the workflow
In this section we'll provide you with a way to run the workflow.

### Change the run-template.config to start your analysis

Always keep these lines in your run.config file:
```
includeConfig "${projectDir}/configs/process.config"
includeConfig "${projectDir}/configs/nextflow.config"
includeConfig "${projectDir}/configs/resources.config"
```
All of the parameters in the params section can also be supplied on the commandline or can be pre-filled in the run.config file.
```
params {

  run {
    snvs = true
    QC = false
    svs = false
    indels = false
    cnvs = false
  }

  // TRAINING
  train {
    version = '2.0.0'
  }
  pta_vcfs_dir = ''
  nopta_vcfs_dir = ''
  // END TRAINING

  // TESTING
  input_vcfs_dir = ''
  bams_dir = ''
  // END TESTING

  out_dir = ''
  bulk_names = [
    ['donor_id', 'sample_id'],
  ]

  snvs {
    rf_rds = "${projectDir}/resources/hg38/snvs/randomforest/randomforest_v1.0.0.rds"
  }

  indels {
    rf_rds = ''
    excludeindellist = "${projectDir}/resources/hg38/indels/excludeindellist/PTA_Indel_ExcludeIndellist_normNoGTrenamed.vcf.gz"
  }
  optional {

    germline_vcfs_dir = ''

    short_variants {
      somatic_vcfs_dir = ''
      walker_vcfs_dir = ''
      phased_vcfs_dir = ''
      ab_tables_dir = ''
      context_beds_dir = ''
      features_beds_dir = ''
    }

    snvs {
      rf_tables_dir = ''
      ptato_vcfs_dir = ''
    }

    indels {
      rf_tables_dir = ''
      ptato_vcfs_dir = ''
    }

    qc {
      wgs_metrics_dir = ''
      alignment_summary_metrics_dir = ''
    }

    svs {
      gridss_driver_vcfs_dir = ''
      gridss_unfiltered_vcfs_dir = ''
      gripss_somatic_filtered_vcfs_dir = ''
      gripss_filtered_files_dir = ''
      integrated_sv_files_dir = ''
    }

    cnvs {
      cobalt_ratio_tsv_dir = ''
      cobalt_filtered_readcounts_dir = ''
      baf_filtered_files_dir = ''
    }
  }


}
```
### Starting the full workflow
Create the run.config file to look like this:
```
params {

  run {
    snvs = true
    QC = true
    svs = true
    indels = true
    cnvs = true
  }

  // TRAINING
  train {
    version = '2.0.0'
  }
  pta_vcfs_dir = ''
  nopta_vcfs_dir = ''
  // END TRAINING

  // TESTING
  input_vcfs_dir = '/path/to/vcfs_dir/'
  bams_dir = '/path/to/bams_dir/'
  // END TESTING

  out_dir = ''
  bulk_names = [
    ['Donor_1', 'mycontrol'],
  ]

  snvs {
    rf_rds = "${projectDir}/resources/hg38/snvs/randomforest/randomforest_v1.0.0.rds"
  }

  indels {
    rf_rds = ''
    excludeindellist = "${projectDir}/resources/hg38/indels/excludeindellist/PTA_Indel_ExcludeIndellist_normNoGTrenamed.vcf.gz"
  }
  optional {

    germline_vcfs_dir = ''

    short_variants {
      somatic_vcfs_dir = ''
      walker_vcfs_dir = ''
      phased_vcfs_dir = ''
      ab_tables_dir = ''
      context_beds_dir = ''
      features_beds_dir = ''
    }

    snvs {
      rf_tables_dir = ''
      ptato_vcfs_dir = ''
    }

    indels {
      rf_tables_dir = ''
      ptato_vcfs_dir = ''
    }

    qc {
      wgs_metrics_dir = ''
      alignment_summary_metrics_dir = ''
    }

    svs {
      gridss_driver_vcfs_dir = ''
      gridss_unfiltered_vcfs_dir = ''
      gripss_somatic_filtered_vcfs_dir = ''
      gripss_filtered_files_dir = ''
      integrated_sv_files_dir = ''
    }

    cnvs {
      cobalt_ratio_tsv_dir = ''
      cobalt_filtered_readcounts_dir = ''
      baf_filtered_files_dir = ''
    }
  }


}
```

Run the workflow on slurm :
```
nextflow run ptato.nf -c run.config --out_dir /processed_data/ptato/ -profile slurm -resume
```

## Input folder structure
```
/path/to/vcfs_dir
  ./Donor_1
    ./myfile.vcf(.gz)
    
/path/to/bams_dir
  ./Donor_1
    ./mycontrol.bam
    ./mysample1.bam
    ./mysample2.bam
    ...
```

