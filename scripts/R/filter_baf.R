.getSourceDir <- function() {
  cmdArgs <- commandArgs(trailingOnly=FALSE)
  fileArg <- "--file="
  match <- grep(fileArg, cmdArgs)
  if (length(match) > 0) {
      sourcedir <- dirname(gsub(fileArg, "", cmdArgs[match]))
  }
  return( sourcedir )
}
sourcedir <- .getSourceDir()

library(VariantAnnotation)
library(gtools)
library(copynumber)
source(paste(sourcedir,"/functions/Filter_BAF.R",sep=""), chdir = T)
source(paste(sourcedir,"/functions/df_to_bedpe.R",sep=""), chdir = T)
options(scipen = 999) 

# Get arguments
args = commandArgs(trailingOnly=TRUE)
SNV_VCF <- args[1]
SAMPLE <- args[2]
NORMAL <- args[3]
CENTROMERE_file <- args[4]
CYTOBAND_file <- args[5]
ref_genome <- args[6]
OUTPUT_dir <- args[7]
Centromere_extension <- as.numeric(args[8])
MIN_DEPTH_SAMPLE <- as.numeric(args[9])
MIN_DEPTH_NORMAL <- as.numeric(args[10])

library( ref_genome, character.only = TRUE )

Genome <- seqlengths(SeqinfoForBSGenome(genome = ref_genome))
names(Genome) <- gsub(pattern = "chr", replacement = "", names(Genome))
Genome <- Genome[nchar(names(Genome)) < 3] # Remove the decoy chromosomes

CENTROMERES <- read.table(CENTROMERE_file)
CENTROMERES_g <- GRanges(seqnames = gsub("chr", "", CENTROMERES[,2]), IRanges(start = CENTROMERES[,3], end = CENTROMERES[,4]))
# Extend centromeres to remove bins around centromeres
start(CENTROMERES_g) <- start(CENTROMERES_g) - Centromere_extension
end(CENTROMERES_g) <- end(CENTROMERES_g) + Centromere_extension
CENTROMERES_g <- reduce(CENTROMERES_g) # merge overlapping positions

if(file.exists(CYTOBAND_file) == T){
  CYTOBAND <- read.delim(CYTOBAND_file, header = FALSE,  col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
  CYTOBAND_g <- GRanges(seqnames = gsub("chr", "", CYTOBAND$chrom), IRanges(start = CYTOBAND$chromStart, end =CYTOBAND$chromEnd), arm = ifelse(grepl(pattern = "p", CYTOBAND$name) == T, "p", "q"))
}

BAFs <- Filter_BAF(SNV_VCF_File = SNV_VCF, OUTPUT_dir = OUTPUT_dir,
                   SAMPLE = SAMPLE, NORMAL = NORMAL,
                   CENTROMERES = CENTROMERES_g,
                   CYTOBAND = CYTOBAND_g,
                   CHR_SIZES = Genome,
                   MIN_DEPTH_SAMPLE = MIN_DEPTH_SAMPLE,
                   MIN_DEPTH_NORMAL = MIN_DEPTH_NORMAL)
