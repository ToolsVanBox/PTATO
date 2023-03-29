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
library(copynumber)
source(paste(sourcedir,"/functions/Filter_COBALT.R",sep=""), chdir = T)
source(paste(sourcedir,"/functions/df_to_bedpe.R",sep=""), chdir = T)
options(scipen = 999)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
COBALT_file <- args[1]
CENTROMERE_file <- args[2]
CYTOBAND_file <- args[3]
COBALT_PON_FILE <- args[4]
ref_genome <- args[5] # Name of BSgenome package (eg "BSgenome.Hsapiens.UCSC.hg38")
OUTPUT_dir <- args[6]
Centromere_extension <- as.numeric(args[7]) # region around the centromeres that should be excluded (in bp)
COBALT_MIN_COV <- as.numeric(args[8]) # the lower percentile of readcount bins that should be excluded
COBALT_MAX_COV <- as.numeric(args[9]) # the upper percentile of readcount bins that should be excluded

library( ref_genome, character.only = TRUE )

if(identical(COBALT_file, character(0))){
  stop(paste("!!! COBALT_file not found"))
} else if (file.exists(COBALT_file) == FALSE) {
  stop("!!! COBALT_file not found")
} else if ( file.exists(COBALT_file) == TRUE) {
  print(paste("# COBALT_file = ", COBALT_file, sep = ""))
}

Genome <- seqlengths(SeqinfoForBSGenome(genome = ref_genome))
names(Genome) <- gsub(pattern = "chr", replacement = "", names(Genome))
Genome <- Genome[nchar(names(Genome)) < 3] # Remove the decoy chromosomes

# Read the centromere and cytoband data (used for both ReadCouns and BAF filtering)
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

print(paste("# Reading: ", COBALT_PON_FILE, sep = ""))
COBALT_PON <- read.delim(COBALT_PON_FILE, check.names = F)

print(paste("# Reading: ", COBALT_file, sep = ""))
COBALT_raw <- read.delim(COBALT_file, check.names = F)

ReadCounts <- FilterCobalt(COBALT_input = COBALT_raw, COBALT_PON = COBALT_PON,
                                      OUTPUT_dir = OUTPUT_dir,
                     CENTROMERES = CENTROMERES_g,
                     CYTOBAND = CYTOBAND_g,
                     COBALT_MIN_COV = COBALT_MIN_COV,
                     COBALT_MAX_COV = COBALT_MAX_COV,
                     CHR_SIZES = Genome)
