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
library(StructuralVariantAnnotation)

source(paste(sourcedir,"/functions/get_circos_config.R",sep=""), chdir = T)
options(scipen = 999)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
CNVs_file <- args[1]
ReadCounts_1mb_file <- args[2]
BAF_1mb_file <- args[3]
ReadCounts_segments_file <- args[4]
BAF_segments_file <- args[5]
SVs_file <- args[6]
OUTPUT_dir <- args[7]
Circos_Path <- args[8]

Make_SV_Circos(CNVs_file = CNVs_file,
               ReadCounts_1mb_file = ReadCounts_1mb_file,
               BAF_1mb_file = BAF_1mb_file,
               ReadCounts_segments_file = ReadCounts_segments_file,
               BAF_segments_file = BAF_segments_file,
               SVs_file = SVs_file,
               OUTPUT_dir = OUTPUT_dir,
               Circos_Path = Circos_Path)
