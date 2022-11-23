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

library(ggplot2)
library(StructuralVariantAnnotation)
source(paste(sourcedir,"/functions/SV_Plots.R",sep=""), chdir = T)
options(scipen = 999)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
ReadCounts_100kb_file <- args[1]
ReadCounts_1mb_file <- args[2]
ReadCounts_segments_file <- args[3]
BAF_binned_file <- args[4]
BAF_segments_file <- args[5]
CNV_file <- args[6]
OUTPUT_dir <- args[7]
SAMPLE <- args[8]

Generate_Plots(ReadCounts_100kb_file = ReadCounts_100kb_file,
               ReadCounts_1mb_file = ReadCounts_1mb_file,
               ReadCounts_segments_file = ReadCounts_segments_file,
               BAF_binned_file = BAF_binned_file,
               BAF_segments_file = BAF_segments_file,
               CNV_file = CNV_file,
               OUTPUT_dir = OUTPUT_dir,
               SAMPLE = SAMPLE)
