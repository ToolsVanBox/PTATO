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
library(StructuralVariantAnnotation)
source(paste(sourcedir,"/functions/Filter_SV_Breakends.R",sep=""), chdir = T)
options(scipen = 999)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
SV_VCF_file <- args[1]
BAFs_file <- args[2]
ReadCounts_file <- args[3]
GRIPSS_PON_FILE <- args[4]
OUTPUT_dir <- args[5]
max_dist <- as.numeric(args[6])
MaxBreakendCov <- as.numeric(args[7])

GRIPSS_PON <- read.delim(GRIPSS_PON_FILE)
BAFs <- read.delim(BAFs_file,header=T)
ReadCounts <- read.delim(ReadCounts_file, header=T)

Breakends <- Filter_Breakends(GRIPSS_PON = GRIPSS_PON,
                 SV_VCF_file = SV_VCF_file,
                 BAFs = BAFs,
                 ReadCounts = ReadCounts,
                 Output_dir = OUTPUT_dir,
                 max_dist = max_dist,
                 MaxBreakendCov = MaxBreakendCov,
                 max_SV_size_annotation = 10e6)
