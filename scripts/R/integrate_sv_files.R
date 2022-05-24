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
source(paste(sourcedir,"/functions/SV_Integration.R",sep=""), chdir = T)
options(scipen = 999) 

# Get arguments
args = commandArgs(trailingOnly=TRUE)
BAFs_file <- args[1]
ReadCounts_file <- args[2]
BAF_segments.bedpe <- args[3]
ReadCounts_segments.bedpe <- args[4]
Breakends_file <- args[5]
OUTPUT_dir <- args[6]
SAMPLE <- args[7]

BAFs <- read.delim(BAFs_file,header=T)
ReadCounts <- read.delim(ReadCounts_file, header=T)
Breakends <- readVcf(Breakends_file)

BAF_segments_g <- pairs2breakpointgr(rtracklayer::import(BAF_segments.bedpe, format = "bedpe"))
ReadCounts_segments_g <- pairs2breakpointgr(rtracklayer::import(ReadCounts_segments.bedpe, format = "bedpe"))
BAF_g <- GRanges(seqnames = BAFs$chromosome, IRanges(start = BAFs$position, end = BAFs$position),
                 BAF = BAFs$VAF_normalized)

ReadCounts_f <- ReadCounts[ReadCounts$FILTER == "PASS",]
ReadCounts_f$CopyNumber <- ReadCounts_f$filteredReadCount*2
ReadCounts_g <- GRanges(seqnames = ReadCounts_f$chromosome[which(ReadCounts_f$FILTER == "PASS")],
                    IRanges(start = ReadCounts_f$position[which(ReadCounts_f$FILTER == "PASS")],
                            end = ReadCounts_f$position[which(ReadCounts_f$FILTER == "PASS")]),
                    RD = ReadCounts_f$CopyNumber)

CNVs <- Merge_CNVs(BAF = BAF_g,
                   RD = ReadCounts_g,
                   BAF_segments = BAF_segments_g,
                   RD_segments = ReadCounts_segments_g,
                   SVs = Breakends, Output_dir = OUTPUT_dir,
                   SAMPLE = SAMPLE)
