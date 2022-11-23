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
source(paste(sourcedir,"/functions/get_context.R",sep=""), chdir = T)

args = commandArgs(trailingOnly=TRUE)

# test if there is an argument: if not, return an error
if (length(args)==0) {
  stop("Input file must be supplied", call.=FALSE)
} else {
  vcf_file = args[1]
  out_file = args[2]
  ref_genome = args[3]
}

library( ref_genome, character.only = TRUE )

chroms <- seqlevels(SeqinfoForBSGenome(genome = ref_genome))
vcf = readVcf(vcf_file)

if ( length(vcf) == 0 ) {
  df <- data.frame()
  write.table(df, file=out_file, quote=F, sep="\t", row.names=F, col.names=F)
  quit()
}

if ( grepl('chr', chroms[1], ignore.case = T) & !grepl('chr', as.character(seqnames(vcf)[1]), ignore.case = T) ) {
  seqlevels(vcf) <- paste0("chr", seqlevels(vcf))
}
if ( !grepl('chr', chroms[1], ignore.case = T) & grepl('chr', as.character(seqnames(vcf)[1]), ignore.case = T) ) {
  seqlevels(vcf) <- gsub("chr", "", seqlevels(vcf))
}

vcf <- vcf[which(seqnames(vcf) %in% chroms),]
gr = granges(vcf)

get_context(gr, out_file)
