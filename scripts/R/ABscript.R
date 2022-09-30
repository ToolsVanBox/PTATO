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
source(paste(sourcedir,"/functions/create_empty_ab_table.R",sep=""), chdir = T)
source(paste(sourcedir,"/functions/AB_somatic.R",sep=""), chdir = T)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
somatic_vcf_fname <- args[1]
germline_vcf_fname <- args[2]
phased_germline_vcf_fname <- args[3]
chrom <- args[4]
output_file <- args[5]
ref_genome <- args[6]
flank <- 200000

somatic_vcf <- readVcf(somatic_vcf_fname)
somatic_vcf_chr <- somatic_vcf[seqnames(somatic_vcf) == chrom,]

if ( length(somatic_vcf_chr) == 0 ) {
  ab_table <- create_empty_ab_table()
} else {
  ab_table <- list()
  for ( i in c(1:length(somatic_vcf_chr))) {
    ab_table_tmp <- AB_somatic( somatic_vcf_chr[i,] )
    if ( length(ab_table_tmp) == 0 ) {
      ab_table_tmp <- create_empty_ab_table()
    } else if ( ab_table_tmp == 0 ) {
      ab_table_tmp <- create_empty_ab_table()
    } else if ( ncol(ab_table_tmp) <= 1 ) {
      ab_table_tmp <- create_empty_ab_table()
    }
    ab_table[[i]] <- ab_table_tmp
  }
  ab_table <- do.call(rbind, ab_table)
  if ( length(ab_table) == 0 ) {
    ab_table <- create_empty_ab_table()
  }
}
write.table(ab_table, output_file, quote = FALSE, col.names = TRUE, row.names = FALSE)
