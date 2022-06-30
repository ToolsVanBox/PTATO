# Script to add pta probabilities from single-sample vcfs to a single multi-sample smurf vcf.
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
source(paste(sourcedir,"/functions/add_ptaprob_sample.R",sep=""), chdir = T)
options(scipen = 999)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
input_vcf_fname <- args[1]
ptato_vcf_fnames <- args[2]
out_vcf_fname <- args[3]

# Read multi-sample smurf vcf.
input_vcf <- readVcf(input_vcf_fname)

# Create empty list for one sample
empty_pta = c(".", ".", ".")
empty_l = lapply(rownames(input_vcf), function(x){return(empty_pta)})
names(empty_l) = rownames(input_vcf)

# Replicate the list for all samples and turn into a matrix
samples = samples(header(input_vcf))
empty_l_l = lapply(seq_along(samples), function(x){return(empty_l)})
empty_m = do.call(cbind, empty_l_l)
colnames(empty_m) = samples

# Set empty pta probabilities for vcf
geno(header(input_vcf))["PTAprob",] = list("3","Float","PTA probability values")
geno(input_vcf)$PTAprob <- empty_m

# Add pta probabilities per sample
trash = lapply(strsplit(ptato_vcf_fnames,",")[[1]], add_ptaprob_sample)

# Write final vcf
writeVcf(input_vcf, out_vcf_fname)
