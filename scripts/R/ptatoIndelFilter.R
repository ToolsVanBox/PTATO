library(VariantAnnotation)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
ptato_vcf_fname <- args[1]
out_vcf_fname <- args[2]

ptato_vcf <- readVcf( ptato_vcf_fname )

ptato_gr <- granges(ptato_vcf)
ptato_gr_df <- as.data.frame(ptato_gr)
PTAprob1 <- as.numeric(geno(ptato_vcf)$PTAprob[,,1])
PTAprob2 <- as.numeric(geno(ptato_vcf)$PTAprob[,,2])
PTAprob3 <- as.numeric(geno(ptato_vcf)$PTAprob[,,3])
PTAprob <- PTAprob1

ptato_gr_df$PTAprob <- PTAprob

PTAprob_cutoff <- 0.5
PTAprob_cutoff_conf <- "NA - NA"
myfilteredrows <- which(ptato_gr_df$PTAprob <= PTAprob_cutoff)

ptatofilter_vcf <- ptato_vcf[myfilteredrows,]

meta(header(ptatofilter_vcf))$PTAprobCutoff <- DataFrame("Value" = paste0('"', PTAprob_cutoff, ' (', PTAprob_cutoff_conf, ')', '"'), row.names = "PTAprobCutoff")

outvcf <- file(out_vcf_fname, open="a")
writeVcf(ptatofilter_vcf, outvcf)
