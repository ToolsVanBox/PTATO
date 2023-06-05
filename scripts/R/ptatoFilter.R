library(VariantAnnotation)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
ptato_vcf_fname <- args[1]
walker_vcf_fname <- args[2]
ptaprob_cutoff_fname <- args[3]
out_vcf_fname <- args[4]

ptato_vcf <- readVcf( ptato_vcf_fname )
walker_vcf <- readVcf( walker_vcf_fname )

ptato_gr <- granges(ptato_vcf)
ptato_gr_df <- as.data.frame(ptato_gr)
PTAprob1 <- as.numeric(geno(ptato_vcf)$PTAprob[,,1])
PTAprob2 <- as.numeric(geno(ptato_vcf)$PTAprob[,,2])
PTAprob3 <- as.numeric(geno(ptato_vcf)$PTAprob[,,3])
PTAprob <- PTAprob1
PTAprob[which(PTAprob == 0 | is.na(PTAprob))] <- PTAprob2[which(PTAprob == 0 | is.na(PTAprob))]
PTAprob[which(PTAprob == 0 | is.na(PTAprob))] <- PTAprob3[which(PTAprob == 0 | is.na(PTAprob))]

ptato_gr_df$PTAprob <- PTAprob

walker_gr <- granges(walker_vcf)
walker_gr_df <- as.data.frame(walker_gr)
walker_gr_df$WS <- as.numeric(geno(walker_vcf)$WS)

PTAprob_cutoff_values <- read.table(ptaprob_cutoff_fname)
PTAprob_cutoff <- as.numeric(PTAprob_cutoff_values[1])
if (is.na(PTAprob_cutoff)) {
  PTAprob_cutoff <- 0.5
}
PTAprob_cutoff_conf <- paste(as.numeric(PTAprob_cutoff_values[2]),'-',as.numeric(PTAprob_cutoff_values[3]))

low_ws_rownames <- rownames(walker_gr_df[which(walker_gr_df$WS < 1),])
myfilteredrows <- which(ptato_gr_df$PTAprob <= PTAprob_cutoff)

ptatofilter_vcf <- ptato_vcf[myfilteredrows,]
ptatofilter_vcf[!rownames(ptatofilter_vcf) %in% low_ws_rownames,]

meta(header(ptatofilter_vcf))$PTAprobCutoff <- DataFrame("Value" = paste0('"', PTAprob_cutoff, ' (', PTAprob_cutoff_conf, ')', '"'), row.names = "PTAprobCutoff")

outvcf <- file(out_vcf_fname, open="a")
writeVcf(ptatofilter_vcf, outvcf)
