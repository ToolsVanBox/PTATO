library(VariantAnnotation)
library(ggplot2)
library(MutationalPatterns)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
ptato_vcf_fname <- args[1]
ref_genome <- args[2]
out_vcf_fname = args[3]
library( ref_genome, character.only = TRUE )

ptato_vcf <- readVcf( ptato_vcf_fname )
chroms <- seqlevels(SeqinfoForBSGenome(genome = ref_genome))
if ( grepl('chr', chroms[1], ignore.case = T) & !grepl('chr', as.character(seqnames(ptato_vcf)[1]), ignore.case = T) ) {
  seqlevels(ptato_vcf) <- paste0("chr", seqlevels(ptato_vcf))
}
if ( !grepl('chr', chroms[1], ignore.case = T) & grepl('chr', as.character(seqnames(ptato_vcf)[1]), ignore.case = T) ) {
  seqlevels(ptato_vcf) <- gsub("chr", "", seqlevels(ptato_vcf))
}

ptato_vcf <- ptato_vcf[which(seqnames(ptato_vcf) %in% chroms),]
ptato_gr <- rowRanges(ptato_vcf)
PTAprob1 <- as.numeric(geno(ptato_vcf)$PTAprob[,,1])
PTAprob2 <- as.numeric(geno(ptato_vcf)$PTAprob[,,2])
PTAprob3 <- as.numeric(geno(ptato_vcf)$PTAprob[,,3])
PTAprob <- PTAprob1
PTAprob[which(PTAprob == 0 | is.na(PTAprob))] <- PTAprob2[which(PTAprob == 0 | is.na(PTAprob))]
PTAprob[which(PTAprob == 0 | is.na(PTAprob))] <- PTAprob3[which(PTAprob == 0 | is.na(PTAprob))]
ptato_gr$PTAprob <- PTAprob

# Makes a list with two mutation matrices. Each matrix contains the mutation profiles at different PTAprobsCutoffs. One matrix contains the mutations passing the filter (eg PTAprob of variant <= Cutoff), the other contains the mutations failing the filter (eg PTAprob of variant > Cutoff)
# Additionally, the original PTAprobCutoff as determined by PTATO is also included in the list
calc_cosim_mutmat <- function(gr, cutoffs =  seq(0.1, 0.8, 0.025)){
  GenomeInfoDb::genome(gr) <- 'hg38'
  seqlevels(gr) <- chroms
  grl_cosim <- list()
  grl_removed <- list()

  for(cutoff in cutoffs){
    #print(cutoff)
    grl_cosim[[as.character(cutoff)]] <- gr[which(gr$PTAprob <= cutoff)]
    grl_removed[[as.character(cutoff)]] <- gr[which(gr$PTAprob > cutoff)]
  }
  mut_mat <- mut_matrix(vcf_list = grl_cosim, ref_genome = ref_genome)
  mut_mat_removed <- mut_matrix(vcf_list = grl_removed, ref_genome = ref_genome)
  mut_mat_sample <- list(PASS = mut_mat, FAIL = mut_mat_removed, PTAprobsCutoff = unique(gr$PTAprobCutoff))
  return(mut_mat_sample)
}

mut_mat_cosim <- calc_cosim_mutmat(gr = ptato_gr)

cos_sim_matrix <- cos_sim_matrix(mut_mat_cosim[["PASS"]], mut_mat_cosim[["PASS"]])

### CLustering (code from MutationalPatterns)
hc.sample <- hclust(dist(cos_sim_matrix[7:nrow(cos_sim_matrix),
                                        7:ncol(cos_sim_matrix)]), method = "complete")

# This function can be used to get the samples in two clusters
sub_grp <- cutree(hc.sample, k = 2)
PTAprob_cutoff <- max(names(sub_grp[sub_grp == 1]))

myfilteredrows <- which(ptato_gr$PTAprob <= PTAprob_cutoff)

ptatofilter_vcf <- ptato_vcf[myfilteredrows,]

# meta(header(ptatofilter_vcf))$PTAprob_cutoff = mean(pred_rec$PTAprob_cutoff)
meta(header(ptatofilter_vcf))$PTAprobCutoff <- DataFrame("Value" = paste0('"', mean(pred_rec$PTAprob_cutoff), ' ', '"'), row.names = "PTAprobCutoff")

outvcf <- file(out_vcf_fname, open="a")
writeVcf(ptatofilter_vcf, outvcf)
