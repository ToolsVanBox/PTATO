add_ptaprob_sample = function(filt_vcf_fname){

    # Read vcf with pta probabilities
    filt_vcf <- readVcf(filt_vcf_fname)
    sample_name = samples(header(filt_vcf))

    # Get pta probabilities
    pta_val_l = lapply(1:3, function(x){
        geno(filt_vcf)$PTAprob[,,x]
    })
    pta_val_combined_l = mapply(c, pta_val_l[[1]], pta_val_l[[2]], pta_val_l[[3]], SIMPLIFY = FALSE)

    # Turn pta probabilities into a matrix
    m = as.matrix(pta_val_combined_l)
    colnames(m) = sample_name

    # Add PTA probabilities to main vcf
    sample_i <- which(rownames(input_vcf) %in% rownames(filt_vcf))
    geno(input_vcf)$PTAprob[sample_i, sample_name] <<-  m
    return(0)
}
