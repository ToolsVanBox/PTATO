library(VariantAnnotation)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
indels_vcf_fname <- args[1]
context_fname <- args[2]
excludeindellist_fname <- args[3]
out_vcf_fname <- args[4]

indels_vcf <- readVcf(indels_vcf_fname)
rownames(indels_vcf) <- paste(as.character(seqnames(indels_vcf)),start(indels_vcf),sep=":")

excludelist_vcf <- readVcf(excludeindellist_fname)
rownames(excludelist_vcf) <- paste(as.character(seqnames(excludelist_vcf)),start(excludelist_vcf),sep=":")

context_df <- read.table(context_fname)
PTAprob <- rep(0,nrow(indels_vcf))
PTAprob[which(rownames(indels_vcf) %in% rownames(excludelist_vcf))] <- 1

muttype_df <- context_df[,1:3]
colnames(muttype_df) <- c('seqnames','start','end')
muttype_df$muttype <- gsub(".+\\((.+)\\).+","\\1",context_df$V4)

indels_vcf_df <- as.data.frame(granges(indels_vcf))
indels_vcf_df$start <- indels_vcf_df$start-1
indels_vcf_df <- merge(indels_vcf_df,muttype_df,by=c('seqnames','start','end'), all.x = TRUE, sort = FALSE)
PTAprob[which(indels_vcf_df$muttype == "C_insertion-5+" | indels_vcf_df$muttype == "T_insertion-5+")] <- 1

PTAprob <- as.list(as.data.frame(mapply(c,PTAprob,PTAprob,PTAprob)))

meta(header(indels_vcf))$PTArandomforest <- DataFrame("Value" = paste0('"', args[1], ' ', '"'), row.names = "PTArandomforest")
geno(header(indels_vcf))["PTAprob",] = list("3","Float","PTA probability values")
geno(indels_vcf)$PTAprob <- matrix(PTAprob,ncol=1)

writeVcf(indels_vcf, out_vcf_fname)
