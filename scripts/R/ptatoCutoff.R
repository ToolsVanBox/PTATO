# This script is used to calculate sample-specific PTA probability cutoffs
# The final cutoffs are determined based on two complementary methods: the linked read (walker) cutoff and the cosine similarity cutoff

library(VariantAnnotation)
library(MutationalPatterns)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
ptato_vcf_fname <- args[1]
walker_vcf_fname = args[2]
ref_genome = args[3]
prec_cutoff = args[4]
out_table_fname = args[5]

library( ref_genome, character.only = TRUE )
genome <- gsub(ref_genome, pattern = ".*\\.", replacement = "") # assuming the ref_genome is BSgenome object (eg BSgenome.Hsapiens.UCSC.hg38)

ptato_vcf <- readVcf( ptato_vcf_fname, genome = genome)
walker_vcf <- readVcf( walker_vcf_fname, genome = genome)

ptato_gr <- granges(ptato_vcf)
ptato_gr_df <- as.data.frame(ptato_gr)
# For some variants, not all RF features can be determined. Therefore, PTATO calculates three PTAprobs for each variants: 1 =  RF with all features, 2 = RF without allelic imbalance and 3 = RF without repliseq and allelic imbalance. If PTAprob 1 cannot be determined, PTAprob 2 will be used (and PTAprob 3 if 2 can also not be determined)
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

merged_df <- merge(ptato_gr_df, walker_gr_df, by=c("seqnames","start"), sort = FALSE)
merged_df_scores <- na.omit(merged_df[,c("seqnames","start","PTAprob","WS")])

true_ws <- 1000
false_ws <- 1

# the linked read (walker) analysis determined a small set of true and false positive variants
# a range of PTAprob cutoff values is made, and for each cutoff the number of included and excluded linked-read based true and false positives is calculated (with precision, recall, specificity and accuracy)
# the cutoff with the smallest difference between precision and recall is taken as the final linked-read (walker) cutoff
myresult <- matrix(ncol=10)[-1,]
for (p in seq(0,1,0.01)) {
  TP <- nrow(merged_df_scores[merged_df_scores$WS >= true_ws & merged_df_scores$PTAprob <= p,])
  FP <- nrow(merged_df_scores[merged_df_scores$WS < false_ws & merged_df_scores$PTAprob <= p,])
  FN <- nrow(merged_df_scores[merged_df_scores$WS >= true_ws & merged_df_scores$PTAprob > p,])
  TN <- nrow(merged_df_scores[merged_df_scores$WS < false_ws & merged_df_scores$PTAprob > p,])
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  accuracy <- (TP + TN) / (TP + FP + FN + TN)
  F1 = 2 * (precision * recall) / (precision + recall)
  myresult <- rbind(myresult, c(p,TP,FP,FN,TN,precision,recall,specificity,accuracy,F1))
}
colnames(myresult) <- c("PTAprob_cutoff","TP","FP","FN","TN","precision","recall","specificity","accuracy","F1")
myresult <- as.data.frame(myresult)
write.table( myresult, out_table_fname, quote=F, col.names=T, row.names=F)

if (nrow(merged_df) == 0) {
  cat("NA", "NA")
} else {
  maxf1 <- myresult[which(myresult$p == min(myresult[which(myresult$F1 == max(na.omit(myresult$F1))),]$p)),]
  if (nrow(maxf1) == 0) {
    cat("NA", "NA")
  } else {
    myresult2 <- myresult[which( abs(myresult$recall-myresult$precision) > 0 ),]
    prec_rec <- myresult2[which(abs(myresult2$recall-myresult2$precision) == min(na.omit(abs(myresult2$recall-myresult2$precision)))),]
    # the prec-recall can be the same at multiple cutoffs, therefore take the mean of these cutoffs
    PTAprob_cutoff_walker <- mean(prec_rec$PTAprob_cutoff)
    if ( length(grep("indels.ptato",ptato_vcf_fname)) > 0 ) {
      PTAprob_cutoff <- c(as.numeric(PTAprob_cutoff_walker),'NA')
      cat( PTAprob_cutoff_walker, PTAprob_cutoff )
    } else {
      chroms <- seqlevels(SeqinfoForBSGenome(genome = ref_genome))
      if ( grepl('chr', chroms[1], ignore.case = T) & !grepl('chr', as.character(seqnames(ptato_vcf)[1]), ignore.case = T) ) {
        seqlevels(ptato_vcf) <- paste0("chr", seqlevels(ptato_vcf))
      }
      if ( !grepl('chr', chroms[1], ignore.case = T) & grepl('chr', as.character(seqnames(ptato_vcf)[1]), ignore.case = T) ) {
        seqlevels(ptato_vcf) <- gsub("chr", "", seqlevels(ptato_vcf))
      }

      # the chroms in the ptato_gr should be the same as the chroms in the ref_genome
      ptato_vcf <- ptato_vcf[which(seqnames(ptato_vcf) %in% chroms),]
      ptato_gr <- rowRanges(ptato_vcf)
      seqlevels(ptato_gr) <- seqlevels(ptato_gr)[seqlevels(ptato_gr) %in% chroms]
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

      # Calculate the cosine similarities between all the mutation profiles (of the mutations < than that specific cutoff) at each different cutoff
      cos_sim_matrix <- cos_sim_matrix(mut_mat_cosim[["PASS"]], mut_mat_cosim[["PASS"]])

      # Cluster the mutation profiles of each different cutoff (code from MutationalPatterns)
      # The first profiles are excluded, because they contain too few variants
      hc.sample <- hclust(dist(cos_sim_matrix[8:nrow(cos_sim_matrix),
                                              8:ncol(cos_sim_matrix)]), method = "complete")

      # Divide the cosine similarity matrix in four clusters (1 = true pos (low number, low PTAprobs), 2 = mostly true pos, 3 = mix of true+false pos, 4 = false pos, high PTAprobs)
      sub_grp <- cutree(hc.sample, k = 4)
      # The mean value of the third cluster (containing a mix of true and false positives) is taken as the Cosine Similarity cutoff
      PTAprob_cutoff_cluster <- mean(as.numeric(names(sub_grp[sub_grp == 3])))

      # Calculate the final PTAprob cutoff.
      # If the precision-recall of the linked-read (walker) cutoff is too low (< prec_cutoff, which is normally 0.5), only the Cosine Similarity cutoff is taken as the final cutoff. Else the mean of the linked-read and cosine similarity cutoff is taken as the final PTAprob cutoff
      if ( mean(prec_rec$precision) < prec_cutoff ) {
        PTAprob_cutoff_conf <- c(PTAprob_cutoff_walker,PTAprob_cutoff_cluster)
        PTAprob_cutoff <- PTAprob_cutoff_cluster
      } else {
        PTAprob_cutoff_conf <- as.numeric(c(PTAprob_cutoff_walker,PTAprob_cutoff_cluster))
        PTAprob_cutoff <- mean(PTAprob_cutoff_conf)
      }

      cat( PTAprob_cutoff, PTAprob_cutoff_conf )
    }
  }
}
