library(VariantAnnotation)

# Get arguments
args = commandArgs(trailingOnly=TRUE)
ptato_vcf_fname <- args[1]
walker_vcf_fname = args[2]
out_vcf_fname = args[3]
out_table_fname = args[4]

ptato_vcf <- readVcf( ptato_vcf_fname )
walker_vcf <- readVcf( walker_vcf_fname )

ptato_gr <- granges(ptato_vcf)
ptato_gr_df <- as.data.frame(ptato_gr)
PTAprob1 <- as.numeric(geno(ptato_vcf)$PTAprob[,,1])
PTAprob2 <- as.numeric(geno(ptato_vcf)$PTAprob[,,2])
PTAprob3 <- as.numeric(geno(ptato_vcf)$PTAprob[,,3])
PTAprob <- PTAprob1
PTAprob[which(PTAprob == 0)] <- PTAprob2[which(PTAprob == 0)]
PTAprob[which(PTAprob == 0)] <- PTAprob3[which(PTAprob == 0)]

ptato_gr_df$PTAprob <- PTAprob

walker_gr <- granges(walker_vcf)
walker_gr_df <- as.data.frame(walker_gr)
walker_gr_df$WS <- as.numeric(geno(walker_vcf)$WS)

merged_df <- merge(ptato_gr_df, walker_gr_df, by=c("seqnames","start"), sort = FALSE)
merged_df_scores <- na.omit(merged_df[,c("seqnames","start","PTAprob","WS")])

true_ws <- 1000
false_ws <- 1

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

maxf1 <- myresult[which(myresult$p == min(myresult[which(myresult$F1 == max(na.omit(myresult$F1))),]$p)),]
myresult2 <- myresult[which( abs(myresult$recall-myresult$precision) > 0 ),]
pred_rec <- myresult2[which(abs(myresult2$recall-myresult2$precision) == min(na.omit(abs(myresult2$recall-myresult2$precision)))),]

# myresult2 <- myresult[which( abs(myresult$specificity-myresult$accuracy) > 0 ),]
# pred_rec <- myresult2[which(abs(myresult2$specificity-myresult2$accuracy) == min(na.omit(abs(myresult2$specificity-myresult2$accuracy)))),]

low_ws_rownames <- rownames(walker_gr_df[which(walker_gr_df$WS < 1),])
myfilteredrows <- which(ptato_gr_df$PTAprob <= mean(pred_rec$PTAprob_cutoff))

ptatofilter_vcf <- ptato_vcf[myfilteredrows,]
ptatofilter_vcf[!rownames(ptatofilter_vcf) %in% low_ws_rownames,]

# meta(header(ptatofilter_vcf))$PTAprob_cutoff = mean(pred_rec$PTAprob_cutoff)
meta(header(ptatofilter_vcf))$PTAprobCutoff <- DataFrame("Value" = paste0('"', mean(pred_rec$PTAprob_cutoff), ' ', '"'), row.names = "PTAprobCutoff")

outvcf <- file(out_vcf_fname, open="a")
writeVcf(ptatofilter_vcf, outvcf)
