library(BSgenome.Hsapiens.UCSC.hg38)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(VariantAnnotation)
args = commandArgs(trailingOnly=TRUE)

# test if there is an argument: if not, return an error
if (length(args)==0) {
  stop("Input file must be supplied", call.=FALSE)
} else {
  vcf_file = args[1]
  out_file = args[2]
}
#' Get indel contexts
#'
#' @details
#' Determines the COSMIC context from a GRanges or GRangesList object containing Indel mutations.
#' It applies the get_indel_context_gr function to each gr in the input.
#' It searches for repeat units both to the left and right of the indel.
#'
#' @param vcf_list GRanges or GRangesList object containing Indel mutations.
#' The mutations should be called similarly to HaplotypeCaller.
#' @param ref_genome BSgenome reference genome object
#'
#' @return A modified version of the input grl. In each gr two columns have been added.
#' "muttype" showing the main indel type and "muttype_sub" which shows the subtype.
#' The subtype is either the number of repeats or the microhomology length.
#'
#' @examples
#'
#' ## Get a GRangesList or GRanges object with only indels.
#' ## See 'read_vcfs_as_granges' or 'get_mut_type' for more info on how to do this.
#' indel_grl <- readRDS(system.file("states/blood_grl_indel.rds",
#'   package = "MutationalPatterns"
#' ))
#'
#' ## Load the corresponding reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' ## Get the indel contexts
#' get_indel_context(indel_grl, ref_genome)
#' @family Indels
#' @importFrom magrittr %>%
#' @seealso
#' \code{\link{read_vcfs_as_granges}}, \code{\link{get_mut_type}}
#'
#' @export
#'
get_indel_context <- function(vcf_list, ref_genome) {
  # Check that the seqnames of the gr and ref_genome match
  # .check_chroms(vcf_list, ref_genome)

  # Convert list to grl if necessary
  if (inherits(vcf_list, "list")) {
    vcf_list <- GenomicRanges::GRangesList(vcf_list)
  }

  if (inherits(vcf_list, "CompressedGRangesList")){

    # Unlist the GrangesList into a single GRanges object.
    gr <- BiocGenerics::unlist(vcf_list, use.names = FALSE)

    # Link samples to mutations
    nr_muts <- S4Vectors::elementNROWS(vcf_list)
    gr$INTERNAL_SAMPLENAME <- rep(names(nr_muts), times = nr_muts)

    # Get mutation context
    gr <- .get_indel_context_gr(gr, ref_genome)

    # Split the GRanges back into a list.
    sample_indx <- gr$INTERNAL_SAMPLENAME
    gr$INTERNAL_SAMPLENAME <- NULL
    grl <- S4Vectors::split(gr, sample_indx)
    return(grl)

  } else if (inherits(vcf_list, "GRanges")) {
    gr <- .get_indel_context_gr(vcf_list, ref_genome)
    return(gr)
  } else {
    .not_gr_or_grl(vcf_list)
  }
}

#' Get indel contexts from a single gr
#'
#' @details
#' Determines the COSMIC context from a GRanges object containing Indel mutations.
#' It throws an error if there are any variants with multiple alternative alleles or SNVs.
#'
#' @param gr GRanges object containing Indel mutations.
#' The mutations should be called similarly to HaplotypeCaller.
#' @param ref_genome BSgenome reference genome object
#'
#' @return A modified version of the input gr. In the gr two columns have been added.
#' "muttype" showing the main indel type and "muttype_sub" which shows the subtype.
#' The subtype is either the number of repeats or the microhomology length.
#'
#'
#' @noRd
#'
#' @importFrom magrittr %>%
#'
.get_indel_context_gr <- function(gr, ref_genome) {

  # Check that no snvs are present.
  # .check_no_substitutions(gr)

  # Calculate indel size to determine main category
  # ref_sizes <- gr %>%
  #   .get_ref() %>%
  #   width()
  # alt_sizes <- gr %>%
  #   .get_alt() %>%
  #   unlist() %>%
  #   width()
  ref_sizes <- unlist(lapply(gr$REF,length))
  alt_sizes <- unlist(lapply(unlist(gr$ALT),length))

  mut_size <- alt_sizes - ref_sizes

  # For the main indel categories, determine their sub categories.
  # (Also split the big deletion categorie into repeat and micro homology.)
  gr_1b_dels <- .get_1bp_dels(gr, mut_size, ref_genome)
  gr_1b_ins <- .get_1bp_ins(gr, mut_size, ref_genome)
  gr_big_dels <- .get_big_dels(gr, mut_size, ref_genome)
  gr_big_ins <- .get_big_ins(gr, mut_size, ref_genome)

  gr <- sort( c(gr_1b_dels, gr_1b_ins, gr_big_dels, gr_big_ins) )
  #   BiocGenerics::sort()
  # gr <- c(gr_1b_dels, gr_1b_ins, gr_big_dels, gr_big_ins) %>%
  #   BiocGenerics::sort()
  return(gr)
}

#' Get contexts from 1bp deletions
#'
#' @details
#' Determines the COSMIC context for the 1bp deletions in a GRanges object containing Indel mutations.
#' This function is called by get_indel_context_gr.
#'
#'
#' @param gr GRanges object containing Indel mutations.
#' The mutations should be called similarly to HaplotypeCaller.
#' @param mut_size A double vector containing the size of each Indel.
#' @param ref_genome BSgenome reference genome object
#'
#' @return A modified version of the input gr.
#' All variants that are not a 1bp deletion are removed.
#' In each gr two columns have been added.
#' "muttype" showing the main indel type and "muttype_sub" which shows the subtype.
#' The subtype is the number of repeats.
#'
#'
#' @importFrom magrittr %>%
#' @noRd
#'
#'
.get_1bp_dels <- function(gr, mut_size, ref_genome) {

  # Select mutations
  gr <- gr[mut_size == -1]
  if (length(gr) == 0) {
    return(gr)
  }

  # Get the deleted bases
  changed_bases <- unlist(lapply(as.character(gr$REF),function(x) { substring(x,2,2) }))
  # changed_bases <- gr %>%
  #   .get_ref() %>%
  #   as.character() %>%
  #   substring(2)

  # Remove any potential names from the ALT column so str_replace_all will work
  names(changed_bases) <- NULL

  # Get homopolymer length
  homopolymer_length <- .indel_get_n_repeats(gr, 19, ref_genome, changed_bases, "deletion")
  changed_bases[changed_bases == "A"] <- "T"
  changed_bases[changed_bases == "G"] <- "C"

  # Return the results
  gr$muttype <- stringr::str_c(changed_bases, "_deletion")
  gr$muttype_sub <- homopolymer_length
  return(gr)
}

#' Get contexts from 1bp insertions
#'
#' @details
#' Determines the COSMIC context for the 1bp insertions in a GRanges object
#' containing Indel mutations. This function is called by get_indel_context_gr.
#'
#'
#' @param gr GRanges object containing Indel mutations.
#' The mutations should be called similarly to HaplotypeCaller.
#' @param mut_size A double vector containing the size of each Indel.
#' @param ref_genome BSgenome reference genome object
#'
#' @return A modified version of the input gr.
#' All variants that are not a 1bp insertion are removed.
#' In each gr two columns have been added.
#' "muttype" showing the main indel type and "muttype_sub" which shows the subtype.
#' The subtype is the number of repeats.
#'
#'
#' @importFrom magrittr %>%
#' @noRd
#'
.get_1bp_ins <- function(gr, mut_size, ref_genome) {

  # Select mutations
  gr <- gr[mut_size == 1]
  if (length(gr) == 0) {
    return(gr)
  }

  # Get inserted bases.
  changed_bases <- unlist(lapply(as.character(unlist(gr$ALT)),function(x) { substring(x,2,2) }))
  # changed_bases <- gr %>%
  #   .get_alt() %>%
  #   unlist() %>%
  #   as.character() %>%
  #   substring(2)

  # Remove any potential names from the ALT column so str_replace_all will work
  names(changed_bases) <- NULL



  # Get homopolymer length
  homopolymer_length <- .indel_get_n_repeats(gr, 20, ref_genome, changed_bases, "insertion")

  changed_bases[changed_bases == "A"] <- "T"
  changed_bases[changed_bases == "G"] <- "C"

  # Return the results
  gr$muttype <- stringr::str_c(changed_bases, "_insertion")
  gr$muttype_sub <- homopolymer_length
  return(gr)
}

#' Get contexts from bigger inserions
#'
#' @details
#' Determines the COSMIC context for insertions larger than 1bp in a GRanges
#' object containing Indel mutations. This function is called by
#' get_indel_context_gr.
#'
#'
#' @param gr GRanges object containing Indel mutations.
#' The mutations should be called similarly to HaplotypeCaller.
#' @param mut_size A double vector containing the size of each Indel.
#' @param ref_genome BSgenome reference genome object
#'
#' @return A modified version of the input gr.
#' All variants that are not insertions larger than 1bp are removed.
#' In each gr two columns have been added.
#' "muttype" showing the main indel type and "muttype_sub" which shows the subtype.
#' The subtype is the number of repeats.
#'
#'
#' @importFrom magrittr %>%
#' @noRd
#'
.get_big_ins <- function(gr, mut_size, ref_genome) {

  # Select mutations
  gr <- gr[mut_size > 1]
  if (length(gr) == 0) {
    return(gr)
  }
  mut_size <- mut_size[mut_size > 1]

  # Get inserted bases
  changed_bases <- unlist(lapply(as.character(unlist(gr$ALT)),function(x) { substring(x,2,2) }))

  # changed_bases <- gr %>%
  #   .get_alt() %>%
  #   unlist() %>%
  #   as.character() %>%
  #   substring(2)
  biggest_ins <- max(nchar(changed_bases))
  # biggest_ins <- changed_bases %>%
  #   nchar() %>%
  #   max()
  flank_dist <- biggest_ins * 20

  # Remove any potential names from the ALT column so str_replace_all will work
  names(changed_bases) <- NULL

  # Get number of repeats
  n_repeats <- .indel_get_n_repeats(gr, flank_dist, ref_genome, changed_bases, "insertion")


  # Return results
  gr$muttype <- stringr::str_c(mut_size, "bp_insertion")
  gr$muttype_sub <- n_repeats
  return(gr)
}

#' Get contexts from larger deletions
#'
#' @details
#' Determines the COSMIC context for deletions larger than 1bp in a GRanges object
#' containing Indel mutations.
#' This function is called by get_indel_context_gr.
#' The function determines if there is microhomology for deletions that are not in repeat regions.
#'
#' @param gr GRanges object containing Indel mutations.
#' The mutations should be called similarly to HaplotypeCaller.
#' @param mut_size A double vector containing the size of each Indel.
#' @param ref_genome BSgenome reference genome object
#'
#' @return A modified version of the input gr.
#' All variants that are not a deletion larger than 1bp are removed.
#' In each gr two columns have been added.
#' "muttype" showing the main indel type and "muttype_sub" which shows the subtype.
#' The subtype is the number of repeats or the microhomology length.
#'
#'
#' @importFrom magrittr %>%
#' @noRd
#'
.get_big_dels <- function(gr, mut_size, ref_genome) {

  # Select mutations
  gr <- gr[mut_size < -1]
  if (length(gr) == 0) {
    return(gr)
  }
  mut_size <- mut_size[mut_size < -1]

  # Get deleted bases
  changed_bases <- unlist(lapply(as.character(gr$REF),function(x) { substring(x,2,2) }))

  # changed_bases <- gr %>%
  #   .get_ref() %>%
  #   as.character() %>%
  #   substring(2)

  biggest_dels <- max(nchar(changed_bases))

  # biggest_dels <- changed_bases %>%
  #   nchar() %>%
  #   max()
  flank_dist <- biggest_dels * 20

  # Remove any potential names from the ALT column so str_replace_all will work
  names(changed_bases) <- NULL

  # Get number of repeats
  n_repeats <- .indel_get_n_repeats(gr, flank_dist, ref_genome, changed_bases, "deletion")

  gr$muttype <- stringr::str_c(abs(mut_size), "bp_deletion")
  gr$muttype_sub <- n_repeats


  # Determine if there is microhomology for deletions that are not in repeat regions.
  # There is always at least 1 'repeat', because of the deleted bases themselves.
  pos_mh <- gr$muttype_sub == 1

  gr_repeat <- gr[!pos_mh]
  gr_mh <- gr[pos_mh]

  if (length(gr_mh) == 0) {
    return(gr_repeat)
  }

  mut_size_mh <- mut_size[pos_mh]
  del_bases_mh <- changed_bases[pos_mh]
  del_bases_s <- strsplit(del_bases_mh, "")

  seq <- .get_extended_sequence(gr_mh, biggest_dels, ref_genome, "right")
  seq_s <- strsplit(as.character(seq), "")


  # Also check for microhomology to the left of the deletion.
  # For this take the reverse sequence to the left
  # of the deletion and the reverse deleted bases.
  rev_del_bases <- Biostrings::reverse(del_bases_mh)
  rev_l_del_bases_s <- strsplit(rev_del_bases, "")

  seq_left <- .get_extended_sequence(gr_mh, biggest_dels, ref_genome, "left")
  rev_l_seq_s <- strsplit(as.character(seq_left), "")

  # For each mutation determine how many bases show hm
  nr_pos_mh <- length(del_bases_s)
  nr_mh <- vector("list", nr_pos_mh)
  for (i in seq_len(nr_pos_mh)) {
    del_bases_sample <- del_bases_s[[i]]
    seq_s_sample <- seq_s[[i]][seq_len(length(del_bases_sample))]
    same <- del_bases_sample == seq_s_sample

    # Determine how many bases are the same before the first difference.
    # na.rm is for when a sequence has been trimmed.
    # r_nr_mh_sample <- cumprod(same) %>%
      # sum(na.rm = TRUE)
    r_nr_mh_sample <- sum(na.omit(cumprod(same)))

    l_del_bases_sample <- rev_l_del_bases_s[[i]]
    l_seq_s_sample <- rev_l_seq_s[[i]][seq_len(length(l_del_bases_sample))]
    l_same <- l_del_bases_sample == l_seq_s_sample
    # l_nr_mh_sample <- cumprod(l_same) %>%
    #   sum(na.rm = TRUE)
    l_nr_mh_sample <- sum(na.omit(cumprod(l_same)))

    nr_mh_sample <- max(r_nr_mh_sample, l_nr_mh_sample)

    nr_mh[[i]] <- nr_mh_sample
  }
  nr_mh <- unlist(nr_mh)

  # Update gr when mh is indeed present
  mh_f <- nr_mh > 0
  gr_mh$muttype[mh_f] <- stringr::str_c(
    abs(mut_size_mh[mh_f]),
    "bp_deletion_with_microhomology"
  )
  gr_mh$muttype_sub[mh_f] <- nr_mh[mh_f]

  # Combine muts with and without mh
  # gr <- c(gr_mh, gr_repeat) %>%
  #   BiocGenerics::sort()
  gr <- sort( c(gr_mh, gr_repeat) )

  return(gr)
}


#' Count how often an indel is repeated on both the left
#' and the right side of the indel. Takes the side with
#' the maximum number of repeats.
#'
#' @param gr GRanges object containing Indel mutations.
#' The mutations should be called similarly to HaplotypeCaller.
#' @param flank_dist A numeric vector of length one containing
#' the number of flanking base pairs.
#' @param ref_genome BSgenome reference genome object
#' @param changed_bases String of the mutated bases.
#' @param type Indel type. Either deletion or insertion.
#'
#' @return double vector containing the number of repeats per indel.
#' @noRd
#'
.indel_get_n_repeats = function(gr, flank_dist, ref_genome, changed_bases, type){

  # Determine the number of repeats for both the left and the right flank.
  n_repeats_right <- .indel_get_n_repeats_right(gr, flank_dist, ref_genome, changed_bases, type)
  n_repeats_left <- .indel_get_n_repeats_left(gr, flank_dist, ref_genome, changed_bases, type)

  n_repeats <- pmax(n_repeats_right, n_repeats_left)
  return(n_repeats)
}

#' Count how often an indel is repeated to the right
#'
#' @param gr GRanges object containing Indel mutations.
#' The mutations should be called similarly to HaplotypeCaller.
#' @param flank_dist A numeric vector of length one containing
#' the number of flanking base pairs.
#' @param ref_genome BSgenome reference genome object
#' @param changed_bases String of the mutated bases.
#' @param type Indel type. Either deletion or insertion.
#'
#' @return double vector containing the number of repeats per indel.
#' @noRd
#'
.indel_get_n_repeats_right = function(gr, flank_dist, ref_genome, changed_bases, type){
  # Find extended sequence
  seq <- .get_extended_sequence(gr, flank_dist, ref_genome, direction = "right")

  n_repeats <- .indel_count_n_repeats(seq, changed_bases, type)
  return(n_repeats)
}

#' Count how often an indel is repeated to the left.
#'
#' @param gr GRanges object containing Indel mutations.
#' The mutations should be called similarly to HaplotypeCaller.
#' @param flank_dist A numeric vector of length one containing
#' the number of flanking base pairs.
#' @param ref_genome BSgenome reference genome object
#' @param changed_bases String of the mutated bases.
#' @param type Indel type. Either deletion or insertion.
#'
#' @return double vector containing the number of repeats per indel.
#' @noRd
#'
.indel_get_n_repeats_left = function(gr, flank_dist, ref_genome, changed_bases, type){
  # Find extended sequence
  seq <- .get_extended_sequence(gr, flank_dist, ref_genome, direction = "left")

  # Reverse the bases, to look in the left direction.
  changed_bases <- Biostrings::reverse(changed_bases)


  n_repeats <- .indel_count_n_repeats(seq, changed_bases, type)
  return(n_repeats)
}

#' Retreive the flanking sequence from a GRanges object
#'
#' @details
#' This function retreives the flanking sequence from a GRanges object.
#' The size of the flanking sequence is supplied as an argument to the function.
#' This function works by first extending the ranges and then retreiving the sequence.
#'
#'
#'
#' @param gr GRanges object containing Indel mutations.
#' The mutations should be called similarly to HaplotypeCaller.
#' @param flank_dist A numeric vector of length one containing the number of flanking base pairs.
#' @param ref_genome BSgenome reference genome object
#' @param direction Whether to flank to the right or left of the indel
#' @return A DNAStringSet containing the flanking bases.
#'
#' @noRd
#'
.get_extended_sequence <- function(gr, flank_dist, ref_genome, direction = c("right", "left")) {

  # Match argument
  direction <- match.arg(direction)


  if (direction == "right"){
    # Flank the granges object on the right, to get a sequence, that can be
    # searched for repeats. This can result in a warning message, when the
    # flanked range extends beyond the chrom length. This message is suppressed.
    withCallingHandlers(
      {
        gr_extended <- GenomicRanges::flank(gr, flank_dist, start = FALSE)
      },
      warning = function(w) {
        if (grepl("out-of-bound range located on sequence", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  } else{
    # Flank the granges object on the left, to get a sequence, that can be
    # searched for repeats. This can result in a warning message, when the
    # flanked range extends beyond the chrom length. This message is suppressed.
    withCallingHandlers(
      {
        gr_extended <- GenomicRanges::flank(gr, flank_dist-1, start = TRUE)
        BiocGenerics::end(gr_extended) <- BiocGenerics::end(gr_extended) + 1
      },
      warning = function(w) {
        if (grepl("out-of-bound range located on sequence", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  # Trim the ranges that are extended beyond the actual length of the chromosome.
  gr_extended <- GenomicRanges::trim(gr_extended)
  seq <- Biostrings::getSeq(BSgenome::getBSgenome(ref_genome), gr_extended)

  #Reverse the sequence if its retreived at the left.
  if (direction == "left"){
    seq <- Biostrings::reverse(seq)
  }
  return(seq)
}

#' Counts how often the mutated bases of an indel are repeated
#' in the surrounding context.
#'
#' @param seq A DNAStringSet containing the flanking bases.
#' @param changed_bases String of the mutated bases.
#' @param type Indel type. Either deletion or insertion.
#'
#' @return double vector containing the number of repeats per indel.
#' @noRd
#'
.indel_count_n_repeats = function(seq, changed_bases, type = c("deletion", "insertion")){

  # Match argument
  type <- match.arg(type)

  # Determine nr. repeats.
  # For each mut replace the deleted basetype in the flanking sequence with Zs.
  seq_z <- stringr::str_replace_all(
    as.character(seq),
    changed_bases,
    rep("Z", length(seq))
  )

  # Remove all bases after the Zs and count how many bases are left.
  # Add +1 for the deleted bases itself
  n_repeats <- nchar(gsub("[^Z].*", "", as.character(seq_z)))

  # n_repeats <- gsub("[^Z].*", "", as.character(seq_z)) %>%
    # nchar()

  if (type == "deletion"){
    n_repeats <- n_repeats + 1
  }
  return(n_repeats)
}

remove_mult_alts = function(gr) {

  mult_alts = elementNROWS(gr$ALT) > 1
  nr_mult_alts = sum(mult_alts)
  if (nr_mult_alts > 0){
    gr = gr[!mult_alts]
    warning(paste0("Removed ", nr_mult_alts, " variant(s), because they had multiple alternative alleles."), call. = F)
  }
  return(gr)
}

select_snvs = function(mutation_gr) {
  # select only snvs in the indel
  mutation_gr = remove_mult_alts(mutation_gr)
  mutation_gr_snv = mutation_gr[nchar(mutation_gr$REF) == 1 &
                                  nchar(unlist(mutation_gr$ALT)) == 1 &
                                  !is.na(seqnames(mutation_gr)) & start(ranges(mutation_gr)) != -1 & seqnames(mutation_gr) != "chrNA" &
                                  mutation_gr$REF != "-" &
                                  unlist(mutation_gr$ALT) != "-", ]

  strand(mutation_gr_snv) = ifelse(mutation_gr_snv$REF == "C" | mutation_gr_snv$REF == "T", "+", "-")

  return(mutation_gr_snv)
}

select_indels = function(mutation_gr) {
  # select only snvs in the indel
  mutation_gr = remove_mult_alts(mutation_gr)
  mutation_gr_indel = mutation_gr[(nchar(mutation_gr$REF) != 1 | nchar(unlist(mutation_gr$ALT)) != 1 ) &
                                  !is.na(seqnames(mutation_gr)) & start(ranges(mutation_gr)) != -1 & seqnames(mutation_gr) != "chrNA" &
                                  mutation_gr$REF != "-" &
                                  unlist(mutation_gr$ALT) != "-", ]

  # strand(mutation_gr_indel) = ifelse(mutation_gr_indel$REF == "C" | mutation_gr_indel$REF == "T" | unlist(mutation_gr_indel$ALT) == "C" | unlist(mutation_gr_indel$ALT) == "T", "+", "-")
  strand(mutation_gr_indel) <- "+"

  return(mutation_gr_indel)
}

select_context_snv = function(gr) {
  gr = remove_mult_alts(gr)
  if ( length(gr) == 0 ) {
    return( gr )
  }
  chr = seqnames(gr)
  start1 = as.integer( start(gr) - 10 )
  end1 = start(gr)-1
  start2 = end(gr)+1
  end2 = as.integer( end(gr) + 10 )
  strand = as.character( strand(gr) )
  context1 = getSeq(Hsapiens, names = chr, start =start1 , end = end1 ,strand = strand, as.character = T )
  context2 = getSeq(Hsapiens, names = chr, start =start2 , end = end2 ,strand = strand, as.character = T )

  context_tmp <- context1[which(strand == '-')]
  context1[which(strand == '-')] <- context2[which(strand == '-')]
  context2[which(strand == '-')] <- context_tmp

  refs = gr$REF
  alts = as.character(unlist(gr$ALT))
  refs[which(strand == '-' & refs == "G")] <- "C"
  refs[which(strand == '-' & refs == "A")] <- "T"
  alts2 <- alts
  alts2[which(strand == '-' & alts == "A")] <- "T"
  alts2[which(strand == '-' & alts == "C")] <- "G"
  alts2[which(strand == '-' & alts == "G")] <- "C"
  alts2[which(strand == '-' & alts == "T")] <- "A"
  alts <- alts2

  gr$context = paste(context1,"(",paste(refs,alts,sep=">"),")",context2,sep="")

  return(gr)
}

select_context_indel = function(gr) {
  gr = remove_mult_alts(gr)
  if ( length(gr) == 0 ) {
    return( gr )
  }
  chr = seqnames(gr)
  start1 = as.integer( start(gr) - 9 )
  end1 = start(gr)
  start2 = end(gr)+1
  end2 = as.integer( end(gr) + 10 )
  strand = as.character( strand(gr) )
  context1 = getSeq(Hsapiens, names = chr, start =start1 , end = end1 ,strand = strand, as.character = T )
  context2 = getSeq(Hsapiens, names = chr, start =start2 , end = end2 ,strand = strand, as.character = T )

  gr$context = paste(context1,"(",gr$indeltype,")",context2,sep="")

  return(gr)
}

rename_indel_types = function( gr ) {
  for( i in c(1:length(gr)) ) {
    muttype_split <- strsplit(gr[i]$muttype,"bp")[[1]]
    if ( length(muttype_split) > 1 ) {
      if ( as.numeric(muttype_split[1]) >= 5 ) {
        muttype_split[1] <- "5+"
        gr[i]$muttype <- paste(muttype_split[1],muttype_split[2],sep="bp")
      }
    }
  }
  ins_mh_rows <- which(grepl("insertion|microhomology", gr$muttype))
  del_rows <- c(1:length(gr))[-ins_mh_rows]
  if ( length(which(grepl("insertion|microhomology", gr$muttype) & as.numeric(gr$muttype_sub) >= 5)) > 0 ) {
    gr[which(grepl("insertion|microhomology", gr$muttype) & as.numeric(gr$muttype_sub) >= 5),]$muttype_sub <- "5+"
  }
  if ( length(which(grepl("deletion", gr$muttype) & as.numeric(gr$muttype_sub) >= 6)) ) {
    gr[which(grepl("deletion", gr$muttype) & as.numeric(gr$muttype_sub) >= 6),]$muttype_sub <- "6+"
  }
  # gr[which(grepl("insertion|microhomology", gr$muttype) & gr$muttype_sub >= 5),]$muttype_sub <- "5+"
  # gr[which(grepl("deletion", gr$muttype) & gr$muttype_sub >= 6),]$muttype_sub <- "6+"

  gr$indeltype <- paste(gr$muttype,gr$muttype_sub,sep="-")

  return(gr)
}

get_context = function(gr_in, output_file){
  gr_snvs = select_snvs(gr_in)
  gr_snvs_context = select_context_snv(gr_snvs)

  gr_indels = select_indels(gr_in)
  if ( length(gr_indels) > 0 ) {
    gr_indels_muttypes = get_indel_context( gr_indels, ref_genome )
    gr_indels_muttypes = rename_indel_types( gr_indels_muttypes )
    gr_indels_context = select_context_indel( gr_indels_muttypes )
  }

  if ( length(gr_snvs) > 0 ) {
    df_snvs <- data.frame(seqnames=gsub("chr","",seqnames(gr_snvs)),
                   starts=start(gr_snvs)-1,
                   ends=end(gr_snvs),
                   names=gr_snvs_context$context,
                   scores=c(rep(".", length(gr_snvs))),
                   strands=strand(gr_snvs)
    )
  } else {
    df_snvs <- data.frame()
  }
  if ( length(gr_indels) > 0 ) {
    df_indels <- data.frame(seqnames=gsub("chr","",seqnames(gr_indels)),
                     starts=start(gr_indels)-1,
                     ends=end(gr_indels),
                     names=gr_indels_context$context,
                     scores=c(rep(".", length(gr_indels))),
                     strands=strand(gr_indels)
    )
  } else {
    df_indels <- data.frame()
  }
  df <- rbind(df_snvs, df_indels)

  write.table(df, file=output_file, quote=F, sep="\t", row.names=F, col.names=F)
}

# get_context = function(grl, output_file, n_cores = 1){
#
#   grl_snvs = mclapply(grl, FUN = select_snvs, mc.cores = n_cores)
#   grl_snvs_context = mclapply(grl_snvs, FUN= select_context_snv, mc.cores = n_cores)
#
#   grl_indels = mclapply(grl, FUN = select_indels, mc.cores = n_cores)
#
#   grl_indels_muttypes = get_indel_context( grl_indels, ref_genome )
#   grl_indels_muttypes = mclapply(grl_indels_muttypes, FUN= rename_indel_types, mc.cores = n_cores)
#   grl_indels_context = mclapply(grl_indels_muttypes, FUN= select_context_indel, mc.cores = n_cores)
#
#   for (name in names(grl)) {
#     gr_snvs = grl_snvs_context[[name]]
#     if (!is.null(gr_snvs) ) {
#       df_snvs <- data.frame(seqnames=gsub("chr","",seqnames(gr_snvs)),
#                      starts=start(gr_snvs)-1,
#                      ends=end(gr_snvs),
#                      names=gr_snvs$context,
#                      scores=c(rep(".", length(gr_snvs))),
#                      strands=strand(gr_snvs)
#       )
#     } else {
#       df_snvs <- data.frame()
#     }
#     gr_indels = grl_indels_context[[name]]
#     if (!is.null(gr_indels) ) {
#       df_indels <- data.frame(seqnames=gsub("chr","",seqnames(gr_indels)),
#                        starts=start(gr_indels)-1,
#                        ends=end(gr_indels),
#                        names=gr_indels$context,
#                        scores=c(rep(".", length(gr_indels))),
#                        strands=strand(gr_indels)
#       )
#     } else {
#       df_indels <- data.frame()
#     }
#     df <- rbind(df_snvs, df_indels)
#     write.table(df, file=output_file, quote=F, sep="\t", row.names=F, col.names=F)
#   }
# }

vcf_name <- gsub(".+/(.+).vcf","\\1",vcf_file)

chroms <- c(1:22,"X","Y")

vcf = readVcf(vcf_file)
vcf <- vcf[which(seqnames(vcf) %in% chroms),]
gr = granges(vcf)

#seqlevelsStyle(gr) = "UCSC"
seqlevels(gr) = paste0("chr", seqlevels(gr))

get_context(gr, out_file)
