
fetch_range_sequences <- function(fasta, ranges) {
  if(!"peak_id" %in% names(ranges)) {
    ranges <- ranges |>
      dplyr::mutate(peak_id = stringr::str_glue("{chr}:{start}-{end}"))
  }

  if(!"strand" %in% names(ranges)) {
    ranges <- ranges |>
      dplyr::mutate(strand = "+")
  }

  fa <- if(inherits(fasta, "FaFile")) {
    fasta
  } else {
    Rsamtools::FaFile(fasta)
  }

  seqlevels_fa <- Rsamtools::scanFaIndex(fa) |>
    GenomeInfoDb::seqnames() |>
    as.character()

  d <- ranges |>
    dplyr::filter(chr %in% seqlevels_fa) |>
    dplyr::mutate(
      strand = as.character(strand),
      strand = dplyr::case_when(
        strand %in% c("-", "-1", "-1L") ~ "-",
        TRUE ~ "+"
      )
    )

  missing_chr <- setdiff(unique(d$chr), seqlevels_fa)
  if(length(missing_chr) > 0) {
    stop("Chromosomes not found in FASTA: ", paste(missing_chr, collapse = ", "))
  }

  gr <- GenomicRanges::makeGRangesFromDataFrame(d)
  seqs <- Rsamtools::scanFa(fa, param = gr)

  rev <- d$strand == "-"
  seqs[rev] <- Biostrings::reverseComplement(seqs[rev])

  d |>
    dplyr::mutate(sequence = as.character(seqs)) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(
      gc = Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence), "GC", as.prob = TRUE)[, 1]
    ) |> 
    dplyr::ungroup()
}


get_jaspar_matrix <- function() {
  jaspar <- JASPAR2024::JASPAR2024()
  sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), JASPAR2024::db(jaspar))
  opts <- list(
    collection = "CORE",
    tax_group = "vertebrates",
    species = 9606,
    all_versions = FALSE,
    matrixtype = "PWM"
  )
  mtx <- TFBSTools::getMatrixSet(sq24, opts)
  RSQLite::dbDisconnect(sq24)
  universalmotif::convert_motifs(mtx)
}

# seq_tb needs a unique peak_id and sequence
search_sequences_for_motifs <- function(seq_tb, motifs) {
  seqs <- Biostrings::DNAStringSet(seq_tb$sequence)
  names(seqs) <- paste0(seq_tb$peak_id)

  universalmotif::scan_sequences(
    motifs,
    seqs,
    RC = TRUE,
    threshold = 0.001,
    threshold.type = "pvalue",
    calc.pvals = TRUE,
    calc.qvals = TRUE,
    return.granges = TRUE,
    no.overlaps = TRUE
  ) |> 
    as_tibble() |> 
    select(
      peak_id = seqnames,
      motif_start = start,
      motif_end = end,
      motif_width = width,
      motif_strand = strand,
      motif_name = motif,
      match_seq = match,
      motif_score = score,
      motif_min_score = min.score,
      motif_max_score = max.score,
      motif_p_value = pvalue,
      motif_fdr = qvalue
    ) |> 
    arrange(motif_fdr)
}


remove_overlapping_peaks <- function(consensus_peaks, target_peaks, minoverlap = 1L, ignore.strand = TRUE) {
  stopifnot(all(c("chr", "start", "end") %in% names(consensus_peaks)))
  stopifnot(all(c("chr", "start", "end") %in% names(target_peaks)))

  universe_gr <- GenomicRanges::makeGRangesFromDataFrame(consensus_peaks)
  target_gr <- GenomicRanges::makeGRangesFromDataFrame(target_peaks)

  hits <- GenomicRanges::findOverlaps(
    query = universe_gr,
    subject = target_gr,
    minoverlap = minoverlap,
    ignore.strand = ignore.strand
  )

  consensus_peaks[-unique(S4Vectors::queryHits(hits)), , drop = FALSE]
}

make_matched_background <- function(target, universe, n_per_target = 1, seed = 42) {
  set.seed(seed)

  target_ids <- unique(target$peak_id)

  bg_pool <- universe |>
    filter(!peak_id %in% target_ids) |>
    mutate(
      width_bin = ntile(width, 10),
      gc_bin = ntile(gc, 10),
      bg_stratum = interaction(width_bin, gc_bin, drop = TRUE)
    )

  target_strata <- target |>
    mutate(
      width_bin = ntile(width, 10),
      gc_bin = ntile(gc, 10),
      bg_stratum = interaction(width_bin, gc_bin, drop = TRUE)
    ) |>
    count(bg_stratum, name = "n_target")

  bg_pool |>
    inner_join(target_strata, by = "bg_stratum") |>
    group_by(bg_stratum) |>
    group_split() |>
    purrr::map_dfr(function(x) {
      n_slice <- min(nrow(x), x$n_target[1] * n_per_target)
      dplyr::slice_sample(x, n = n_slice)
    })
}


db_motif_enrichment <- function(dec, consensus_peaks, genome_file, fdr_limit = 0.01, logfc_limit = 3) {
  target_peaks <- dec |> 
    filter(fdr < fdr_limit & abs(log_fc) > logfc_limit)
  universe_peaks <- remove_overlapping_peaks(consensus_peaks, target_peaks, minoverlap = 50)

  target_seq <- fetch_range_sequences(genome_file, target_peaks) |> 
    select(peak_id, chr, start, end, gc, sequence) |> 
    mutate(width = end - start + 1)
  universe_seq <- fetch_range_sequences(genome_file, universe_peaks) |> 
    select(peak_id, chr, start, end, gc, sequence) |> 
    mutate(width = end - start + 1)

}