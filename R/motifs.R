#' Make a named DNAStringSet from peak sequences
#'
#' Converts a table of peak sequences into a `DNAStringSet` and names each
#' sequence with its peak identifier.
#'
#' @param tb A data frame or tibble with `sequence` and `peak_id` columns.
#'
#' @return A `Biostrings::DNAStringSet` named by `peak_id`.
make_dna_string_set <- function(tb) {
  s <- Biostrings::DNAStringSet(tb$sequence)
  names(s) <- paste0(tb$peak_id)
  s
}


#' Fetch genomic sequences for ranges
#'
#' Extracts sequences from a FASTA file for genomic ranges, reverse-complements
#' sequences on the negative strand, and adds sequence width and GC content.
#'
#' @param fasta Path to an indexed FASTA file or an `Rsamtools::FaFile` object.
#' @param ranges A data frame or tibble with `chr`, `start`, and `end` columns.
#'   Optional `peak_id` values are used as sequence identifiers; if absent, they
#'   are created from the genomic coordinates. Optional `strand` values are used
#'   to orient sequences; missing or non-negative values are treated as `"+"`.
#'
#' @return A tibble containing the input ranges with `peak_id`, `strand`,
#'   `sequence`, `width`, and `gc` columns.
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

  missing_chr <- setdiff(unique(ranges$chr), seqlevels_fa)
  if(length(missing_chr) > 0) {
    stop("Chromosomes not found in FASTA: ", paste(missing_chr, collapse = ", "))
  }

  d <- ranges |>
    dplyr::filter(chr %in% seqlevels_fa) |>
    dplyr::mutate(
      strand = as.character(strand),
      strand = dplyr::case_when(
        strand %in% c("-", "-1", "-1L") ~ "-",
        TRUE ~ "+"
      )
    )

  gr <- GenomicRanges::makeGRangesFromDataFrame(d)
  seqs <- Rsamtools::scanFa(fa, param = gr)

  rev <- d$strand == "-"
  seqs[rev] <- Biostrings::reverseComplement(seqs[rev])

  d |>
    dplyr::mutate(
      sequence = as.character(seqs),
      width = end - start + 1
    ) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(
      gc = Biostrings::letterFrequency(Biostrings::DNAStringSet(sequence), "GC", as.prob = TRUE)[, 1]
    ) |> 
    dplyr::ungroup()
}


#' Get human vertebrate JASPAR CORE PWMs
#'
#' Loads the JASPAR2024 database, selects the latest human vertebrate CORE
#' position weight matrices, and converts them to `universalmotif` objects.
#'
#' @return A list of `universalmotif` motif objects.
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

#' Search target peaks for motif matches
#'
#' Extracts target peak sequences and scans them for motif matches using
#' `universalmotif::scan_sequences()`. Target peaks should usually be filtered
#' for statistically significant peaks before calling this function.
#'
#' @param target_peaks A data frame or tibble with `peak_id`, `chr`, `start`,
#'   `end`, and `strand` columns.
#' @param motifs A motif list compatible with `universalmotif::scan_sequences()`.
#' @param genome_file Path to an indexed FASTA file or an `Rsamtools::FaFile`
#'   object.
#'
#' @return A tibble of motif hits, sorted by motif FDR, with peak identifiers,
#'   motif identifiers, match coordinates, strand, sequence, score, p-value, and
#'   FDR columns.
search_peaks_for_motifs <- function(target_peaks, motifs, genome_file) {
  
  targets <- target_peaks |>
    dplyr::distinct(peak_id, chr, start, end, strand)
  target_seqs <- fetch_range_sequences(genome_file, targets)
  seqs <- make_dna_string_set(target_seqs)

  universalmotif::scan_sequences(
    motifs,
    seqs,
    RC = TRUE,
    threshold = 0.001,
    threshold.type = "pvalue",
    calc.pvals = TRUE,
    calc.qvals = TRUE,
    return.granges = TRUE,
    no.overlaps = TRUE,
    verbose = 2
  ) |> 
    as_tibble() |> 
    select(
      peak_id = seqnames,
      motif_id = motif.i,
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


#' Remove peaks overlapping target peaks
#'
#' Removes rows from a consensus peak set that overlap any target peak.
#'
#' @param consensus_peaks A data frame or tibble with `chr`, `start`, and `end`
#'   columns.
#' @param target_peaks A data frame or tibble with `chr`, `start`, and `end`
#'   columns.
#' @param minoverlap Minimum number of overlapping bases required to count as an
#'   overlap. Default is 1.
#' @param ignore.strand Logical; whether to ignore strand when finding overlaps.
#'   Default is `TRUE`.
#'
#' @return `consensus_peaks` with overlapping rows removed.
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

#' Make a GC- and width-matched background peak set
#'
#' Samples background peaks from a universe after stratifying target and universe
#' peaks by width and GC-content bins.
#'
#' @param targets A data frame or tibble with `peak_id`, `chr`, `start`, `end`,
#'   `strand`, `width`, and `gc` columns.
#' @param universe A data frame or tibble with `peak_id`, `chr`, `start`, `end`,
#'   `strand`, `width`, and `gc` columns.
#' @param seed Random seed used for reproducible background sampling. Default is
#'   42.
#' @param width_bins Number of quantile bins to use for peak width. Default is
#'   10.
#' @param gc_bins Number of quantile bins to use for GC content. Default is 10.
#' @param n_per_target Number of background peaks to sample per target peak
#'   within each stratum. Default is 1.
#'
#' @return A list with `background`, the sampled background peaks, and
#'   `diagnostics`, a tibble reporting available and sampled background counts
#'   for each stratum.
make_matched_background <- function(targets, universe, seed = 42,
                                    width_bins = 10, gc_bins = 10,
                                    n_per_target = 1) {
  set.seed(seed)

  stopifnot(all(c("peak_id", "width", "gc") %in% names(targets)))
  stopifnot(all(c("peak_id", "width", "gc") %in% names(universe)))

  target_ids <- unique(targets$peak_id)

  targets2 <- targets |>
    dplyr::ungroup() |>
    dplyr::select(peak_id, chr, start, end, strand, width, gc) |>
    dplyr::mutate(.set = "targets")

  universe2 <- universe |>
    dplyr::ungroup() |>
    dplyr::filter(!peak_id %in% target_ids) |>
    dplyr::select(peak_id, chr, start, end, strand, width, gc) |>
    dplyr::mutate(.set = "background")

  binned <- dplyr::bind_rows(targets2, universe2) |>
    dplyr::mutate(
      width_bin = dplyr::ntile(width, width_bins),
      gc_bin = dplyr::ntile(gc, gc_bins),
      bg_stratum = interaction(width_bin, gc_bin, drop = TRUE)
    )

  targets_strata <- binned |>
    dplyr::filter(.set == "targets") |>
    dplyr::count(bg_stratum, name = "n_targets")

  bg_pool <- binned |>
    dplyr::filter(.set == "background") |>
    dplyr::inner_join(targets_strata, by = "bg_stratum")

  bg_counts <- binned |>
    dplyr::filter(.set == "background") |>
    dplyr::count(bg_stratum, name = "n_background")

  diagnostics <- targets_strata |>
    dplyr::left_join(bg_counts, by = "bg_stratum") |>
    dplyr::mutate(
      n_background = dplyr::coalesce(n_background, 0L),
      n_requested = n_targets * n_per_target,
      n_sampled = pmin(n_background, n_requested),
      enough_background = n_background >= n_requested
    )

  background <- bg_pool |>
    dplyr::group_by(bg_stratum) |>
    dplyr::group_split() |>
    purrr::map_dfr(function(x) {
      n_slice <- min(nrow(x), x$n_targets[1] * n_per_target)
      dplyr::slice_sample(x, n = n_slice)
    }) |>
    dplyr::select(-.set, -n_targets)

  list(
    background = background,
    diagnostics = diagnostics
  )
}


#' Test motif enrichment in target peaks
#'
#' Runs motif enrichment for each contrast direction by comparing target peak
#' sequences with a GC- and width-matched background. Target peaks should usually
#' be filtered for statistically significant peaks before calling this function.
#'
#' @param target_peaks A data frame or tibble with `peak_id`, `chr`, `start`,
#'   `end`, `strand`, and `contrast_direction` columns.
#' @param consensus_seqs A data frame or tibble containing candidate background
#'   peaks with `peak_id`, `chr`, `start`, `end`, `strand`, `width`, `gc`, and
#'   `sequence` columns.
#' @param genome_file Path to an indexed FASTA file or an `Rsamtools::FaFile`
#'   object.
#' @param motifs A motif list compatible with `universalmotif::enrich_motifs()`.
#' @param min_target Minimum number of target peaks required to run enrichment
#'   for a contrast. Default is 50.
#' @param hit_qval_limit Motif-hit q-value threshold passed to
#'   `universalmotif::enrich_motifs()`. Default is 0.01.
#' @param enr_qval_limit Maximum enrichment q-value returned by
#'   `universalmotif::enrich_motifs()`. Default is 0.01.
#'
#' @return A tibble of motif enrichment results for all contrasts that pass
#'   `min_target`, with motif identifiers and a `contrast_direction` column.
motif_enrichment <- function(target_peaks, consensus_seqs, genome_file, motifs, min_target = 50,
                             hit_qval_limit = 0.01, enr_qval_limit = 0.01) {

  all_targets <- fetch_range_sequences(genome_file, target_peaks)

  contrasts <- unique(as.character(target_peaks$contrast_direction))
  map(contrasts, function(ctr) {
    message("Enrichment for contrast ", ctr)
    targets <- all_targets |> 
      filter(contrast_direction == ctr)
    if(nrow(targets) < min_target) {
      message("Not enough targets. Skipping contrast: ", ctr)
      return(NULL)
    }

    universe <- remove_overlapping_peaks(consensus_seqs, targets, minoverlap = 1)
    matched <- make_matched_background(targets, universe)
    not_enough <- matched$diagnostics |> filter(!enough_background)
    if(nrow(not_enough) > 0) {
      warning("Not enough background peaks for contrast: ", ctr)
      print(not_enough)
    }

    background <- matched$background |> 
      select(peak_id) |> 
      left_join(universe, by = join_by(peak_id)) # sequence is dropped in matching

    n_target <- nrow(targets)
    n_background <- nrow(background)

    if(n_background < n_target * 0.8) {
      warning(
        "Matched background for ", ctr, " contains only ",
        n_background, " peaks for ", n_target, " targets."
      )
    }

    target_seqs <- make_dna_string_set(targets)
    background_seqs <- make_dna_string_set(background)

    universalmotif::enrich_motifs(
      motifs = motifs,
      sequences = target_seqs,
      bkg.sequences = background_seqs,
      threshold = hit_qval_limit,
      threshold.type = "qvalue",
      max.q = enr_qval_limit,
      mode = "seq.hits",
      RC = TRUE,
      verbose = 2
    ) |> 
      as_tibble() |> 
      janitor::clean_names() |> 
      rename(
        motif_name = motif,
        motif_id = motif_i
      ) |> 
      add_column(contrast_direction = ctr)
  }) |> 
    list_rbind()
}
