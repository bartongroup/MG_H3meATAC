read_metadata <- function(meta_file) {
  read_csv(meta_file, show_col_types = FALSE) |> 
    mutate(replicate = as.character(replicate)) |> 
    separate_wider_delim(group, delim = "_", names = c("time_point", "treatment"), cols_remove = FALSE) |> 
    unite(sample, c(time_point, treatment, replicate), remove = FALSE) |> 
    mutate(sample = str_remove(sample, "DB1")) |> 
    mutate(
      time_point = fct_relevel(time_point, "iPS"),
      treatment = fct_relevel(treatment, "Ctrl"),
      replicate = fct_reorder(replicate, as.integer(replicate))
    ) |> 
    arrange(time_point, treatment, replicate) |> 
    mutate(sample = as_factor(sample), group = as_factor(group))
}


write_dec <- function(dec, name) {
  fname = file.path("tab", str_glue("{name}.csv"))
  dec |> 
    mutate(across(where(is.double), ~signif(.x, 4))) |> 
    write_csv(fname)
}

write_bed <- function(pks, fname, extra_cols = NULL) {
  pks |> 
    mutate(chr = paste0("chr", chr)) |>
    select(chr, start, end, all_of(extra_cols)) |> 
    write_tsv(fname, col_names = FALSE)
}


make_gr_regions <- function(df, chromosomes = c(1:22, "X", "Y")) {
  df |> 
    distinct(chr, start, end) |> 
    filter(chr %in% chromosomes) |> 
    mutate(chr = paste0("chr", chr)) |> 
    GenomicRanges::makeGRangesFromDataFrame()
}


random_overlap_test <- function(reg1, reg2, n_perm = 300) {
  pt <- regioneR::overlapPermTest(
    A = r1,
    B = r2,
    ntimes = n_perm,
    force.parallel = TRUE
  ) |> 
    suppressMessages() |> 
    suppressWarnings()
}

num_overlaps_once <- function(A, B) {
  regioneR::numOverlaps(A, B, count.once = TRUE)
}

rmsk_stats <- function(rmsk, pks, group_var = "rep_fam", n_perm = 100, n_top = 20) {
  # This is used implicitly by regioneR::randomizeRegions()
  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
    stop("Package BSgenome.Hsapiens.UCSC.hg38 is required for rmsk_stats().")
  }
  gr_pk <- make_gr_regions(pks)
  rmsk <- rmsk |> 
    mutate(grp = get(group_var))

  groups <- rmsk$grp |> 
    unique()

  rc <- map(groups, function(group) {
    gr_rm <- rmsk |> 
      filter(grp == group) |> 
      make_gr_regions()
    ov <- GenomicRanges::countOverlaps(gr_rm, gr_pk)
    n_over <- sum(ov > 0)
    tibble(
      group = group,
      n_rmask = length(gr_rm),
      n_over = n_over,
      perc = 100 * n_over / length(gr_rm),
      median_width = median(width(gr_rm))
    )
  }) |> 
    list_rbind() |> 
    arrange(desc(perc)) 

  perm <- map(1:n_top, function(i) {
    r <- rc[i, ]
    gr_rm <- rmsk |> 
      filter(grp == r$group) |> 
      make_gr_regions()

    print(paste("Permutation for", r$group))
    pt <- regioneR::permTest(
      A = gr_rm,
      B = gr_pk,
      randomize.function = function(A, ...) {
        regioneR::randomizeRegions(A = A, genome = "hg38")
      },
      evaluate.function = num_overlaps_once,
      alternative = "greater",
      ntimes = n_perm,
      force.parallel = TRUE
    ) |> 
      suppressMessages() |> 
      suppressWarnings()

    tibble(
      group = r$group,
      p_value = pt$num_overlaps_once$pval,
      perm_over = pt$num_overlaps_once$permuted |> mean()
    )
  }) |> 
    list_rbind()

  rc |> 
    left_join(perm, by = join_by(group))
}




filter_dec <- function(dec, ctrs, fdr_limit = 0.01, logfc_limit = 0.76, direction = c("up", "down")) {
  direction <- match.arg(direction)
  dec |> 
    mutate(dir = case_when(
      log_fc > 0 ~ "up",
      log_fc <= 0 ~ "down"
    )) |>
    filter(contrast %in% ctrs & fdr < fdr_limit & abs(log_fc) > logfc_limit & dir == direction)
}



make_random_regions <- function(peaks, chrom_sizes, n_regions = 100000,
                                chromosomes = c(as.character(1:22), "X", "Y")) {
  cs <- chrom_sizes |>
    dplyr::filter(chr %in% chromosomes) |>
    dplyr::select(chr, chrom_size)

  widths <- peaks |>
    dplyr::filter(chr %in% chromosomes) |>
    dplyr::slice_sample(n = n_regions) |> 
    dplyr::mutate(width = end - start) |>
    dplyr::pull(width)
    
  # For each sampled width, choose a chromosome with probability proportional
  # to the number of valid start positions available on that chromosome.
  chr_id <- purrr::map_int(widths, \(w) {
    weights <- pmax(cs$chrom_size - w + 1, 0)
    sample.int(nrow(cs), size = 1, prob = weights)
  })

  max_start <- cs$chrom_size[chr_id] - widths

  tibble::tibble(
    peak_id = sprintf("peak_%d", 1:length(widths)),
    chr = cs$chr[chr_id],
    width = widths,
    start = purrr::map_int(max_start, \(m) sample.int(m + 1, size = 1) - 1L),
    end = start + width
  ) |>
    dplyr::select(peak_id, chr, start, end, width)
}

plot_peak_coverage <- function(peaks, chrom_sizes, meta, annot = NULL) {
  if(!is.null(annot)) {
    peaks <- peaks |> 
      filter(annotation == annot)
  }

  tot <- sum(chrom_sizes$chrom_size)
  
  d <- peaks |> 
    distinct(sample, chr, start, end) |> 
    filter(chr %in% chrom_sizes$chr) |> 
    mutate(width = end - start + 1) |> 
    group_by(sample) |> 
    summarise(peak_coverage = sum(width), .groups = "drop") |> 
    mutate(peak_perc = 100 * peak_coverage / tot) |> 
    left_join(meta, by = "sample")

  env <- new.env(parent = globalenv())
  env$d <- d

  with(env, {
    d |> 
      ggplot(aes(x = sample, y = peak_perc, fill = group)) +
      th() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      ) +
      geom_col() +
      scale_fill_manual(values = okabe_ito_palette) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
      labs(x = NULL, y = "Genome coverage by peaks (%)")
  })
}

get_motif_stats <- function(dec_res, dec_for_motifs) {
  d1 <- dec_for_motifs |> 
    distinct(contrast_direction, peak_id) |> 
    count(contrast_direction, name = "DB peaks")
  d2 <- dec_res |> 
    distinct(contrast_direction, peak_id) |> 
    count(contrast_direction, name = "Peaks with hits")
  d3 <- dec_res |> 
    distinct(contrast_direction, motif_name, peak_id) |> 
    count(contrast_direction, peak_id) |> 
    group_by(contrast_direction) |> 
    summarise(`Motifs per peak` = mean(n), .groups = "drop")
  d4 <- dec_res |> 
    distinct(contrast_direction, motif_name) |> 
    count(contrast_direction, name = "Unique motifs")

  inner_join(d1, d2, by = "contrast_direction") |>
    inner_join(d3, by = "contrast_direction") |>
    inner_join(d4, by = "contrast_direction")
}

get_top_motif_stats <- function(dec_res, n_top = 30) {
  dec_res |> 
    slice_max(order_by = motif_score, n = n_top)
}


tabulate_annotations <- function(an, group_var, annot_var, n_top = NULL, rand_an = NULL, order_annot = FALSE) {
  d <- an |> 
    mutate(grp = get(group_var), ann = get(annot_var))

  # Some data don't have signal
  if(!("signal" %in% names(d))) {
    d <- d |> mutate(signal = 0)
  }

  d <- d |> 
    distinct(grp, peak_id, ann, signal)

  if(!is.null(n_top)) {
    d <- d |> 
      group_by(grp) |> 
      slice_max(signal, n = n_top) |>
      ungroup() |> 
      select(-signal)
  }

  if(!is.null(rand_an)) {
    ra <- rand_an |> 
      mutate(ann = get(annot_var)) |> 
      select(peak_id, ann) |> 
      mutate(grp = "Background")
    d <- d |> 
      bind_rows(ra)
  }

  dc <- d |> 
    count(grp, ann, .drop = FALSE) |> 
    group_split(grp) |> 
    map(function(w) {
      w |> 
        mutate(tot = sum(n)) |> 
        select(grp, ann, n, tot) |> 
        nest(data = c(n, tot)) |> 
        mutate(
          tst = map(data, ~prop.test(.x$n, .x$tot)),
          tidied = map(tst, tidy)
        ) |> 
        unnest(tidied) |> 
        unnest(data) |> 
        select(-tst)
    }) |> 
    list_rbind() |> 
    rename(prop = estimate) |> 
    mutate(grp = fct_relevel(grp, "Background", after = Inf))

  if(order_annot)
    dc <- dc |> mutate(ann = fct_reorder(ann, -n))
  
  if(annot_var == "rep_category") {
    dc <- dc |>
      mutate(ann = factor(ann, levels = c("SINE", "LINE", "LTR", "Simple_repeat", "DNA", "Low_complexity", "Satellite", "Retroposon", 
"Other")))
  }

  dc
}


# Fisher's test of each annotation and group against pooled others
# dc created by tabulate_annotations, must contain columns 'ann', 'grp' and 'n'
annot_group_test <- function(dc) {
  dc <- dc |>
    mutate(
      ann = as.character(ann),
      grp = as.character(grp)
    )

  crossing(
    ann = unique(dc$ann),
    grp = unique(dc$grp)
  ) |>
    mutate(
      tidied = map2(ann, grp, function(ant, gr) {
        d <- dc |>
          summarise(
            a = sum(n[ann == ant & grp == gr]),
            b = sum(n[ann != ant & grp == gr]),
            c = sum(n[ann == ant & grp != gr]),
            d = sum(n[ann != ant & grp != gr])
          )
        x <- d |> 
          as.numeric() |> 
          matrix(nrow = 2, byrow = TRUE)

        depleted <- x[1, 1] / sum(x[1, ]) < x[2, 1] / sum(x[2, ])
        n_obs <- x[1, 1]

        fisher.test(x, alternative = "two.sided") |>
          broom::tidy() |>
          mutate(
            #estimate = if_else(depleted, -estimate, estimate),
            estimate = estimate,
            n_ann_grp = d$a,
            n_notann_grp = d$b,
            n_ann_notgrp = d$c,
            n_notann_notgrp = d$d
          )
      })
    ) |>
    unnest(tidied) |>
    select(
      ann, grp,
      starts_with("n_"),
      odds_ratio = estimate,
      ci_lo = conf.low,
      ci_up = conf.high,
      p_value = p.value
    ) |>
    mutate(
      fdr = p.adjust(p_value, method = "BH"),
      log_odds_ratio = if_else(odds_ratio == 0, 0, log2(odds_ratio)),
      p_in = n_ann_grp / (n_ann_grp + n_notann_grp),
      p_out = n_ann_notgrp / (n_ann_notgrp + n_notann_notgrp)
    )
}


annotation_fisher_example <- function(tst, annot = "Promoter", group = "48h_Ctrl-iPS_Ctrl_up") {
  d <- tst |> 
    filter(ann == annot & grp == group) |> 
    select(starts_with("n_"))

  tribble(
    ~., ~`With annotation`, ~`Without annotation`,
    "In group", d$n_ann_grp, d$n_notann_grp,
    "Not in group", d$n_ann_notgrp, d$n_notann_notgrp
  )
}


add_unannotated_peaks <- function(peaks, annotated, annot_var,
                                  none = "none",
                                  key = c("peak_id", "contrast")) {
  stopifnot(
    all(key %in% names(peaks)),
    all(key %in% names(annotated)),
    annot_var %in% names(annotated)
  )

  # The unannotated DE table should contain one row per peak and contrast.
  if(anyDuplicated(peaks[key])) {
    stop("`peaks` contains duplicated peak/contrast keys.")
  }

  annotated <- annotated |>
    mutate("{annot_var}" := coalesce(as.character(.data[[annot_var]]), none))

  annotated_keys <- annotated |>
    select(all_of(key)) |>
    distinct()

  unannotated <- peaks |>
    anti_join(annotated_keys, by = key) |>
    mutate("{annot_var}" := none)

  bind_rows(annotated, unannotated)
}


annot_de_test <- function(dec, fdr_limit, logfc_limit, annot_var) {
  dc <- dec |> 
    mutate(
      ann = get(annot_var),
      sig = fdr < fdr_limit & abs(log_fc) > logfc_limit
    ) |> 
    select(peak_id, chr, start, end, width, strand, contrast, contrast_direction, sig, ann)

  # "none" indicates a peak missing an annotation
  anns <- dc |>
    distinct(ann) |>
    filter(!is.na(ann), ann != "none") |>
    pull(ann) |>
    as.character()

  dct <- dc |> distinct(contrast, contrast_direction)

  map2(dct$contrast, dct$contrast_direction, function(ctr, ctrd) {
    d_sig <- dc |>
      filter(contrast_direction == ctrd & sig)
    d_nsig <- dc |>
      filter(contrast == ctr & !sig)

    if(n_distinct(d_sig$peak_id) < 10)
      return(NULL)

    # peaks have multiple annotations, now look at peaks only
    # we ignore GC content here, hence fake gc = 1 for make_matched_background()
    pks_sig <- d_sig |> distinct(peak_id, chr, start, end, width, strand) |> add_column(gc = 1)
    pks_nsig <- d_nsig |> distinct(peak_id, chr, start, end, width, strand) |> add_column(gc = 1)
    bkg <- make_matched_background(pks_sig, pks_nsig, gc_bins = 1)

    d_bkg <- d_nsig |> 
      filter(peak_id %in% bkg$background$peak_id)

    n_sig <- n_distinct(d_sig$peak_id)
    n_bkg <- n_distinct(d_bkg$peak_id)

    map(anns, function(an) {
      a <- d_sig |>
        filter(ann == an) |>
        distinct(peak_id) |>
        nrow()
      
      c <- d_bkg |>
        filter(ann == an) |>
        distinct(peak_id) |>
        nrow()
      
      d <- tibble(
        a = a,
        b = n_sig - a,
        c = c,
        d = n_bkg - c
      )

      x <- d |> unlist() |> matrix(nrow = 2, byrow = TRUE)

      n_obs <- x[1, 1]

      fisher.test(x, alternative = "two.sided") |>
        broom::tidy() |>
        mutate(
          estimate = estimate,
          n_ann_grp = d$a,
          n_notann_grp = d$b,
          n_ann_notgrp = d$c,
          n_notann_notgrp = d$d,
          ann = an,
          grp = ctrd
        )
    }) |>
    list_rbind() |> 
    select(
      ann,
      grp,
      starts_with("n_"),
      odds_ratio = estimate,
      ci_lo = conf.low,
      ci_up = conf.high,
      p_value = p.value
    ) |>
    mutate(
      fdr = p.adjust(p_value, method = "BH"),
      log_odds_ratio = log2(odds_ratio),
      p_in = n_ann_grp / (n_ann_grp + n_notann_grp),
      p_out = n_ann_notgrp / (n_ann_notgrp + n_notann_notgrp)
    )
  }) |> 
    list_rbind()
}


# boths decs contain only significant peaks
compare_cut_atac_ <- function(dec_atac, dec_cut) {
  
  gr_atac <- GenomicRanges::makeGRangesFromDataFrame(dec_atac)
  tot_atac <- dec_atac |> 
    count(contrast_direction, name = "n_atac")
  
  map(c("Ctrl", "SETDB1"), function(treat) {
    map(c("down", "up"), function(dir) {
      ctr <- str_glue("48h_{treat}-iPS_{treat}_{dir}")
      gr_cut <- dec_cut |> 
        filter(contrast_direction == ctr) |> 
        GenomicRanges::makeGRangesFromDataFrame()
      hits <- GenomicRanges::findOverlaps(
        query = gr_cut,
        subject = gr_atac,
        ignore.strand = TRUE
      )
      cut_hits <- S4Vectors::queryHits(hits) |> unique()
      atac_hits <- S4Vectors::subjectHits(hits) |> unique()

      dec_atac[atac_hits, ] |> 
        count(contrast_direction, name = "n_cut_atac") |> 
        add_column(cut_contrast = ctr) |> 
        left_join(tot_atac, by = "contrast_direction") |> 
        mutate(
        )
    }) |> 
      list_rbind()
  }) |> 
    list_rbind() |> 
    rename(atac_contrast = contrast_direction) |> 
    relocate(cut_contrast, .before = 1)
}




compare_cut_atac <- function(dec_atac, dec_cut) {
  
  use_contrasts <- c(
    "48h_Ctrl-iPS_Ctrl_down",
    "48h_Ctrl-iPS_Ctrl_up",
    "48h_SETDB1-iPS_SETDB1_down",
    "48h_SETDB1-iPS_SETDB1_up"
  )
  
  map(use_contrasts, function(cut_ctr) {
    cut <- dec_cut |> 
      filter(contrast_direction == cut_ctr & sig) |> 
      distinct(chr, start, end)
    gr_cut <- GenomicRanges::makeGRangesFromDataFrame(cut)
    
    map(use_contrasts, function(atac_ctr) {
      atac <- dec_atac |> 
        filter(contrast_direction == atac_ctr) |> 
        distinct(chr, start, end, sig)
      gr_atac <- GenomicRanges::makeGRangesFromDataFrame(atac)
      cut_overlap <- GenomicRanges::countOverlaps(
        gr_atac,
        gr_cut,
        ignore.strand = TRUE
      ) > 0
      
      d <- tibble(
        n_cut_atac = sum(atac$sig & cut_overlap),
        n_atac_ncut = sum(atac$sig & !cut_overlap),
        n_ncut_atac = sum(!atac$sig & cut_overlap),
        n_natac_ncut = sum(!atac$sig & !cut_overlap)
      )
      
      x <- d |> 
        as.numeric() |> 
        matrix(nrow = 2, byrow = TRUE)
        cutover_prop = x[1, 1] / sum(x[, 1])
        ncutover_prop = x[1, 2] / sum(x[, 2])
      
      fisher.test(x, alternative = "greater") |> 
      broom::tidy() |> 
      mutate(
        cut_contrast = cut_ctr,
        atac_contrast = atac_ctr,
        n_cut = nrow(cut),
        n_atac = sum(atac$sig),
        cutover_prop,
        ncutover_prop,
        atac_status = case_when(
          cutover_prop > ncutover_prop ~ "enriched",
          cutover_prop < ncutover_prop ~ "depleted",
          TRUE ~ "unchanged"
        ),
        x = list(x)
      )
    }) |> 
    list_rbind()
  }) |>
  list_rbind() |> 
    select(
      cut_contrast, atac_contrast,
      starts_with("n_"),
      ends_with("_prop"),
      atac_status,
      x,
      odds_ratio = estimate,
      ci_lo = conf.low,
      ci_up = conf.high,
      p_value = p.value
    ) |> 
  mutate(fdr = p.adjust(p_value, method = "BH"))
}

get_peak_overlaps <- function(pks1, pks2, suffix1 = "_atac", suffix2 = "_cut") {

  pks1 <- pks1 |> mutate(idx = row_number())
  pks2 <- pks2 |> mutate(idx = row_number())
  over <- GenomicRanges::findOverlaps(
    query = GenomicRanges::makeGRangesFromDataFrame(pks1),
    subject = GenomicRanges::makeGRangesFromDataFrame(pks2),
    ignore.strand = TRUE
  )

  over |> 
    as_tibble() |> 
    left_join(pks1, by = join_by(queryHits == idx)) |> 
    left_join(pks2, by = join_by(subjectHits == idx), suffix = c(suffix1, suffix2))
}


get_cut_atac_overlaps <- function(dec_atac, dec_cut, ctr_atac, ctr_cut) {
  cut <- dec_cut |> 
    filter(contrast == ctr_cut)
  atac <- dec_atac |> 
    filter(contrast == ctr_atac)
 
  get_peak_overlaps(atac, cut)
}

count_cut_atac_overlaps <- function(dec_atac, dec_cut) {
  cut <- dec_cut |> 
    distinct(chr, start, end, strand)
  atac <- dec_atac |> 
    distinct(chr, start, end, strand)
  over <- GenomicRanges::findOverlaps(
    query = GenomicRanges::makeGRangesFromDataFrame(atac),
    subject = GenomicRanges::makeGRangesFromDataFrame(cut),
    ignore.strand = TRUE
  )

  tibble(
    n_cut = nrow(cut),
    n_atac = nrow(atac),
    n_over = length(over)
  )
}

# Do the promoters that gain accessibility upon SETDB1 depletion
# (in contrast 48h_SETDB1-48h_Ctrl_up) loose accessibility upon SETDB1
# depletion in CUT&Tag?
question_1 <- function(dec_atac, dec_cut) {
  d_atac <- dec_atac |> 
    filter(contrast_direction == "48h_SETDB1-48h_Ctrl_up" & sig & annotation == "Promoter")
  d_cut <- dec_cut |> 
    filter(contrast == "48h_SETDB1-48h_Ctrl")
  ov <- get_peak_overlaps(d_atac, d_cut) |> 
    select(peak_id_atac, log_fc_atac, fdr_atac, sig_atac, log_fc_cut, fdr_cut, sig_cut, gene_id, gene_symbol, gene_name) |>
    arrange(fdr_cut)
}

# Do LTR rmasker regions that gain accessibility upon SETDB1 depletion
# (in contrast 48h_SETDB1-48h_Ctrl_up) loose accessibility upon SETDB1
# depletion in CUT&Tag?
question_2 <- function(dec_rmsk_groups_libs_narrow, dec_cut) {
  d_atac <- dec_rmsk_groups_libs_narrow |> 
    mutate(sig = fdr < 0.01 & abs(log_fc) > 1) |> 
    filter(contrast_direction == "48h_SETDB1-48h_Ctrl_up" & sig & rep_category == "LTR")
  d_cut <- dec_cut |> 
    filter(contrast == "48h_SETDB1-48h_Ctrl")
  ov <- get_peak_overlaps(d_atac, d_cut) |> 
    select(peak_id_atac, log_fc_atac, fdr_atac, sig_atac, log_fc_cut, fdr_cut, sig_cut) |>
    arrange(fdr_cut)
}

# Look at motifs that are enriched in 48h_Ctrl-iPS_Ctrl_up.
# Is there a change in CUT&Tag level for these motifs in 48h_Ctrl-iPS_Ctrl?
# In other words, does the opening of a motif target regions coincide with loss of H3K9me3?
question_3 <- function(dec_atac, tf_enr, tf_scan, dec_cut) {
  scn <- tf_scan |> 
    filter(motif_fdr < 0.01) |> 
    distinct(motif_id, peak_id) |> 
    separate_wider_delim(peak_id, delim = ":", names = c("chr", "startend"), cols_remove = FALSE)  |> 
    separate_wider_delim(startend, delim = "-", names = c("start", "end"), cols_remove = TRUE) |> 
    mutate(start = as.integer(start), end = as.integer(end)) 
  atac_sig <- dec_atac |> 
    filter(contrast_direction == "48h_Ctrl-iPS_Ctrl_up" & sig) |> 
    distinct(peak_id)
  enr <- tf_enr |>
    filter(contrast_direction == "48h_Ctrl-iPS_Ctrl_up") |> 
    select(motif_id, motif_name, motif_consensus, target_enrichment) |> 
    inner_join(scn, by = join_by(motif_id)) |> 
    group_by(motif_consensus, chr, start, end, peak_id) |> 
    summarise(motif_name = str_c(sort(motif_name), collapse = ";"), target_enrichment = mean(target_enrichment), .groups = "drop") |> 
    filter(peak_id %in% atac_sig$peak_id)
  d_cut <- dec_cut |> 
    filter(contrast == "48h_Ctrl-iPS_Ctrl")
  ov <- get_peak_overlaps(enr, d_cut) |> 
     select(motif_name, motif_consensus, target_enrichment, peak_id_atac, peak_id_cut, log_fc_cut = log_fc, p_value_cut = p_value, fdr_cut = fdr, sig_cut = sig) |> 
     arrange(fdr_cut)
}

# Monte-Carlo to assess logFC distribution of CUT&Tag for motifs that are enriched in 48h_Ctrl-iPS_Ctrl_up.
mc_ov_cut <- function(dec_cut, ov) {
  bg <- dec_cut |>
    filter(contrast == "48h_Ctrl-iPS_Ctrl") |>
    pull(log_fc)

  obs <- ov |>
    group_by(motif_name) |>
    summarise(
      n = n(),
      obs_mean = mean(log_fc_cut),
      .groups = "drop"
    )

  set.seed(42)

  boot_means <- function(n, B = 10000) {
    replicate(B, mean(sample(bg, size = n, replace = FALSE)))
  }

  res <- obs |>
    mutate(
      null = map(n, boot_means),
      null_mean = map_dbl(null, mean),
      p_greater = map2_dbl(null, obs_mean, \(x, m) mean(x >= m)),
      p_less = map2_dbl(null, obs_mean, \(x, m) mean(x <= m)),
      p_two_sided = pmin(1, 2 * pmin(p_greater, p_less))
    ) |>
    mutate(fdr = p.adjust(p_two_sided, method = "BH"))
}

plot_question_3 <- function(ov, dec_cut) {
  M <- median(dec_cut$log_fc)
  md <- ov |> 
    group_by(motif_name) |> 
    summarise(m = mean(log_fc_cut), .groups = "drop") |> 
    mutate(motif_name = fct_reorder(motif_name, m))
  ov |>
    mutate(motif_name = factor(motif_name, levels = levels(md$motif_name))) |>
    ggplot(aes(x = log_fc_cut, y = motif_name, colour = sig_cut)) +
    th() +
    theme(legend.position = "none") +
    geom_quasirandom() +
    #stat_summary(fun = mean, geom = "crossbar", orientation = "y", linewidth = 3, colour = "red") +
    stat_summary(fun = mean, geom = "point", shape = 124, size = 6, colour = "red") +
    geom_vline(xintercept = 0, colour = "brown", alpha = 0.3) +
    scale_colour_manual(values = c("grey70", "black")) +
    labs(x = "log2FC in CUT&Tag (48h_Ctrl-iPS_Ctrl)", y = "Motif name")
}

# Identify genes whose promoters show changes in accessibility in ATAC&seq.
question_4 <- function(dec_chipseeker_groups_libs_narrow) {
  dec_chipseeker_groups_libs_narrow |> 
    filter(fdr < 0.01 & abs(log_fc) > 1 & annotation == "Promoter") |>
    distinct(contrast, log_fc, gene_symbol, gene_name, gene_biotype) |> 
    arrange(contrast, log_fc)
}

# RNA-seq data from https://www.nature.com/articles/s41586-020-2144-9
read_matsuda_data <- function() {
  path <- "info/41586_2020_2144_MOESM3_ESM.xlsx"
  stopifnot(file.exists(path))

  fill_right <- function(x) {
    for (i in seq_along(x)[-1]) {
      if (is.na(x[i])) x[i] <- x[i - 1]
    }
    x
  }

  simplify_condition <- function(x) {
    sub("^D[0-9]+_([^_]+)_([^＿]+)＿FPKM$", "\\1_\\2", x)
  }

  hdr <- read_excel(path, sheet = "FPKM", col_names = FALSE, n_max = 2, .name_repair = "minimal")
  condition <- as.character(hdr[1, ])
  replicate  <- as.character(hdr[2, ])

  condition <- fill_right(condition)
  simple_condition <- simplify_condition(condition)

  col_names <- c(
    condition[1:3],
    paste(simple_condition[-(1:3)], replicate[-(1:3)], sep = "_")
  )

  readxl::read_excel(
    path,
    sheet = "FPKM",
    skip = 2,
    col_names = col_names,
    .name_repair = "minimal"
  ) |> 
    rename(group = Group) |> 
    mutate(across(-c(gene_id, gene, group), as.numeric)) |> 
    pivot_longer(-c(gene_id, gene, group)) |> 
    separate_wider_delim(name, names = c("condition", "donor", "replicate"), delim = "_")
}


de_matsuda <- function(matsuda, dnr) {
  d <- matsuda |>
    filter(condition %in% c("iPS", "PS", "PSM") & donor == dnr) |> 
    unite(sample, c(condition, replicate), remove = FALSE)
  meta <- d |> 
    distinct(sample, condition, replicate) |> 
    mutate(condition = as_factor(condition))

  tab <- d |> 
    mutate(value = log2(value + 1)) |> 
    pivot_wider(id_cols = gene_id, names_from = sample, values_from = value) |> 
    column_to_rownames("gene_id") |> 
    as.matrix()
  tab <- tab[, meta$sample]

  design <- model.matrix(~ 0 + condition, data = meta)
  colnames(design) <- levels(meta$condition)

  fit <- limma::lmFit(tab, design)

  contrast <- limma::makeContrasts(
    PS_PSM_vs_iPS = (PS + PSM) / 2 - iPS,
    levels = design
  )

  fit2 <- limma::contrasts.fit(fit, contrast)
  fit2 <- limma::eBayes(fit2)

  limma::topTable(fit2, coef = "PS_PSM_vs_iPS", number = Inf) |> 
    as_tibble(rownames = "gene_id")
}


save_matsuda <- function(matsuda) {
  i2g <- matsuda |> distinct(gene_id, gene)

  matsuda |> 
    de_matsuda("201B7") |>
    left_join(i2g, by = join_by(gene_id)) |>
    relocate(gene, .after = gene_id) |> 
    write_csv("tab/matsuda_de_201B7.csv")

  matsuda |> 
    de_matsuda("1231A3") |>
    left_join(i2g, by = join_by(gene_id)) |>
    relocate(gene, .after = gene_id) |> 
    write_csv("tab/matsuda_de_1231A3.csv")

}