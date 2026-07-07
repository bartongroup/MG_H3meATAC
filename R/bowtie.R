read_one_bowtie_summary <- function(sfile) {
  if(!file.exists(sfile)) {
    message(str_glue("File {sfile} does not exist."))
    stop()
  }
  r <- read_lines(sfile) |> 
    str_trim()
  nums <- r[1:5] |> 
    str_remove("\\s.+$") |> 
    as.numeric()
  tibble(
    total_reads = nums[1],
    paired = nums[2],
    not_aligned = nums[3],
    aligned_once = nums[4],
    aligned_multi = nums[5]
  ) |> 
    mutate(align_rate = (aligned_once + aligned_multi) / total_reads)
}


read_bowtie_summaries <- function(config, meta, bam_dir = NULL) {
  if(is.null(bam_dir)) {
    bam_dir <- config$dir$bam
  }
  map2(meta$raw_sample, meta$sample, function(rsmpl, smpl) {
    f <- file.path(config$data_path, bam_dir, str_glue("{rsmpl}.summary.txt"))
    read_one_bowtie_summary(f) |> 
      add_column(sample = smpl)
  }) |> 
    bind_rows()
}

merge_two_bowtie_summaries <- function(sum_human, sum_dm) {
  bind_rows(
    sum_human |> add_column(species = "human"),
    sum_dm |> add_column(species = "drosophila")
  ) |> 
    rename(once = aligned_once, multi = aligned_multi) |> 
    mutate(
      once_perc = 100 * once / total_reads,
      multi_perc = 100 * multi / total_reads
    ) |> 
    pivot_wider(id_cols = c(sample, total_reads, paired), names_from = species, values_from = c(once, multi, once_perc, multi_perc))
}

plot_bowtie_summaries <- function(bsum, meta) {
  bsum |> 
    mutate(
      once_perc = 100 * aligned_once / total_reads,
      multi_perc = 100 * aligned_multi / total_reads,
      total_reads_m = total_reads / 1e6
    ) |> 
    rename(
      `Total reads (M)` = total_reads_m,
      `% aligned once` = once_perc,
      `% aligned multi` = multi_perc
    ) |> 
    select(sample, `Total reads (M)`, `% aligned once`, `% aligned multi`) |> 
    pivot_longer(-sample)  |> 
    left_join(meta, by = join_by(sample)) |> 
  ggplot(aes(x = group, y = value, fill = replicate)) +
    th() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_beeswarm(size = 2, cex = 3, shape = 21, colour = "grey50") +
    scale_fill_manual(values = okabe_ito_palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) +
    facet_wrap(~name, scales = "free_y", ncol = 1)    
}

read_one_picard_duplication <- function(dfile) {
  r <- read_lines(dfile)
  i0 <- which(str_detect(r, "^LIBRARY"))
  names <- r[i0] |> str_split_1("\\t")
  values <- r[i0 + 1] |> str_split_1("\\t")
  tibble(name = names, value = values)
}


read_picard_duplication <- function(config, meta) {
  map2(meta$raw_sample, meta$sample, function(rsmpl, smpl) {
    f <- file.path(config$data_path, config$dir$bam, str_glue("{rsmpl}_dedup.metrics.txt"))
    read_one_picard_duplication(f) |> 
      add_column(sample = smpl)
  }) |> 
    bind_rows()
}


plot_picard_duplication <- function(pd, meta) {
  pd |> 
    filter(name %in% c("READ_PAIRS_EXAMINED", "ESTIMATED_LIBRARY_SIZE", "PERCENT_DUPLICATION", "READ_PAIR_OPTICAL_DUPLICATES")) |> 
    mutate(value = as.numeric(value)) |> 
    pivot_wider(id_cols = sample) |> 
    mutate(
      `% optical duplicates` = 100 * READ_PAIR_OPTICAL_DUPLICATES / READ_PAIRS_EXAMINED,
      PERCENT_DUPLICATION = 100 * PERCENT_DUPLICATION,
      ESTIMATED_LIBRARY_SIZE = ESTIMATED_LIBRARY_SIZE / 1e6
    ) |> 
    rename(
      `% duplication` = PERCENT_DUPLICATION,
      `Estimated lib size (M)` = ESTIMATED_LIBRARY_SIZE
    ) |> 
    select(-c(READ_PAIRS_EXAMINED, READ_PAIR_OPTICAL_DUPLICATES)) |> 
    pivot_longer(-sample) |> 
    left_join(meta, by = join_by(sample)) |> 
  ggplot(aes(x = group, y = value, fill = replicate)) +
    th() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_beeswarm(size = 2, cex = 3, shape = 21, colour = "grey50") +
    scale_fill_manual(values = okabe_ito_palette) +   
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) +
    facet_wrap(~name, scales = "free_y", ncol = 1)    
}


plot_spikein_ratios <- function(bowtie_sum_merged, meta, brks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)) {
  bowtie_sum_merged |> 
    mutate(
      all_dp = once_drosophila + multi_drosophila,
      all_human = once_human + multi_human,
      dp_prop = 100 * all_dp / (all_dp + all_human)
    ) |> 
    left_join(meta, by = join_by(sample)) |> 
  ggplot(aes(x = group, y = dp_prop, fill = replicate)) +
    th() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_beeswarm(size = 2, cex = 3, shape = 21, colour = "grey50") +
    scale_fill_manual(values = okabe_ito_palette) +   
    scale_y_log10(breaks = brks, labels = as.character(brks)) +
    #scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA), trans = "sqrt", breaks = c(1, 2, 5, 10, 20, 40, 60)) +
    labs(x = "Group", y = "% drosophila hits")
}