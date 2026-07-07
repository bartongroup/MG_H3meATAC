read_one_peak_file <- function(fname, format = c("narrow", "broad")) {
  format <- match.arg(format)
  macs_colnames <- c(
    "chr",
    "start",
    "end",
    "peak_id",
    "score",
    "strand",
    "signal",
    "p_value",
    "q_value"
  )
  if(format == "narrow") {
    macs_colnames <- c(macs_colnames, "summit")
  }
  d <- read_tsv(fname, col_names = macs_colnames, show_col_types = FALSE) |> 
    mutate(
      across(c(start, end, score), as.integer),
      across(c(p_value, q_value), function(x) {10^-x})
    )
  if(format == "narrow") {
    d <- d |> 
      mutate(summit = as.integer(summit))
  }
  d
}

read_peak_files <- function(config, meta, format = c("narrow", "broad")) {
  format <- match.arg(format)
  ext <- ifelse(format == "narrow", "narrowPeak", "broadPeak")
  path <- file.path(config$data_path, config$dir[[str_glue("macs_{format}")]])

  pk <- map2(meta$raw_sample, meta$sample, function(rsmpl, smpl) {
    f <- file.path(path, str_glue("{rsmpl}_peaks.{ext}"))
    read_one_peak_file(f, format) |> 
      add_column(sample = smpl)
  }) |> 
    bind_rows()

  mt <- meta |> 
    select(sample, group)

  # Add group (condition)
  pk |> 
    left_join(mt, by = join_by(sample))
}


plot_peak_distribution_sample <- function(npeaks, smpl, chrom_sizes, bin_size = 1, bar_size = 1.2, y_scale = 0.8, chromosomes = c(1:22, "X", "Y")) {  
  cs <- chrom_sizes |> 
    filter(chr != "MT") |> 
    mutate(
      chr = factor(chr, levels = chromosomes) |> fct_rev(),
      y = chr |> as.integer(),
      chrom_size = chrom_size / 1e6
    )
  np <- npeaks |> 
    filter(sample == smpl) |>
    mutate(
      across(c(start, end), \(x) {x / 1e6}),
      mid = (start + end) / 2,
      chr = factor(chr, levels = chromosomes) |> fct_rev(),
      y = chr |> as.integer()
    ) |> 
    mutate(bin = cut(mid, breaks = seq(0, max(cs$chrom_size), by = bin_size)))  |> 
    group_by(bin, chr, y) |> 
    summarise(
      mid = mean(mid),
      n = n(),
      s = log(mean(signal))
    ) |> 
    ungroup() |> 
    mutate(
      s = s / max(s) * y_scale,
      n = n / max(n) * y_scale
    )
  
  ggplot() +
    th() +
    geom_segment(data = cs, aes(x = 0, xend = chrom_size, y = chr, yend = chr), colour = "grey80") +
    geom_segment(data = np, aes(x = mid, xend = mid, y = y, yend = y + n, group = y), linewidth = bar_size) +
    labs(x = "Position on chromosome (Mb)", y = "Normalized peak count")
}

# Normfac is a tibble with sample and norm columns
plot_peak_distribution_chrom <- function(npeaks, chrom, meta, normfac = NULL, bin_size = 1, bar_size = 1.2) {  
  np <- npeaks |> 
    filter(chr == chrom & str_detect(sample, "Me3")) |>
    mutate(
      across(c(start, end), \(x) {x / 1e6}),
      mid = (start + end) / 2,
    ) |> 
    mutate(bin = cut(mid, breaks = seq(0, max(mid), by = bin_size)))  |> 
    group_by(sample, bin) |> 
    summarise(
      mid = mean(mid),
      n = n(),
      s = log(mean(signal))
    ) |> 
    ungroup() |> 
    mutate(
      s = s / max(s) * 0.8,
      n = n / max(n) * 0.8
    ) |> 
    left_join(meta, by = join_by(sample)) |> 
    select(sample, mid, n, group)
  
  if(!is.null(normfac)) {
    np <- np |> 
      left_join(normfac, by = join_by(sample)) |> 
      mutate(n = n / norm) |> 
      select(-norm)
  }

  np |> 
    ggplot() +
    th() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    geom_col(aes(x = mid, y = n, fill = group), width = bin_size) +
    scale_fill_manual(values = okabe_ito_palette) +
    facet_wrap(~sample, ncol = 1, strip.position = "right") +
    labs(x = "Position on chromosome (Mb)", y = "Normalized peak count")
}


