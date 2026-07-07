plot_region_ggbio <- function(samples, config, meta, region, peaks_broad = NULL, peaks_narrow = NULL,
                               normfac = NULL, rel_scale = 5) {
  env <- new.env(parent = globalenv())
  
  chr <- GenomeInfoDb::seqlevels(region) |> as.character()

  bam_dat <- map(samples, function(smpl) {
    rsmpl <- meta |> 
      filter(sample == smpl) |> 
      pull(raw_sample) |> 
      as.character()
    bam_file <- file.path(config$data_path, config$dir$bam, str_glue("{rsmpl}.noMT.bam"))

    an <- tibble(text = smpl, x = -Inf, y = Inf, hj = -0.1, vj = 2)
    
    tr <- ggbio::autoplot(bam_file, which = region, stat = "coverage") 
    
    dat <- ggplot_build(tr@ggplot)$data[[1]] |> 
      as_tibble() |> 
      select(x = x, y = y) |> 
      add_column(sample = smpl) |> 
      left_join(meta, by = join_by(sample))

    # Scale data
    if(!is.null(normfac)) {
      scale_fac <- normfac |> 
        filter(sample == smpl) |> 
        pull(norm)

      dat <- dat |> 
        mutate(y = y / scale_fac)
    }
    dat
  }) |> 
    list_rbind() |> 
    mutate(sample = factor(sample, levels = levels(meta$sample)))

  y_max <- max(bam_dat$y)
  y_min <- -y_max / rel_scale
  y1 <- y_min / 3
  y2 <- 2 * y1

  if(is.null(peaks_broad)) peaks_broad <- tibble()
  if(is.null(peaks_narrow)) peaks_narrow <- tibble()
  
  peak_dat <- bind_rows(
    peaks_broad |> add_column(type = "broad", y = y1),
    peaks_narrow |> add_column(type = "narrow", y = y2)
  ) |>
    filter(sample %in% samples & chr == GenomeInfoDb::seqlevels(region) |> as.character() & start < end(region) & end > start(region))

  env$chr <- chr
  env$y_min <- y_min
  env$bam_dat <- bam_dat
  env$peak_dat <- peak_dat

  with(env, {
    bam_dat |> 
      ggplot() +
      th() +
      geom_area(aes(x = x / 1e6, y = y, fill = group)) +
      geom_hline(yintercept = 0, alpha = 0.2) +
      facet_wrap(~ sample, ncol = 1, strip.position = "right") +
      scale_fill_manual(values = okabe_ito_palette) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(y_min, NA)) +
      geom_segment(data = peak_dat, aes(x = start / 1e6, xend = end / 1e6, y = y, yend = y, colour = type), linewidth = 1.2) +
      scale_colour_brewer(palette = "Paired") +
      labs(x = str_glue("Position on chr {chr} (Mb)"), y = "Normalised count", fill = "Group", colour = "MACS")
  })
  
}


plot_region_ggbio_peak_id <- function(dec, peakid, config, meta, peaks_broad, peaks_narrow, normfac = NULL, f = 2) {
  if(is.numeric(peakid)) {
    peakid <- dec$peak_id[peakid]
  }

  pk <- dec |> 
    filter(peak_id == peakid)
  stopifnot(nrow(pk) == 1)

  #sample <- pk$sample[1] |> as.character()
  #samples <- meta |> filter(as.character(replicate) == "1") |> pull(sample) |> as.character()

  samples <- meta$sample

  region <- make_region(pk$chr, pk$start, pk$end) |> 
    expand_region(f = f)

  plot_region_ggbio(samples, config, meta, region, peaks_broad, peaks_narrow, normfac)
}


# This sometimes works, sometimes doesn't.
# ggbio objects are nasty
save_ggbio_png <- function(gg, file) {
  png(file, width = 1000, height = 700)
  gg
  dev.off()
}