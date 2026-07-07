okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#AAAAAA", "#000000")

annot_palette <- c("#F0E442", "#33A02C", "#A6CEE3", "#1F78B4", "#B2DF8A", "#f0dfd4", "#D55E00")

th <- function() {
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid = ggplot2::element_blank())
}

scientific_10 <- function(v, limit = 3.99) {
  sc <- map_chr(v, function(x) {
    if (is.na(x)) {
      as.character(x)
    } else if (x == 0 | abs(log10(x)) < limit) {
      format(x, scientific = FALSE)
    } else {
      format(x, scientific = TRUE) |> 
        str_replace("e", "%*%10^") |> 
        str_replace("\\+", "") |> 
        str_replace("\\^0+", "^")
    }
  })
  parse(text = sc)
}


gs <- function(gg, name, width, height) {
  ggplot2::ggsave(filename = file.path("doc/images", paste0(name, ".png")), plot = gg, device = "png",
         width = width, height = height, dpi = 300, limitsize = FALSE)
}

gp <- function(gg, name, width, height) {
  ggplot2::ggsave(filename = file.path("fig", paste0(name, ".pdf")), plot = gg, device = "pdf",
                  width = width, height = height)
}



plot_volma <- function(res, point_size, point_alpha, scales = "fixed") {
  r <- res |>
    select(x, y, sel, group)
  r_sel <- r |> filter(sel)
  r_nsel <- r |> filter(!sel)

  env <- new.env(parent = globalenv())
  env$point_size <- point_size
  env$point_alpha <- point_alpha
  env$scales <- scales
  env$r_sel <- r_sel
  env$r_nsel <- r_nsel

  with(env, {
    g <- ggplot() +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        legend.position = "none"
      ) +
      geom_point(data = r_nsel, aes(x = x, y = y), colour = "grey70",
                 size = point_size, alpha = point_alpha) +
      facet_wrap( ~ group, scale = scales)

    if (nrow(r_sel) > 0) {
      g <- g + geom_point(data = r_sel, aes(x = x, y = y), colour = "black",
                          size = point_size, alpha = point_alpha)
    }
    g
  })


}

plot_ma <- function(res, expr_var = log_cpm, fc_var = log_fc, p_var = p_value, fdr_var = fdr, 
                    group_var = contrast, sel_genes = NULL, fdr_limit = 0.01, logfc_limit = 0,
                    point_size = 0.5, point_alpha = 0.5) {
  
   if(is.null(sel_genes)) {
    res <- res |> mutate(sel = {{ fdr_var }} < fdr_limit & abs({{ fc_var }}) >= logfc_limit)
  } else {
    res <- res |> mutate(sel = id %in% sel_genes)
  }
  
  res |>
    mutate(
      group = {{ group_var }},
      x = {{ expr_var }},
      y = {{ fc_var }},
    ) |>
    plot_volma(point_size, point_alpha) +
    geom_hline(yintercept = 0, linewidth = 0.1, alpha = 0.5) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) +
    labs(x = expression(log[10]~Intensity), y = expression(log[2]~FC))
}

plot_volcano <- function(res, fc_var = log_fc, p_var = p_value, fdr_var = fdr, group_var = contrast,
                         sel_genes = NULL, fdr_limit = 0.01, logfc_limit = 1, point_size = 0.5,
                         point_alpha = 0.5, scales = "fixed") {
  if(is.null(sel_genes)) {
    res <- res |> mutate(sel = {{ fdr_var }} < fdr_limit & abs({{ fc_var }}) >= logfc_limit)
  } else {
    res <- res |> mutate(sel = id %in% sel_genes)
  }

  res |>
    mutate(
      group = {{ group_var }},
      x = {{ fc_var }},
      y = -log10({{ p_var}})
    ) |> 
    plot_volma(point_size, point_alpha, scales) +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA))
}

plot_pdist <- function(res, p_var = p_value, group_var = contrast, n_bins = 50) {
  brks <- seq(0, 1, length.out = n_bins)
  r <- res |>
    mutate(p = {{ p_var }}, grp = {{ group_var }}) |>
    select(p, grp)
  rm(res)
  ggplot(r, aes(x = p, y = after_stat(density))) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_histogram(breaks = brks) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    facet_wrap(~grp)
}

plot_vmp <- function(res, expr_var = log_cpm, fc_var = log_fc, p_var = p_value, fdr_var = fdr,
                     group_var = contrast, sel_genes = NULL, fdr_limit = 0.01, logfc_limit = 0, point_size = 0.5, point_alpha = 0.5, n_bins = 50){
  list(
    volcano = plot_volcano(res, {{fc_var}}, {{p_var}}, {{fdr_var}}, {{group_var}}, sel_genes, fdr_limit, logfc_limit, point_size, point_alpha, scales = "free_x"),
    ma = plot_ma(res, {{expr_var}}, {{fc_var}}, {{p_var}}, {{fdr_var}}, {{group_var}}, sel_genes, fdr_limit, logfc_limit, point_size, point_alpha),
    pdist = plot_pdist(res, {{p_var}}, {{group_var}}, n_bins)
  )
}


plot_up_down <- function(dec, fdr_limit = 0.01, logfc_limit = 0) {
  env <- new.env(parent = globalenv())

  env$d <- dec |> 
    filter(fdr < fdr_limit & abs(log_fc) > logfc_limit) |>
    mutate(direction = case_when(
      log_fc > 0 ~ "Up",
      log_fc <= 0 ~ "Down"
    )) |> 
    count(contrast, direction) |> 
    mutate(n = if_else(direction == "Down", -n, n)) 

  with(env, {
    d |> 
      ggplot(aes(x = contrast, y = n, fill = direction)) +
      th() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_col() +
      geom_hline(yintercept = 0, linewidth = 0.1, alpha = 0.5) +
      scale_fill_manual(values = c("Up" = "#0072B2", "Down" = "#D55E00")) +
      #scale_y_continuous(labels = scientific_10) +
      labs(x = NULL, y = "Number of differentially bound peaks")
  })
}

plot_de_peak_distribution <- function(dec, chrom_sizes, fdr_limit = 0.01, chromosomes = c(1:22, "X", "Y")) {
  cs <- chrom_sizes |> 
    filter(chr != "MT") |> 
    mutate(
      chr = factor(chr, levels = chromosomes) |> fct_rev(),
      y = chr |> as.integer(),
      chrom_size = chrom_size / 1e6
    ) 
  np <- dec |> 
    filter(fdr < fdr_limit) |>
    mutate(
      across(c(start, end), \(x) {x / 1e6}),
      mid = (start + end) / 2,
      chr = factor(chr, levels = chromosomes) |> fct_rev(),
      y = chr |> as.integer(),
      log_fc = log_fc / max(abs(log_fc)),
      sgn = sign(log_fc) |> as_factor()
    ) 
  
  ggplot() +
    th() +
    theme(legend.position = "none") +
    geom_segment(data = cs, aes(x = 0, xend = chrom_size, y = chr, yend = chr), colour = "grey80") +
    geom_segment(data = np, aes(x = mid, xend = mid, y = y, yend = y + log_fc, group = y, colour = sgn), linewidth = 0.5) +
    scale_colour_manual(values = c("1" = "blue", "-1" = "red", "0" = "grey50")) +
    facet_wrap(~contrast, ncol = 1, strip.position = "right") +
    labs(x = "Position on chromosome (Mb)", y = "Normalized peak count")
}


plot_number_of_peaks <- function(peaks, meta) {
  env <- new.env(parent = globalenv())

  env$d <- peaks |> 
    count(sample) |> 
    left_join(meta, by = "sample")

  with(env, {
    d |> 
      ggplot(aes(x = sample, y = n, fill = group)) +
      th() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      geom_col() +
      scale_fill_manual(values = okabe_ito_palette) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
      labs(x = NULL, y = "Number of peaks")
  })
}

plot_group_peak_tss <- function(an, group_var = "group", n_top = 30000, rand_an = NULL, max_dist = 300) {
  d <- an |> 
    mutate(
      grp = get(group_var)
    ) |> 
    distinct(grp, peak_id, distance_to_tss, signal)

  if(!is.null(n_top)) {
    d <- d |> 
      group_by(grp) |> 
      slice_max(signal, n = n_top) |>
      ungroup() |> 
      select(-signal)
  }

  if(!is.null(rand_an)) {
    ra <- rand_an |> 
      select(peak_id, distance_to_tss) |> 
      mutate(grp = "Random")
    d <- d |> 
      bind_rows(ra)
  }

  d <- d |> 
    mutate(tsskb = distance_to_tss / 1000) |> 
    filter(abs(tsskb) <= max_dist)
  
  d |> 
    ggplot(aes(x = grp, y = distance_to_tss / 1000, colour = grp)) +
    th() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    theme(legend.position = "none") +
    geom_quasirandom(size = 0.3) +
    geom_hline(yintercept = 0, linewidth = 0.1, alpha = 0.5) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = "Distance to TSS (kb)")
}


plot_group_peak_annotations <- function(an, group_var = "group", annot_var = "annotation", rand_an = NULL,
                                        n_top = 30000, order_annot = FALSE, palette = 1, stacked = TRUE,
                                        plot_proportions = TRUE) {
  
  dc <- tabulate_annotations(an, group_var, annot_var, n_top, rand_an, order_annot)

  if(palette == 1) {
    pal <- annot_palette
  } else if(palette == 2) {
    pal <- c("#FDBF6F", annot_palette)
  } else if(palette == 3) {
    pal <- RColorBrewer::brewer.pal(12, "Paired")
  }

  if(stacked) {
    pos <- "stack"
    top_margin <- 0
  } else {
    pos <- position_dodge(width = 0.9)
    top_margin <- 0.03
  }

   if(plot_proportions) {
     dc <- dc |> mutate(y = prop)
     ylab <- "Proportion"
   } else {
     dc <- dc |> mutate(y = n)
     top_margin <- 0.03
     ylab <- "Count"
   }

  env <- new.env(parent = globalenv())
  env$dc <- dc
  env$pal <- pal
  env$pos <- pos
  env$stacked <- stacked
  env$top_margin <- top_margin
  env$ylab <- ylab

  with(env, {
    pl <- dc |> 
      ggplot(aes(x = grp, y = y, fill = ann)) +
      th() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      geom_col(position = pos) +
      #scale_fill_manual(values = okabe_ito_palette) +
      scale_fill_manual(values = pal) +
      scale_y_continuous(expand = expansion(mult = c(0, top_margin))) +
      labs(x = NULL, y = ylab, fill = "Annotation")
    if(!stacked) {
      pl <- pl + geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = pos, width = 0.2)
    }
    pl
  })

}

plot_peak_annotatios <- function(an, ncol = 4, group_var = "sample", annot_var = "annotation", scale_y = TRUE) {
  env <- new.env(parent = globalenv())
  env$ncol <- ncol
  env$group_var <- group_var
  env$annot_var = annot_var
  env$scale_y = scale_y

  if(scale_y) {
    sc <- 1000
    env$y_lab <- "Count (thousands)"
  } else {
    sc <- 1
    env$y_lab <- "Count"
  }

  env$d <- an |> 
    mutate(grp = get(group_var), ann = get(annot_var)) |> 
    distinct(grp, peak_id, ann) |> 
    count(grp, ann) |> 
    mutate(n = n / sc)
  with(env, {
    d |> 
      ggplot(aes(x = ann, y = n, fill = ann)) +
      th() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      geom_col() +
      facet_wrap(~grp, ncol = ncol, scales = "free_y") +
      scale_fill_brewer(palette = "Paired") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03)), labels = scientific_10) +
      labs(x = NULL, y = y_lab)
  })
}


plot_peak_rmsk <- function(dec, fdr_limit = NULL, ncol = 1, top_n = 16) {
  env <- new.env(parent = globalenv())
  env$ncol <- ncol

 
  if(!is.null(fdr_limit)) {
    # FDR selection per contrast
    dc <- dec |> 
      filter(fdr < fdr_limit)
  } else {
    # All peaks, fake contrast for display
    ctrs <- dec$contrast |> unique()
    dc <- dec |> 
      filter(contrast == ctrs[1]) |> 
      mutate(contrast = "All peaks")
  }  

  env$dc <- dc |> 
    distinct(contrast, peak_id, rep_fam) |> 
    count(contrast, rep_fam) |> 
    arrange(desc(n)) |> 
    mutate(f = fct_reorder(rep_fam, -n)) |> 
    group_split(contrast)  |> 
    map(function(w) {
      w |> 
        mutate(tot = sum(n)) |> 
        select(contrast, f, n, tot) |> 
        nest(data = c(n, tot)) |> 
        mutate(
          tst = map(data, ~prop.test(.x$n, .x$tot)),
          tidied = map(tst, tidy)
        ) |> 
        unnest(tidied) |> 
        unnest(data) |> 
        select(-tst) |> 
        head(top_n)
    }) |> 
    list_rbind()

  if(nrow(dc) == 0)
    return(NULL)

  with(env, {
    dc |> 
      ggplot(aes(x = f, y = estimate, fill = contrast)) +
      th() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top"
      ) +
      geom_col(width = 0.8, position = position_dodge(width = 0.8)) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.8), width = 0.2) +
      scale_fill_manual(values = okabe_ito_palette) +
      #facet_wrap(~contrast, ncol = ncol, scales = "fixed") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03)), labels = scales::comma) +
      labs(x = NULL, y = "Proportion of counts")
  })
}



plot_db_clustering <- function(db, meta, text_size = 10, colour_var = "group") {
  tab <- db$binding[, 4:ncol(db$binding)]

  dendr <- t(tab) |> 
    dist(method = "euclidean") |> 
    hclust(method = "complete") |>
    dendsort::dendsort() |> 
    ggdendro::dendro_data()
  seg <- ggdendro::segment(dendr)
  meta <- meta |>
    mutate(
      colvar = get(colour_var)
    )
  labs <- left_join(dendr$labels |> mutate(label = as.character(label)), meta, by = c("label" = "sample")) |> 
    mutate(colour = okabe_ito_palette[as_factor(colvar)])
 
 theme.d <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = text_size, colour = labs$colour),
    axis.line.y = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(linewidth = 0.5),
    axis.ticks.y = ggplot2::element_blank()
  )
  
  ggplot() +
    theme.d +
    coord_flip() +
    geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_continuous(breaks = seq_along(labs$label), labels = labs$label) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(seg$y) * 1.03)) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = "Distance")
}


plot_db_biding_region <- function(db, meta, region) {
  env <- new.env(parent = globalenv())
  
  chr <- seqnames(region) |> as.character()
  env$chr <- chr
  env$d <- db$binding |> 
    as_tibble() |>
    mutate(CHR = db$chrmap[CHR]) |> 
    filter(CHR == chr & END >= start(region) & START <= end(region)) |> 
    pivot_longer(-c(CHR, START, END), names_to = "sample", values_to = "value") |> 
    left_join(meta, by = "sample") |> 
    mutate(pos = (START + END) / 2 / 1e6)

  with(env, {
    ggplot(d, aes(x = pos, xend = pos, y = 0, yend = value, colour = group)) + 
      th() +
      theme(
        legend.position = "none",
        strip.text = element_text(size = 7) 
      ) +
      geom_segment() +
      scale_colour_manual(values = okabe_ito_palette) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
      facet_wrap(~sample, ncol = 4) +
      labs(x = str_glue("Position on chromosome {chr} (Mb)"), y = "Binding signal")
  })
}


plot_compare_spikein_normalisation <- function(normfac, spikein_norms) {
  normfac |>
    left_join(spikein_norms, by = join_by(sample)) |> 
    mutate(
      depth_norm = 1 / scale,
      depth_norm = depth_norm / mean(depth_norm)
    ) |> 
  ggplot(aes(x = depth_norm, y = spike_norm)) +
    th() +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, colour = "brown", linetype = "dashed") +
    labs(x = "Spike-in normalisation from spike-in counts", y = "Spike-in normalisation from DiffBind")
}


plot_peak_width_distribution <- function(peaks) {
  env <- new.env(parent = globalenv())
    
  env$d <- peaks |> 
    mutate(width = end - start) |> 
    select(width)

  with(env, {
    d |> 
      ggplot(aes(x = width, y = after_stat(density))) +
      th() +
      geom_histogram(bins = 200) +
      scale_x_log10(labels = scientific_10) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03)), labels = scientific_10) +
      labs(x = "Peak width (bp)", y = "Density")
  })
    
}


plot_total_bindings <- function(db, meta, y_max = NA) {
  env <- new.env(parent = globalenv())
  
  env$y_max <- y_max
  env$d <- extract_bindings(db) |> 
    group_by(sample) |> 
    summarise(tot = sum(value) / 1e6) |> 
    left_join(meta, by = join_by(sample))

  with(env, {
    d |> 
      ggplot(aes(x = group, y = tot, fill = replicate)) +
      th() +
      theme(
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      geom_beeswarm(size = 2, cex = 3, shape = 21, colour = "grey50") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, y_max)) +
      scale_fill_manual(values = okabe_ito_palette) +
      labs(x = NULL, y = "Total normalised abundance")
  })
}

plot_all_total_bindings <- function(db_none, db_libs, db_spike, meta, y_max = 60) {
  pl1 <- plot_total_bindings(db_none, meta, y_max) + ylab(expression(Raw~count~("×"~10^6)))
  pl2 <- plot_total_bindings(db_libs, meta, y_max) + ylab(expression(Library~size~normalised~count~("×"~10^6)))
  pl3 <- plot_total_bindings(db_spike, meta, y_max) + ylab(expression("Spike-in"~normalised~count~("×"~10^6)))

  plot_grid(pl1, pl2, pl3, nrow = 1, labels = c("A", "B", "C"))
}


plot_peak_bindings <- function(db, pids, meta) {
  env <- new.env(parent = globalenv())

  env$d <- db |>
    extract_bindings(pids) |>
    left_join(meta, by = join_by(sample))

  with(env, {
    d |> 
      ggplot(aes(x = group, y = value, fill = replicate)) +
      th() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      geom_beeswarm(size = 2, cex = 3, shape = 21) +
      facet_wrap(~peak_id) +
      scale_fill_manual(values = okabe_ito_palette) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
      labs(x = NULL, y = "Abundance")
  })
}


plot_gene_biotype_distribution <- function(deg, fdr_limit = 0.01) {
  deg |> 
    filter(fdr < fdr_limit) |> 
    group_by(contrast) |> 
    count(gene_biotype) |> 
    mutate(perc = 100 * n / sum(n))  |> 
    ungroup() |> 
  ggplot(aes(x = gene_biotype, y = perc, fill = gene_biotype)) +
    th() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_col(colour = "grey60") +
    scale_fill_manual(values = okabe_ito_palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    facet_wrap(~ contrast) +
    labs(x = "Gene biotype", y = "Percentage")
}


plot_all_genes_biotype <- function() {
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

  entrez_ids <- names(GenomicFeatures::genes(txdb))  # Entrez IDs in this TxDb

  gene_type <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys     = entrez_ids,
    keytype  = "ENTREZID",
    columns  = c("GENETYPE", "SYMBOL", "GENENAME") 
  ) |> 
    as_tibble() |> 
    select(id = ENTREZID, gene_biotype = GENETYPE) |> 
    drop_na()

  gene_type |> 
    distinct(id, gene_biotype) |> 
    add_column(fdr = 0, contrast = "All genes") |> 
    plot_gene_biotype_distribution()
} 


get_pileup_from_peaks <- function(pks, max_dist) {
  regions <- pks |> 
    distinct(grp, start, end, gene_start, distance_to_tss, gene_strand) |> 
    mutate(
      rel_start = if_else(gene_strand == 1, start - gene_start, gene_start - end),
      rel_end = if_else(gene_strand == 1, end - gene_start, gene_start - start)
    ) |> 
    filter(rel_start <= max_dist, rel_end >= -max_dist) |>
    mutate(
      rel_start = pmax(rel_start, -max_dist),
      rel_end   = pmin(rel_end, max_dist),
      end_event = rel_end + 1L
    )
 
  pileup <- bind_rows(
    regions |> transmute(grp, pos = rel_start, delta = 1L),
    regions |> transmute(grp, pos = end_event, delta = -1L)
  ) |>
    group_by(grp, pos) |>
    summarise(delta = sum(delta), .groups = "drop") |>
    group_by(grp) |>
    complete(pos = -max_dist:(max_dist + 1L), fill = list(delta = 0L)) |>
    arrange(pos, .by_group = TRUE) |>
    mutate(coverage = cumsum(delta)) |>
    filter(pos <= max_dist) |>
    ungroup()

  pileup
}

plot_peaks_tss <- function(pks, rand_an = NULL, max_dist = 3000, ncol = 4, group_var = "group") {

  pks <- pks |> 
    mutate(grp = get(group_var))

  if(!is.null(rand_an)) {
    pks <- pks |> 
      bind_rows(rand_an |> mutate(grp = "Random", grp = fct_relevel(grp, "Random", after = Inf)))
  }

  env <- new.env(parent = globalenv())
  env$pileup <- get_pileup_from_peaks(pks, max_dist)
  env$ncol <- ncol
  
  with(env, {
    pileup |>
      ggplot(aes(x = pos, y = coverage)) +
      th() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      geom_vline(xintercept = 0, color = "brown", alpha = 0.3) +
      geom_step() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.04)), limits = c(0, NA)) +
      facet_wrap(~grp, scales = "free_y", ncol = ncol) +
      labs(x = "Position relative to TSS", y = "Coverage")
  })

}


# tst should be created by annot_group_test()
plot_annot_group_test <- function(tst, fdr_limit = 0.01, max_or = 5) {
  sig_tst <- tst |> filter(fdr < fdr_limit)
  mx <- max(abs(tst$odds_ratio))
  if(mx > max_or) mx <- max_or
  tst |> 
    ggplot(aes(x = grp, y = ann, fill = odds_ratio)) +
    th() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    geom_tile() +
    geom_point(data = sig_tst, shape = 8, size = 2, stroke = 0.4) + 
    scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-mx, mx)) +
    labs(x = NULL, y = NULL, fill = "Odds Ratio")
}