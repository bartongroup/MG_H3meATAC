mn_plot_group_peak_annotations <- function(an, group_var = "group", annot_var = "annotation",
                                        n_top = 30000, rand_an = NULL, palette = 1,
                                        order_annot = FALSE, stacked = TRUE) {
  
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
      mutate(grp = "Random")
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
    rename(prop = estimate)

  if(order_annot)
    dc <- dc |> mutate(ann = fct_reorder(ann, -n))
  
  if(annot_var == "rep_category") {
    dc <- dc |>
      mutate(ann = factor(ann, levels = c("SINE", "LINE", "LTR", "Simple_repeat", "DNA", "Low_complexity", "Satellite", "Retroposon", 
"Other")))
  }
    
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

    pl <- dc |> 
      ggplot(aes(x = grp, y = prop, fill = ann)) +
      th() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_col(position = pos) +
      #scale_fill_manual(values = okabe_ito_palette) +
      scale_fill_manual(values = pal) +
      scale_y_continuous(expand = expansion(mult = c(0, top_margin))) +
      labs(x = NULL, y = "Proportion", fill = "Annotation")
    if(!stacked) {
      pl <- pl + geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = pos, width = 0.2)
    }
    
  pl

}


mn_plot_volcano <- function(res, fdr_limit = 0.01, logfc_limit = 0, point_size = 0.5,
                            point_alpha = 0.5, scales = "fixed", ncol = NULL) {
  
  r <- res |>
    mutate(
      x = log_fc,
      y = -log10(p_value),
      group = contrast,
      up = fdr < fdr_limit & log_fc > logfc_limit,
      down = fdr < fdr_limit & log_fc < -logfc_limit,
    ) |>
    select(x, y, up, down, group)

  r_up <- r |> filter(up)
  r_down <- r |> filter(down)
  r_nsel <- r |> filter(!up & !down)
  
  g <- res |>
    ggplot() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(data = r_nsel, aes(x = x, y = y), colour = "grey70",
              size = point_size, alpha = point_alpha) +
    facet_wrap( ~ group, scale = scales, ncol = ncol) +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA))


    if (nrow(r_up) > 0) {
      g <- g + geom_point(data = r_up, aes(x = x, y = y), colour = "black",
                          size = point_size, alpha = point_alpha)
    }

    if (nrow(r_down) > 0) {
      g <- g + geom_point(data = r_down, aes(x = x, y = y), colour = "black",
                          size = point_size, alpha = point_alpha)
    }

  g
}



mn_plot_annot_group_test <- function(tst, fdr_limit = 0.01, log_or_limit = 0.5, max_or = NULL) {
    if("SINE" %in% tst$ann) {
    tst <- tst |>
      mutate(ann = factor(ann, levels = c("SINE", "LINE", "LTR", "Simple_repeat", "DNA", "Low_complexity", "Satellite", "Retroposon", 
"Other")) |> fct_rev())
  }

  sig_tst <- tst |> filter(fdr < fdr_limit & abs(log_odds_ratio) > log_or_limit)
  inf_tst <- tst |> filter(is.infinite(log_odds_ratio))
  mx <- max(abs(tst$log_odds_ratio))
  if(!is.null(max_or)){
    mx <- max_or
    tst <- tst |> 
      mutate(log_odds_ratio = if_else(abs(log_odds_ratio) > mx, sign(log_odds_ratio) * mx, log_odds_ratio))
  } 
  tst |> 
    ggplot(aes(x = grp, y = ann, fill = log_odds_ratio)) +
    th() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    geom_tile() +
    geom_point(data = sig_tst, shape = 8, size = 2, stroke = 0.4) + 
    geom_text(data = inf_tst, label = expression(infinity), size = 5, colour = "white", alpha = 0.5) +
    scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-mx, mx)) +
    labs(x = NULL, y = NULL, fill = expression(log[2]~OR))
}


mn_get_enr_contrasts_logp <-  function(tf_enr, ctr1, ctr2) {
  d <- tf_enr |> 
    mutate(lp = -log10(pval)) |>
    filter(contrast_direction %in% c(ctr1, ctr2)) 
  
  d |> 
    pivot_wider(id_cols = c(motif_name, motif_consensus), names_from = contrast_direction, values_from = lp) |> 
    rename(x = !!sym(ctr1), y = !!sym(ctr2)) |> 
    mutate(
      mrg = is.na(x) | is.na(y),
      across(x:y, ~replace_na(.x, 0))
    ) |> 
    distinct(x, y, motif_name, mrg)
}

mn_plot_enr_contrasts_logp <- function(tf_enr, ctr1, ctr2) {
  mn_get_enr_contrasts_logp(tf_enr, ctr1, ctr2) |>
    ggplot(aes(x = x, y = y, shape = mrg)) +
    th() +
    theme(legend.position = "none") +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = motif_name), size = 2) +
    scale_shape_manual(values = c(16, 1)) +
    labs(x = bquote(-log[10] ~ P ~ (.(ctr1))), y = bquote(-log[10] ~ P ~ (.(ctr2))))
}


merge_matsuda <- function(matsuda, dec_atac, dnr) {
  dem <- de_matsuda(matsuda, dnr)
  m <- matsuda |> 
    filter(donor == dnr & condition %in% c("iPS", "PS", "PSM") & value > 0) |> 
    group_by(gene_id, gene, group, condition) |> 
    summarise(v = mean(value), .groups = "drop") |> 
    pivot_wider(id_cols = c(gene_id, gene, group), names_from = condition, values_from = v) |> 
    mutate(smat_logfc = log2(0.5 * (PS + PSM) / iPS)) |>
    rename(gene_symbol = gene) |> 
    left_join(dem, by = join_by(gene_id)) |> 
    rename(mat_logfc = logFC) |> 
    filter(adj.P.Val < 0.01 & abs(mat_logfc) > 1)
  d <- dec_atac |> 
    filter(fdr < 0.01 & abs(log_fc) > 1 & contrast == "48h_Ctrl-iPS_Ctrl" & annotation %in% c("Promoter", "1-5 kb")) |>
    select(gene_symbol, atac_logfc = log_fc)
  # Decided not to average
    #drop_na() |> 
    #group_by(gene_symbol) |> 
    #summarise(atac_logfc = mean(atac_logfc))

  inner_join(m, d, by = join_by(gene_symbol)) |> 
    select(gene_id, gene_symbol, group, PS, PSM, iPS, mat_logfc, p_value = P.Value, fdr = adj.P.Val, atac_logfc)
}

mn_plot_matsuda_comparison <- function(matsuda, dec_atac, dnr) {
 
  df <- merge_matsuda(matsuda, dec_atac, dnr)

  df |> 
    ggplot(aes(x = mat_logfc, y = atac_logfc, colour = -log10(fdr))) +
    th() +
    geom_point() +
    geom_vline(xintercept = 0, colour = "grey80") +
    geom_hline(yintercept = 0, colour = "grey80") +
    scale_colour_distiller(palette = "Reds", direction = 1) +
    labs(x = expression(log[2]~(0.5*(PS + PSM) / iPS)), y = expression(log[2]~FC~"ATAC-seq"), colour = expression(-log[10]~FDR))
}