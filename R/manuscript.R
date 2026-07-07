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
                            point_alpha = 0.5, scales = "fixed") {
  
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
    facet_wrap( ~ group, scale = scales) +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA))


    if (nrow(r_up) > 0) {
      g <- g + geom_point(data = r_up, aes(x = x, y = y), colour = "#0072B2",
                          size = point_size, alpha = point_alpha)
    }

    if (nrow(r_down) > 0) {
      g <- g + geom_point(data = r_down, aes(x = x, y = y), colour = "#D55E00",
                          size = point_size, alpha = point_alpha)
    }

  g
}
