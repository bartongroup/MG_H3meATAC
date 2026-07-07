read_fastq_length <- function(path, meta) {
  map2_dfr(meta$raw_sample, meta$sample, function(rsam, sam) {
    cfile <- file.path(path, str_glue("{rsam}_count.txt"))
    if(file.exists(cfile)) {
      read_tsv(cfile, col_names = "read_count", show_col_types = FALSE) |> 
        mutate(sample = sam, raw_sample = rsam, .before = 1)
    }
  })
}

read_fastq_lengths <- function(config, meta) {
  trimmed_path <- file.path(config$data_path, config$dir$fastq_trimmed)
  f_trim <- read_fastq_length(trimmed_path, meta) |> 
    add_column(type = "Trimmed")
  clean_path <- file.path(config$data_path, config$dir$fastq_clean)
  f_clean <- read_fastq_length(clean_path, meta) |> 
    add_column(type = "Clean")
  bind_rows(f_trim, f_clean)
}


parse_one_fscreen <- function(file, raw_sample) {
  read_tsv(file, skip = 1, comment = "%Hit", col_types = cols()) |> 
    pivot_longer(-Genome) |> 
    mutate(raw_sample = raw_sample)
}

parse_fscreens <- function(config, meta, suffix = "_R1_screen.txt") {
  path <- file.path(config$data_path, config$dir$fastq_screen)
  meta$raw_sample |> 
    map(function(rs) {
      file <- file.path(path, str_c(rs, suffix))
      parse_one_fscreen(file, rs)
    }) |> 
    list_rbind() |> 
    left_join(select(meta, raw_sample, sample), by = "raw_sample")
}

plot_fscreen_map <- function(fscreen) {
  fscreen |> 
    filter(name == "%One_hit_one_genome") |> 
    mutate(value = na_if(value, 0)) |> 
    mutate_at(c("Genome", "sample"), as_factor) |> 
    ggplot(aes(x = sample, y = Genome, fill = value)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
    geom_tile() +
    scale_fill_viridis_c(option = "cividis", trans = "log10", limits = c(0.01, 100), labels = c(0.01,0.1,1,10,100)) +
    labs(x = NULL, y = "Genome", fill = "% Unique hits")
}


plot_fscreen_sample <- function(fscreen, smpl) {
  fscreen |>
    filter(sample == smpl & str_detect(name, "%") & !str_detect(name, "Unma")) |>
    mutate(Genome = as_factor(Genome) |> fct_reorder(value, sum)) |> 
    mutate(name = str_replace_all(name, "_", " ")) |> 
    ggplot(aes(x = Genome, y = value, fill = name)) +
    geom_col() +
    scale_fill_manual(values = okabe_ito_palette[c(1,2,3,5)]) +
    coord_flip() +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Genome", y = "Percentage", title = smpl, fill = "Quantity")
  
}

# for paired
make_roots_paired <- function(meta, suffix = "_R") {
  roots_1 <- paste0(meta$raw_sample, suffix, "1")
  roots_2 <- paste0(meta$raw_sample, suffix, "2")
  tibble(
    root = c(roots_1, roots_2),
    raw_sample = rep(meta$raw_sample, 2),
    sample = rep(meta$sample, 2),
    pair = c(rep(1, nrow(meta)), rep(2, nrow(meta))) |> as_factor()
  )
}

# for single
make_roots_single <- function(meta) {
  roots <- meta$raw_sample
  tibble(
    root = roots,
    raw_sample = meta$raw_sample,
    sample = meta$sample,
    pair = rep(1, nrow(meta)) |> as_factor()
  )
}

parse_one_qc <- function(zip_file, root) {
  if (!file.exists(zip_file)) {
    warning(paste("File", zip_file, "not found."))
    return(NULL)
  }
  dirname <- str_remove(basename(zip_file), ".zip")
  con <- unz(zip_file, file.path(dirname, "fastqc_data.txt"), open = "rb")
  s <- read_lines(con)
  close(con)
  p1 <- str_which(s, "Per base sequence quality")[1] + 1
  p2 <- str_which(s, "END_MODULE")
  p2 <- min(p2[p2 > p1]) - 1
  st <- s[p1:p2] |> str_remove("#")
  read_tsv(I(st), show_col_types = FALSE, progress = FALSE) |>
    mutate(Base = str_remove(Base, "\\-.*$") |> as.integer()) |> 
    mutate(root = root)
}

parse_qcs <- function(config, meta, paired, suffix = "_fastqc.zip") {
  path <- file.path(config$data_path, config$dir$qc)
  if (paired) {
    roots <-  make_roots_paired(meta)
  } else {
    roots <-  make_roots_single(meta)
  }
  roots$root |> 
    map_dfr(~parse_one_qc(file.path(path, paste0(.x, suffix)), .x)) |> 
    left_join(roots, by = "root")
}


plot_fscreen_slog <- function(fs, slog, genome, fs_var = "%One_hit_one_genome", slog_var = "Uniquely mapped reads %") {
  fsf <- fs |> 
    filter(name == fs_var & Genome == genome) |> 
    group_by(raw_sample) |> 
    summarise(screen = mean(value))
  slogf <- slog |> 
    filter(key == slog_var) |> 
    mutate(value = str_remove(value, "%") |> as.numeric()) |> 
    rename(star = value)
  fsf |> 
    left_join(slogf, by = "raw_sample") |> 
  ggplot(aes(x = star, y = screen)) +
    theme_bw() +
    geom_point() +
    labs(x = slog_var, y = fs_var)
}

plot_ribo_proportion <- function(fl) {
  fl |> 
    pivot_wider(id_cols = sample, names_from = type, values_from = read_count) |> 
    mutate(prop = 1 - Clean / Trimmed) |> 
  ggplot(aes(x = sample, y = prop)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    geom_col() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = NULL, y = "Ribosomal proportion")
}




plot_flow_qc <- function(slog, qc, meta, base = 20) {
  qcf <- qc |> 
    filter(Base == base) |> 
    group_by(raw_sample) |> 
    summarise(quality = mean(Mean)) 
  slogf <- slog |> 
    filter(key == "Uniquely mapped reads %") |> 
    mutate(mapped_prop = str_remove(value, "%") |> as.numeric()) |> 
    select(raw_sample, mapped_prop)
  meta |> 
    left_join(qcf, by = "raw_sample") |> 
    left_join(slogf, by = "raw_sample") |> 
    pivot_longer(cols = c(quality, mapped_prop)) |> 
  ggplot(aes(x = value, y = flow_cell)) +
    theme_bw() +
    geom_beeswarm(groupOnX = FALSE, size = 0.5, cex = 0.5) +
    facet_wrap(~name, scales = "free_x") +
    labs(x = NULL, y = "Flow cell")
}

# read file consisting of first fastq lines, created with
# for file in *gz; do printf "$file "; gunzip -c $file | head -n 1; done > fastq_headers.txt
read_fastq_headers <- function(file) {
  read_delim(file, " ", col_names = c("file_name", "id", "attr"), col_types = "ccc") |> 
    mutate(raw_sample = str_remove(file_name, "_[12].fq.gz")) |> 
    mutate(accession = str_remove(id, "@") |> str_remove("\\.1")) |> 
    separate(attr, c("flow_cell", "lane", "tile", "x", "yr"), sep = ":") |> 
    separate(yr, c("y", "pair"), sep = "/") |> 
    group_by(raw_sample) |> 
    summarise(flow_cell = first(flow_cell) |> as_factor(), accession = first(accession) |> as.factor())
}


read_houskeeping_genes <- function(path, classes = c("Ribosomal", "Citric", "RNA")) {
  map_dfr(classes, function(cl) {
    file <- file.path(path, str_glue("protein_class_{cl}.tsv"))
    stopifnot(file.exists(file))
    read_tsv(file, show_col_types = FALSE) |>
      select(gene_symbol = Gene) |>
      add_column(class = cl, .before = 1)
  })
}

parse_idxstats <- function(config, meta) {
  path <- file.path(config$data_path, config$dir$chr_count)
  map2_dfr(meta$raw_sample, meta$sample, function(rsam, sam) {
    sfile <- file.path(path, glue::glue("{rsam}.txt"))
    if (file.exists(sfile)) {
      read_tsv(sfile, col_names = c("chr", "length", "count", "unmapped"), show_col_types = FALSE) |> 
        add_column(raw_sample = rsam, sample = sam)
    }
  })
}

plot_chrom_proportion <- function(ids, chromosomes, per_length = TRUE) {
  d <- ids |>
    filter(chr %in% chromosomes) 
  
  if(per_length) {
    d <- d |> mutate(count = count / (length / 1000))
    ylab <- "Counts per kb of chromosome"
  } else {
    d <- d |> mutate(count = count / 1e6)
    ylab <- "Counts per chromosome (millions)"
  }
  
  d |> 
    mutate(chr = factor(chr, levels = chromosomes)) |> 
    ggplot(aes(x = chr, y = count)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 7)
    ) +
    geom_col() +
    facet_wrap(~ sample, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Chromosome", y = ylab)
}

plot_mitochondrial_proportion <- function(ids, meta) {
  ids |> 
    mutate(MT = if_else(chr == "MT", "MT", "not_MT")) |> 
    group_by(sample, MT) |> 
    summarise(count = sum(count)) |> 
    ungroup() |> 
    pivot_wider(id_cols = sample, names_from = MT, values_from = count) |> 
    mutate(MT_perc = 100 * MT / (MT + not_MT)) |> 
    left_join(meta, by = join_by(sample)) |> 
  ggplot(aes(x = group, y = MT_perc, fill = replicate)) +
    th() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_beeswarm(size = 2, cex = 3, shape = 21, colour = "grey50") +
    scale_fill_manual(values = okabe_ito_palette) +   
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA), trans = "sqrt", breaks = c(1, 2, 5, 10, 20, 40, 60)) +
    labs(x = "Group", y = "% mitochondrial reads")

}

get_software_versions <- function(ver_file, selection) {
  read_table(ver_file, col_names = c("name", "version", "hash", "source"), comment = "#", show_col_types = FALSE) |>
    filter(name %in% selection)
}



plot_qualities <- function(qc) {
  qc |> 
    ggplot(aes(x = Base, y = Mean, group = root, colour = pair)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_line(alpha = 0.3, linewidth = 0.6) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = "Base", y = "Mean quality score")
}

plot_cluster_qualities <- function(qc, text.size = 10) {
  tab <- qc |> 
    pivot_wider(id_cols = c(pair, Base), names_from = sample, values_from = Mean) |> 
    select(-c(pair, Base)) |> 
    as.matrix()
  
  hc <- t(tab) |> 
    dist() |> 
    hclust()
  
  dendr <- ggdendro::dendro_data(hc)
  seg <- ggdendro::segment(dendr)
  theme.d <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = text.size),
    axis.line.y = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(linewidth = 0.5),
    axis.ticks.y = ggplot2::element_blank()
  )
  ggplot() +
    theme.d +
    coord_flip() +
    geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_continuous(breaks = seq_along(dendr$labels$label), labels = dendr$labels$label) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(seg$y) * 1.03)) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = "Distance")
  
}


read_spikein_norms <- function(config, meta) {
  path <- file.path(config$data_path, config$dir$spike_norm)

  map2(meta$sample, meta$raw_sample, function(smpl, rsmpl) {
    scale_file = file.path(path, str_glue("{rsmpl}.spikein_scale.txt"))
    depth_file = file.path(path, str_glue("{rsmpl}.spikein_depth.txt"))
    tibble(
      sample = smpl,
      depth = read_lines(depth_file) |> as.numeric(),
      scale = read_lines(scale_file) |> as.numeric()
    ) 
  }) |> 
  list_rbind()
}

plot_spikein_norms <- function(sn, meta) {
  sn |> 
    mutate(depth = depth / 1e3) |> 
    rename(
      `Sequencing depth (x1000)` = depth,
      `Scaling factor` = scale
    ) |> 
    pivot_longer(-sample) |> 
    left_join(meta, by = join_by(sample)) |> 
    ggplot(aes(x = group, y = value, fill = replicate)) +
    th() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_beeswarm(size = 2, cex = 3, shape = 21, colour = "grey50") +
    scale_fill_manual(values = okabe_ito_palette) +   
    facet_wrap(~name, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) +
    labs(x = "Group", y = "Value")
}