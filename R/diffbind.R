make_diffbind_metadata <- function(meta, config, type = "broad") {
  if(type == "broad") {
    macs_path <- config$dir$macs_broad
    suffix <- "peaks.broadPeak"
  } else {
    macs_path <- config$dir$macs_narrow
    suffix <- "peaks.narrowPeak"
  }
  
  meta |> 
    transmute(
      SampleID = sample,
      Tissue = time_point,
      Factor = "H3K9me3",
      Condition = group |> str_remove("^Me._"),
      Treatment = treatment,
      Replicate = replicate,
      bamReads = file.path(config$remote_path, config$dir$bam, str_glue("{raw_sample}.noMT.bam")),
      Spikein = file.path(config$remote_path, config$dir$bam_dm, str_glue("{raw_sample}.bam")),
      Peaks = file.path(config$remote_path, macs_path, str_glue("{raw_sample}_{suffix}")),
      PeakCaller = "bed",
      ScoreCol = 5
    )
}

make_diffbind_metadata_for_cluster <-  function(meta, config, type = "broad") {
  config$remote_path <- "."
  make_diffbind_metadata(meta, config, type) |> 
    write_csv(str_glue("diffbind_metadata_{type}.csv"))
}

read_and_count_diffbind <- function(meta_db, min_overlap = 2) {
  
  ## Read all peaks & metadata 
  db <- DiffBind::dba(sampleSheet = meta_db)
  
  ## Build per-condition consensus peaksets
  peakset.cons <- DiffBind::dba.peakset(
    db,
    consensus = DiffBind::DBA_CONDITION,
    minOverlap = min_overlap
  )
  
  ## Restrict to consensus peaks and extract as GRanges
  peakset.cons <- DiffBind::dba(
    peakset.cons,
    mask = peakset.cons$masks$Consensus,
    minOverlap = 1
  )
  peaks.cons <- DiffBind::dba.peakset(peakset.cons, bRetrieve = TRUE)

  # Make sure multimappers are not silently dropped
  db$config$mapQCth <- 0

  db <- DiffBind::dba.count(
    db,
    peaks = peaks.cons,
    summits = FALSE,
    bRemoveDuplicates = FALSE,
    bUseSummarizeOverlaps = TRUE,
    bParallel = FALSE  # alas, parallel keep crashing
  )
  db
}

spikein_normalization <- function(dbc) {
  DiffBind::dba.normalize(dbc, normalize = DiffBind::DBA_NORM_LIB, spikein = TRUE)
}

# AFAIK this doesn't do anything, data seem to be normalised
# to the library size be default
libsize_normalization <- function(dbc) {
  DiffBind::dba.normalize(dbc, normalize = DiffBind::DBA_NORM_LIB)
}

# "fake" normalisation to the total count of 1e7 - effectively removes all normalisations
# and returns an object with original counts per peak
no_normalization <- function(dbc, tot_count = 1e7) {
  ns <- nrow(dbc$samples)
  DiffBind::dba.normalize(dbc, normalize = DiffBind::DBA_NORM_LIB, library = rep(tot_count, ns))
}



# Need spike-in normalized data (slow process)
# Default libsize normalization is calculated on the fly
get_norm_factors <- function(dbcs) {
  
  samples <- DiffBind::dba.show(dbcs)$ID

  nsp <- DiffBind::dba.normalize(dbcs, bRetrieve = TRUE)
  spikein_lib_sizes <- nsp$lib.sizes
  spikein_norm_factors <- nsp$norm.factors

  dbc <- DiffBind::dba.normalize(dbcs, normalize = DiffBind::DBA_NORM_LIB)
  nf <- DiffBind::dba.normalize(dbc, bRetrieve = TRUE)
  lib_sizes <- nf$lib.sizes
  norm_factors <- nf$norm.factors

  tibble(
    sample = samples,
    spikein_lib_size = spikein_lib_sizes,
    spikein_norm_factor = spikein_norm_factors,
    lib_size = lib_sizes,
    norm_factor = norm_factors
  )
}


de_full_model <- function(dbc) {
  dbc <- DiffBind::dba.contrast(dbc, contrast = c("Tissue", "iPS", "48h"))
  dbc <- DiffBind::dba.contrast(dbc, contrast = c("Treatment", "Ctrl", "SETDB1"))
  dbc <- DiffBind::dba.contrast(dbc, design = "~ Tissue + Treatment")
  DiffBind::dba.analyze(dbc, bBlacklist = FALSE, bGreylist = FALSE)
}

de_groups <- function(dbc) {
  dbc <- DiffBind::dba.contrast(dbc, contrast = c("Condition", "iPS_Ctrl", "48h_Ctrl"))
  dbc <- DiffBind::dba.contrast(dbc, contrast = c("Condition", "48h_Ctrl", "48h_SETDB1"))
  dbc <- DiffBind::dba.contrast(dbc, contrast = c("Condition", "iPS_Ctrl", "iPS_SETDB1"))
  dbc <- DiffBind::dba.contrast(dbc, contrast = c("Condition", "iPS_SETDB1", "48h_SETDB1")) 
  # dbc <- DiffBind::dba.contrast(dbc, contrast = c("Condition", "iPS_Ctrl", "48h_SETDB1"))
  # dbc <- DiffBind::dba.contrast(dbc, contrast = c("Condition", "iPS_SETDB1", "48h_Ctrl"))
  DiffBind::dba.analyze(dbc, bBlacklist = FALSE, bGreylist = FALSE)
}

de_model <- function(dbc, model) {
  if(model == "full") {
    de_full_model(dbc)
  } else if (model == "groups") {
    de_groups(dbc)
  } else {
    stop("Unknown model: ", model)
  }
}

extract_consensus <- function(db) {
  db |> 
    DiffBind::dba.peakset(bRetrieve = TRUE) |> 
    as_tibble() |> 
    select(chr = seqnames, start, end, width)
}

extract_bindings <- function(db, pids = NULL) {
  d <- db |> 
    DiffBind::dba.peakset(bRetrieve = TRUE) |> 
    as_tibble() |> 
    mutate(peak_id = str_glue("{seqnames}:{start}-{end}"), .before = 1) |> 
    select(-c(seqnames, start, end, width, strand)) |> 
    pivot_longer(-peak_id, names_to = "sample", values_to = "value")
  if(!is.null(pids)) {
    d <- d |> 
      filter(peak_id %in% pids)
  }
  d
}

extract_de_results <- function(dbc) {
  contrasts <- DiffBind::dba.show(dbc, bContrasts = TRUE) |> 
    as_tibble() |> 
    mutate(contrast = str_glue("{Group2}-{Group}"))
  map(1:nrow(contrasts), function(i) {
    contrast_name <- contrasts$contrast[i]
    DiffBind::dba.report(dbc, contrast = i, th = 1) |> 
      as_tibble() |> 
      rename(conc_1 = 7, conc_2 = 8) |> 
      mutate(contrast = contrast_name)
  }) |> 
    bind_rows() |> 
    rename(
      chr = seqnames,
      conc = Conc,
      log_fc = Fold,
      p_value = p.value,
      fdr = FDR
    ) |> 
    mutate(
      log_fc = -log_fc,  # They provide group1 - group2, but we prefer group2 - group1
      direction = if_else(log_fc > 0, "up", "down"),
      peak_id = str_glue("{chr}:{start}-{end}"), .before = 1
    ) |> 
    unite(contrast_direction, contrast, direction, sep = "_", remove = FALSE)
}

prepare_de_genes <- function(dec, keep_annot = c("Intron", "Exon", "Promoter (<=1kb)", "3' UTR", "5' UTR")) {
  dec |> 
    select(peak_id, id = gene_id, contrast, log_fc, p_value, fdr, annotation, gene_biotype) |>
    filter(!is.na(id) & annotation %in% keep_annot)
}


read_correct_cache <- function(type = "broad", spike = TRUE) {
  s <- ifelse(spike, "_spike", "")
  d <- readr::read_rds(paste0("cache/dbc_", type, s, ".rds"))
  d$class[c("bamRead", "Spikein"), ] <- d$class[c("bamRead", "Spikein"), ] |>
    str_replace("./bam", "/Volumes/mgierlinski/projects/h3mecut/bam")
  d
}


read_cache <- function(type, norm) {
  if(norm == "spike") {
    read_correct_cache(type, spike = TRUE)
  } else if(norm == "libs") {
    read_correct_cache(type, spike = FALSE) |> 
      libsize_normalization() 
  } else if(norm == "none") {
    read_correct_cache(type, spike = FALSE) |> 
      no_normalization()
  } else {
    stop("Huh?")
  }
}