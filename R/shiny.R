write_shiny_file <- function(obj, method = c("qs", "rds")) {
  method <- match.arg(method)
  path <- file.path("shiny", "data")
  if (!dir.exists(path)) dir.create(path)
  obj_name <- deparse(substitute(obj))
  cat(paste("  Writing", obj_name))
  if(method == "qs") {
    file_name <- file.path(path, str_glue("{obj_name}.qs"))
    qs2::qs_save(obj, file_name)
  } else if(method == "rds") {
    file_name <- file.path(path, str_glue("{obj_name}.rds"))
    write_rds(obj, file_name, compress = "xz")
  }
  cat("\n")
}

make_dbc_dexdash <- function(dbc, meta, name) {
  abu <- extract_bindings(dbc) |> 
    rename(id = peak_id, abundance = value)
  dec <- extract_de_results(dbc) |> 
    rename(id = peak_id, expr = conc)

  dexdash::dexdash_set(dec, abu, meta, name)
}

map_peak_to_gene <- function(dbc) {
  # Make peak_id - gene_id mapping
  dec <- extract_de_results(dbc)
  p2g <- dec |> 
    annotate_peaks_chipseeker() |> 
    prepare_de_genes() |> 
    distinct(peak_id, gene_id = id)
  dec |> 
    distinct(peak_id) |> 
    left_join(p2g, by = join_by(peak_id)) |> 
    rename(id = peak_id)
}

make_dbc_annotations <- function(p2g, terms) {
  ontologies <- names(terms)
  pterms <- map(ontologies, function(ont) {
    maps <- p2g |> 
      inner_join(terms[[ont]]$mapping |> rename(gene_id = id), by = join_by(gene_id), relationship = "many-to-many") |> 
      distinct(id, term_id)
    list(
      terms = terms[[ont]]$terms,
      mapping = maps
    )
  }) |> 
    set_names(ontologies)
  
  prepare_terms_fenr(pterms)
}

prepare_for_shiny <- function(dbc_full_libs_broad, dbc_full_spike_broad, dbc_groups_libs_broad, dbc_groups_spike_broad,
                              meta, terms, gns) {

  # Peak positions and IDs are the same in all sets
  p2g <- map_peak_to_gene(dbc_full_libs_broad)
  fterms <- make_dbc_annotations(p2g, terms)

  gn <- gns |> 
    distinct(gene_id = id, name = gene_symbol, description) |> 
    unite(description, c(name, description), sep = ": ")
  features <- p2g |> 
    left_join(gn, by = join_by(gene_id)) |>
    mutate(name = id, .after = id)
  meat <- meta |> 
    select(sample, group, time_point, treatment, replicate)
 
  dexset <- dexdash::dexdash_list(
    make_dbc_dexdash(dbc_full_libs_broad, meta, "Full libsize"),
    make_dbc_dexdash(dbc_full_spike_broad, meta, "Full spikein"),
    make_dbc_dexdash(dbc_groups_libs_broad, meta, "Pairs libsize"),
    make_dbc_dexdash(dbc_groups_spike_broad, meta, "Pairs spikein")
  )

  write_shiny_file(dexset)
  write_shiny_file(features)
  write_shiny_file(fterms)
}