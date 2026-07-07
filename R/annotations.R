get_functional_terms <- function(config, geneids) {
  SPECIES <- list(
    human = list(
      go = "goa_human",
      re = "Homo sapiens",
      kg = "hsa",
      wi = "Homo sapiens"
    ),
    yeast = list(
      go = "sgd",
      re =  "Saccharomyces cerevisiae",
      kg = "sce"
    ),
    rat = list(
      go = "rgd",
      re = "Rattus norvegicus",
      kg = "rno"
    ),
    mouse = list(
      go = "mgi",
      re = "Mus musculus",
      kg = "mmu"
    )
  )
  
  geneids <- geneids |> 
    select(id, gene_symbol) |> 
    distinct()

  sp <- SPECIES[[config$species]]
  cat("Loading GO term data\n")
  go <- fenr::fetch_go(dataset = config$ensembl$dataset)
  cat("Loading Reactome data\n")
  re <- fenr::fetch_reactome(sp$re)
  cat("Loading KEGG data\n")
  kg <- fenr::fetch_kegg(sp$kg)
  #cat("Loading BioPlanet data\n")
  #bp <- fetch_bp()
  
  go$terms <- go$terms |> 
    filter(namespace == "biological_process")
  go$mapping <- go$mapping |> 
    drop_na() |> 
    rename(id = gene_id) |> 
    filter(term_id %in% go$terms$term_id)
  re$mapping <- re$mapping |> 
    rename(id = gene_id)
  kg$mapping <- kg$mapping |> 
    select(-gene_id) |> 
    left_join(geneids, by = "gene_symbol", relationship = "many-to-many")
  #bp$mapping <- bp$mapping |>
  #  left_join(geneids, by = "gene_symbol", relationship = "many-to-many")

  list(
    go = go,
    re = re,
    kg = kg
   # bp = bp
  )
}


prepare_terms_fenr <- function(terms, all_features = NULL, feature_name = "id") {
  ontologies <- names(terms)
  
  map(ontologies, function(ont) {
    trm <- terms[[ont]]
    fenr::prepare_for_enrichment(trm$terms, trm$mapping, all_features = all_features, feature_name = feature_name)
  }) |> 
    set_names(ontologies)
}




fetch_range_sequences <- function(fasta, ranges) {
  if(!"peak_id" %in% names(ranges)) {
    ranges <- ranges |>
      dplyr::mutate(peak_id = stringr::str_glue("{chr}:{start}-{end}"))
  }

  if(!"strand" %in% names(ranges)) {
    ranges <- ranges |>
      dplyr::mutate(strand = "+")
  }

  fa <- if(inherits(fasta, "FaFile")) {
    fasta
  } else {
    Rsamtools::FaFile(fasta)
  }

  seqlevels_fa <- Rsamtools::scanFaIndex(fa) |>
    GenomeInfoDb::seqnames() |>
    as.character()

  d <- ranges |>
    dplyr::filter(chr %in% seqlevels_fa) |>
    dplyr::mutate(
      strand = as.character(strand),
      strand = dplyr::case_when(
        strand %in% c("-", "-1", "-1L") ~ "-",
        TRUE ~ "+"
      )
    )

  missing_chr <- setdiff(unique(d$chr), seqlevels_fa)
  if(length(missing_chr) > 0) {
    stop("Chromosomes not found in FASTA: ", paste(missing_chr, collapse = ", "))
  }

  gr <- GenomicRanges::makeGRangesFromDataFrame(d)
  seqs <- Rsamtools::scanFa(fa, param = gr)

  rev <- d$strand == "-"
  seqs[rev] <- Biostrings::reverseComplement(seqs[rev])

  d |>
    dplyr::mutate(sequence = as.character(seqs))
}