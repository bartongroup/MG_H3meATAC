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




