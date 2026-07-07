' Run fgsea on Term Data
#'
#' This function performs fgsea analysis on a given set of term data and DE results.
#'
#' @param term_data A list containing term-to-feature mappings and term-to-name mappings.
#' @param res A dataframe containing DE results with feature IDs and values.
#' @param min_size Minimum size of a term for fgsea analysis. Default is 3.
#' @param n_perm Number of permutations for fgsea analysis. Default is 1000.
#'
#' @return A tibble with the fgsea results including term IDs, term names, p-values, FDR, NES,
#' size, and leading edge.
fgsea_run <- function(term_data, res, min_size = 3, n_perm = 1000) {
  res <- res |>
    filter(!is.na(value) & !is.na(feature_id)) |> 
    mutate(first_id = str_remove(feature_id, ";.*$"))
  ranks <-  set_names(res$value, res$first_id)
  fgsea::fgsea(pathways = as.list(term_data$term2feature), stats = ranks, nproc = 6, minSize = min_size, eps = 0, nPermSimple = n_perm) |>
    as_tibble() |> 
    rename(term_id = pathway) |> 
    mutate(term_name = as.character(as.list(term_data$term2name)[term_id])) |> 
    arrange(NES) |>
    select(term_id, term_name, p_value = pval, fdr = padj, nes = NES, size, leading_edge = leadingEdge)
}

#' Run fgsea for all ontologies
#'
#' This function performs fgsea analysis for all provided ontologies and DE results.
#'
#' @param de A dataframe containing DE results.
#' @param fterms A list of `fenr_terms` objects containing feature terms for various ontologies.
#' @param filt An expression used to filter the DE results. Default is "TRUE".
#' @param feature_var Column name in `de` containing feature IDs. Default is "id".
#' @param rank_expr An expression to create a rank variable. We suggest "logFC" (default)
#'   or "-sign(logFC) * log10(PValue)".
#' @param group_var Column name in `de` containing group information. Default is "contrast".
#' @param min_size Minimum size of a term for fgsea analysis. Default is 3.
#' @param n_perm Number of permutations for fgsea analysis. Default is 1000.
#'
#' @return A named list of tibbles with fgsea results for each ontology.
fgsea_all_terms <- function(de, fterms, filt = "TRUE", feature_var = "id", rank_expr = "logFC",
                            group_var = "contrast", min_size = 3, n_perm = 1000, seed = 42) {
  ontologies <- names(fterms)
  de <- de |> 
    mutate(
      value = !!rlang::parse_expr(rank_expr),
      feature_id = get(feature_var)
    ) |> 
    filter(!!rlang::parse_expr(filt)) |> 
    group_by(feature_id, contrast) |>
    # special for DiffBind results, multiple peaks per gene
    # here we average the rank expression
    summarise(value = mean(value), .groups = "drop") 
  

  set.seed(seed)
  map(ontologies, function(ont) {
    cat(str_glue("  Computing fgsea for {ont}\n\n"))
    de |>
      group_split(!!sym(group_var)) |> # this is usually contrast
      map(function(w) {
        fgsea_run(fterms[[ont]], w, min_size, n_perm) |>
          mutate(!!group_var := dplyr::first(w[[group_var]]))
      }) |> 
      list_rbind() |> 
      add_column(ontology = ont)
  }) |>
    list_rbind()
}

tabulate_gsea <- function(gse, gns, fdr_limit = 0.05, max_n = 12) {
  g2s <- gns |> 
    distinct(id, gene_symbol)
  gse |> 
    filter(fdr < fdr_limit) |>
    unnest(leading_edge) |> 
    rename(id = leading_edge) |> 
    left_join(g2s, by = "id") |> 
    select(-id) |> 
    group_by(ontology, term_id, term_name, contrast, nes, fdr, size) |>
    summarise(genes = gene_symbol |> unique() |> head(max_n) |> str_c(collapse = ", ")) |> 
    ungroup()
}



#' Plot fgsea Enrichment
#'
#' This function plots the fgsea enrichment for a given term and DE results.
#'
#' @param de A dataframe containing DE results.
#' @param fterms A list of feature terms.
#' @param info A list "contrast", "ontology" and "term_id" to display.
#' @param value_var Column name in `de` containing values. Default is "logFC".
#'
#' @return A ggplot object showing the fgsea enrichment plot.
plot_fgsea_enrichment <- function(de, fterms, info, rank_expr = "logFC") {
  de <- de |>
    filter(contrast == info$contrast) |> 
    mutate(value = !!rlang::parse_expr(rank_expr))
  lst <- as.list(fterms[[info$ontology]])$term2feature[[info$term_id]]
  rnks <- set_names(de$value, de$id)
  fgsea::plotEnrichment(lst, rnks)
}

#' Split Genes by fgsea Results
#'
#' This function splits genes based on fgsea results and DE results.
#'
#' @param se A dataframe containing DE results.
#' @param fg A dataframe containing fgsea results.
#' @param groupvar Column name in `fg` containing group information. Default is "contrast".
#'
#' @return A dataframe with genes split by fgsea results.
split_genes_fgsea <- function(se, fg, groupvar = "contrast") {
  fg |> 
    filter(padj < 0.05) |> 
    group_split(term, !!sym(groupvar)) |> 
    map_dfr(function(w) {
      term <- as.character(w$term)
      gr <- as.character(w[[groupvar]])
      genes <- w$leading_edge[[1]]
      se |> 
        filter(id %in% genes & !!sym(groupvar) == gr) |> 
        mutate(term_id = term, .before = "id")
    })
}

#' Generate GSEA DE Results
#'
#' This function generates GSEA DE results based on provided GSEA and DE data.
#'
#' @param gse A list containing GSEA results for various ontologies.
#' @param de A dataframe containing DE results.
#' @param fdr_limit FDR threshold for selecting significant results. Default is 0.01.
#'
#' @return A dataframe with GSEA DE results.
gsea_de <- function(gse, de, fdr_limit = 0.01) {
  contrasts <- unique(de$contrast) |> 
    as.character()
  ontologies <- unique(gse$ontology) |> 
    as.character()
  
  map_dfr(ontologies, function(ont) {
    map_dfr(contrasts, function(ctr) {
      de_sel <- de |> 
        filter(contrast == ctr)
      gse_sel <-  gse |> 
        filter(contrast == ctr & ontology == ont)
      sig_genes <- de_sel |> 
        filter(FDR < fdr_limit) |> 
        pull(id)
      
      gse_sel |> 
        filter(fdr < fdr_limit) |> 
        unnest(leading_edge) |> 
        filter(leading_edge %in% sig_genes) |> 
        rename(id = leading_edge) |> 
        left_join(de_sel, by = "id") |> 
        select(term_id, term_name, nes, fdr, id, gene_symbol, logFC, logCPM, PValue, FDR) |> 
        add_column(ontology = ont, .before = 1) |> 
        add_column(contrast = ctr)
    })
  })
}

plot_gse_term_gene <- function(gde, fdr_limit = 0.01, logfc_limit = 0) {
  d <- gde |> 
    filter(fdr < fdr_limit & abs(logFC) > logfc_limit) |> 
    group_by(gene_symbol) |> 
    mutate(n_terms = n()) |> 
    ungroup() |> 
    group_by(term_name) |> 
    mutate(n_genes = n()) |> 
    ungroup() |> 
    mutate(
      gene_symbol = fct_reorder(gene_symbol, n_terms),
      term_name = fct_reorder(term_name, desc(n_genes))
    )
  mx <- max(abs(d$logFC))
  d |> 
    ggplot(aes(x = term_name, y = gene_symbol, fill = logFC)) +
    th +
    geom_tile() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
      scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-mx, mx))
}

#' Get Term IDs Matching Query
#'
#' This function retrieves term IDs matching a given query from GSEA results.
#'
#' @param gso A dataframe containing GSEA results.
#' @param query A string query to search for in term names.
#' @param fdr_limit FDR threshold for selecting significant terms. Default is 0.05.
#'
#' @return A vector of term IDs matching the query.
get_terms_str <- function(gso, query, fdr_limit = 0.05) {
  gso |> 
    filter(str_detect(term_name, query) & padj < fdr_limit) |>
    pull(term_id)
}

#' Select GSEA Results by Text
#'
#' This function selects GSEA results matching given text queries.
#'
#' @param gse A tibble containing GSEA results for various ontologies.
#' @param texts A vector of text queries to search for in term names.
#' @param fdr_limit FDR threshold for selecting significant terms. Default is 0.05.
#'
#' @return A dataframe with GSEA results matching the text queries.
gsa_text_selection <- function(gse, texts, fdr_limit = 0.05) {
  map(texts, function(txt) {
    gse |> 
      filter(str_detect(term_name, regex(txt, ignore_case = TRUE)) & fdr < fdr_limit) |> 
      add_column(search = txt)
  }) |> 
    list_rbind()
}

#' Create GSEA Example
#'
#' This function creates a GSEA example based on DE results, terms, and GSEA data.
#'
#' @param de A dataframe containing DE results.
#' @param terms A list of feature terms for various ontologies.
#' @param gse A list containing GSEA results for various ontologies.
#' @param example A list containing example information with ontology, term ID, and contrast.
#'
#' @return A list with details about the GSEA example.
make_gse_example <- function(de, terms, gse, example) {
  genes_used <- de |>
    pull(id) |>
    unique()
  n_genes_used <- length(genes_used)
  n_genes_in_term <- terms[[example$ontology]]$mapping |>
    filter(term_id == example$term_id) |>
    filter(id %in% genes_used) |> 
    nrow()
  n_leading_edge <- gse |> 
    filter(ontology == example$ontology) |> 
    filter(term_id == example$term_id & contrast == example$contrast) |>
    pull(leading_edge) |>
    unlist() |>
    length()
  term_name <- terms[[example$ontology]]$terms |> 
    filter(term_id == example$term_id) |> 
    pull(term_name)
  list(
    term_id = example$term_id,
    term_name = term_name,
    contrast = example$contrast,
    n_genes_used = n_genes_used,
    n_genes_in_term = n_genes_in_term,
    n_leading_edge = n_leading_edge
  )
}

#' Write GSEA Results to CSV
#'
#' This function writes the top GSEA results to a CSV file.
#'
#' @param gse A list containing GSEA results for various ontologies.
#' @param gene_info A dataframe containing gene information.
#' @param file The output CSV file path.
#' @param n_top Number of top leading edge genes to include. Default is 10.
#' @param fdr_limit FDR threshold for selecting significant terms. Default is 0.01.
#'
#' @return None. Writes the results to a CSV file.
write_gse <- function(gse, gene_info, file, n_top = 10, fdr_limit = 0.01) {
  gene_info <- gene_info |> 
    select(id, gene_symbol) |> 
    distinct()

  list_rbind(gse, names_to = "ontology") |> 
    filter(fdr < fdr_limit) |> 
    rename(id = leading_edge) |> 
    unnest(id) |> 
    left_join(gene_info, by = join_by(id)) |> 
    select(-id) |> 
    group_by(across(-gene_symbol)) |> 
    summarise(leading_edge_top_10 = str_c(head(gene_symbol, n_top), collapse = ", ")) |> 
    arrange(nes) |> 
    mutate(across(where(is.numeric), ~signif(.x, 4))) |> 
    write_csv(file)
}
