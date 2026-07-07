#' Create a new Ensembl biomaRt object
#'
#' This function initializes a new Ensembl biomaRt object using the specified
#' configuration.
#'
#' @param config A list containing Ensembl configuration details, including
#'   `biomart`, `dataset`, and `version`.
#'
#' @return An object of class `Mart` from the biomaRt package.
new_mart <- function(config) {
  biomaRt::useEnsembl(
    biomart = config$ensembl$biomart,
    dataset = config$ensembl$dataset,
    version = config$ensembl$version
  )
}


#' Fetch Chromosome Names from Ensembl biomaRt
#'
#' This function retrieves "proper" chromosome names (excluding patches and
#' haplotypes) from an Ensembl biomaRt object.
#'
#' @param mart An object of class `Mart` from the biomaRt package.
#'
#' @return A character vector of chromosome names.
biomart_fetch_chromosomes <- function(mart) {
  chrnames <-  biomaRt::getBM(attributes = "chromosome_name", mart = mart, uniqueRows = TRUE)
  nums <- chrnames |>
    tibble::as_tibble() |>
    dplyr::filter(str_detect(chromosome_name, "^\\d{1,2}$")) |>
    dplyr::pull(chromosome_name) |>
    as.numeric() |>
    sort() |>
    as.character()
  c(nums, "X", "Y", "MT")
}

#' Fetch Gene Information from Ensembl biomaRt
#'
#' This function retrieves gene information from an Ensembl biomaRt object, with
#' options to include NCBI gene IDs.
#'
#' @param mart An object of class `Mart` from the biomaRt package.
#' @param use_cache Logical indicating whether to use cached data. Default is
#'   `TRUE`.
#' @param with_ncbi Logical indicating whether to include NCBI gene IDs. Default
#'   is `FALSE`.
#'
#' @return A tibble containing gene information.
biomart_fetch_genes <- function(mart, use_cache = TRUE, with_ncbi = FALSE) {
  attr <- c(
    chr = "chromosome_name",
    start = "start_position",
    tss = "transcription_start_site",
    end = "end_position",
    strand = "strand",
    gene_biotype = "gene_biotype",
    gc_content = "percentage_gene_gc_content",
    id = "ensembl_gene_id",
    gene_symbol = "external_gene_name",
    description = "description"
  )
  if(with_ncbi)
    attr <- c(attr, c(ncbi_id = "entrezgene_id"))

  biomaRt::getBM(attributes = attr, mart = mart, useCache = use_cache) |>
  dplyr::rename(!!attr) |>
    dplyr::mutate(
      description = str_remove(description, "\\s\\[.*\\]"),
      gene_symbol = na_if(gene_symbol, "")
    ) |>
    tibble::as_tibble() |>
    dplyr::mutate(gene_symbol = if_else(is.na(gene_symbol), id, gene_symbol))
}

biomart_fetch_transcripts <- function(mart, use_cache = TRUE) {
  attr <- c(
    gene_id = "ensembl_gene_id",
    transcript_id = "ensembl_transcript_id",
    transcript_is_canonical = "transcript_is_canonical",
    gene_symbol = "external_gene_name"
  )

  biomaRt::getBM(attributes = attr, mart = mart, useCache = use_cache) |>
    dplyr::rename(!!attr) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      dplyr::across(
        c(gene_symbol),
        ~ dplyr::na_if(.x, "")
      )
    )
}


#' Fetch Gene Structures from Ensembl biomaRt
#'
#' This function retrieves gene structures, including exons and UTRs, from an
#' Ensembl biomaRt object.
#'
#' @param mart An object of class `Mart` from the biomaRt package.
#' @param gene_ids A character vector of Ensembl gene IDs to retrieve structures
#'   for.
#' @param use_cache Logical indicating whether to use cached data. Default is
#'   `TRUE`.
#'
#' @return A tibble containing gene structure information.
biomart_fetch_structures <- function(mart, gene_ids, use_cache = TRUE) {
  attr <-  c(
    chr = "chromosome_name",
    id = "ensembl_gene_id",
    transcript_id = "ensembl_transcript_id",
    exon_id = "ensembl_exon_id",
    gene_symbol = "external_gene_name",
    start = "start_position",
    end = "end_position",
    tss = "transcription_start_site",
    exon_start = "exon_chrom_start",
    exon_end = "exon_chrom_end",
    start_5_utr = "5_utr_start",
    end_5_utr = "5_utr_end",
    start_3_utr = "3_utr_start",
    end_3_utr = "3_utr_end",
    transcript_start = "transcript_start",
    transcript_end = "transcript_end",
    strand = "strand"
  )

  res <- biomaRt::getBM(attributes = attr, mart = mart, filters = "ensembl_gene_id", values = gene_ids, useCache = use_cache) |>
    dplyr::rename(!!attr) |>
    dplyr::mutate(
      gene_symbol = na_if(gene_symbol, "")
    ) |>
    tibble::as_tibble() |>
    dplyr::mutate(gene_symbol = if_else(is.na(gene_symbol), id, gene_symbol))

  # transcript_is_canonical is on a "feature_page" page, unlike all other attributes above, which
  # are on the "structure" page. Need a separate grab.
  attr <- c(
    transcript_id = "ensembl_transcript_id",
    canonical = "transcript_is_canonical"
  )
  can <- biomaRt::getBM(attributes = attr, mart = mart, filters = "ensembl_gene_id", values = gene_ids, useCache = use_cache) |>
    dplyr::rename(!!attr) |> 
    tibble::as_tibble() |>
    dplyr::mutate(canonical = tidyr::replace_na(canonical, 0) |> as.logical())
  
  res |> 
    dplyr::left_join(can, by = join_by(transcript_id))
}

#' Fetch Protein Information from Ensembl biomaRt
#'
#' This function retrieves protein information, including Ensembl peptide IDs,
#' from an Ensembl biomaRt object.
#'
#' @param mart An object of class `Mart` from the biomaRt package.
#' @param gene_ids A character vector of Ensembl gene IDs to retrieve protein
#'   information for.
#' @param use_cache Logical indicating whether to use cached data. Default is
#'   `TRUE`.
#'
#' @return A tibble containing protein information.
biomart_fetch_proteins <- function(mart, gene_ids, use_cache = TRUE) {
  attr <- c(
    id = "ensembl_gene_id",
    gene_symbol = "external_gene_name",
    protein_id = "ensembl_peptide_id"
  )

  biomaRt::getBM(attributes = attr, mart = mart, filters = "ensembl_gene_id", values = gene_ids, useCache = use_cache) |>
    dplyr::rename(!!attr) |>
    tibble::as_tibble() |>
    filter(protein_id != "")
}

#' Fetch Gene Lengths and GC Content from Ensembl biomaRt
#'
#' This function retrieves gene lengths and GC content from an Ensembl biomaRt
#' object, inspired by the `getGeneLengthAndGCContent` function in EDASeq.
#'
#' @param mart An object of class `Mart` from the biomaRt package.
#' @param gene_ids A character vector of Ensembl gene IDs to retrieve lengths
#'   and GC content for.
#' @param use_cache Logical indicating whether to use cached data. Default is
#'   `TRUE`.
#'
#' @return A tibble containing gene lengths and GC content.
biomart_fetch_gene_lengths <- function(mart, gene_ids, use_cache = TRUE) {
  exons <- biomaRt::getBM(attributes = c(
    "ensembl_gene_id",
    "ensembl_exon_id",
    "chromosome_name",
    "exon_chrom_start",
    "exon_chrom_end"),
    mart = mart, filters = "ensembl_gene_id", values = gene_ids, useCache = use_cache
  )
  # find all exons, reduce overlap
  coords <- exons |>
    group_split(ensembl_gene_id) |>
    map(function(w) {
      IRanges::IRanges(w$exon_chrom_start, w$exon_chrom_end) |>
        IRanges::reduce() |>
        as.data.frame() |>
        as_tibble() |>
        mutate(ensembl_gene_id = w$ensembl_gene_id[1], chr = w$chromosome_name[1])
    }, .progress = TRUE) |>
    list_rbind() |> 
    rename(id = ensembl_gene_id, length = width)
  coords |>
    group_by(id, chr) |>
    summarise(length = sum(length))
}

#' Fetch Gene Sequences from Ensembl biomaRt
#'
#' This function retrieves gene sequences, including exons and introns, from an
#' Ensembl biomaRt object in batches.
#'
#' @param mart An object of class `Mart` from the biomaRt package.
#' @param gids A character vector of Ensembl gene IDs to retrieve sequences for.
#' @param batch_size Integer indicating the number of gene IDs to process per
#'   batch. Default is `1000`.
#'
#' @return A tibble containing gene sequences.
biomart_fetch_gene_sequences <- function(mart, gids, batch_size = 1000) {
  attr <- c(
    id = "ensembl_gene_id",
    gene_symbol = "external_gene_name",
    sequence = "gene_exon_intron"
  )

  n_batch <- ceiling(length(gids) / batch_size)
  pb <- progress_bar$new(format = "  Fetching gene sequences [:bar] :percent eta: :eta", total = n_batch)

  tibble(
    id = gids,
    batch = ntile(id, n_batch)
  ) |>
    group_split(batch) |>
    map(function(w) {
      tb <- biomaRt::getBM(attributes = attr, mart = mart, filters = "ensembl_gene_id", values = w$id) |>
        dplyr::rename(!!attr) |>
        tibble::as_tibble()
      pb$tick()
      tb
    }) |>
    list_rbind()
}



