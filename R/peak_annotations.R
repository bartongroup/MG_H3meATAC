annotate_peaks_chipseeker <- function(bed) {
  gr <- GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
  GenomeInfoDb::seqlevels(gr) <- paste0("chr",  GenomeInfoDb::seqlevels(gr))

  pk <- ChIPseeker::annotatePeak(
    gr,
    tssRegion = c(-1000, 0),
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb = "org.Hs.eg.db",
    level = "gene"
  ) |> 
    as_tibble() |> 
    clean_names() |> 
    rename(
      chr = seqnames,
      entrez_id = gene_id,
      gene_symbol = symbol,
      gene_id = ensembl,
      gene_name = genename,
      annot_id = annotation
    ) |> 
    mutate(
      chr = str_remove(chr, "chr"),
      #annotation = annot_id |> str_remove("(?<=^(Intron|Exon))\\s*\\([^)]*\\)"),
      annotation = annot_id |> str_remove("\\s\\(.+\\)"),
      annotation = fct_relevel(
        annotation, 
        c(
          "Promoter",
          "5' UTR",
          "Exon",
          "Intron",
          "3' UTR",
          "Downstream",
          "Distal Intergenic"
        )
      )
    ) |> 
    relocate(peak_id, .before = 1)
  
  # Extract gene biotype from org.Hs.eg.db and add to the result
  gids <- pk$entrez_id
  gene_info <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gids,
    columns = "GENETYPE",
    keytype = "ENTREZID"
  )

  pk |> 
    mutate(gene_biotype = gene_info$GENETYPE)
}


#######################################################

annotate_peaks_annotatr <- function(bed, annotations) {
  gr <- GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)

  an <- annotatr::annotate_regions(gr, annotations) |> 
    as_tibble() |> 
    clean_names()  |> 
    rename(
      chr = seqnames,
      gene_chr = annot_seqnames,
      entrez_id = annot_gene_id,
      transcript_id = annot_tx_id,
      gene_symbol = annot_symbol,
      annotation = annot_type
    ) |> 
    mutate(
      annotation = annotation |> str_remove("hg38_genes_") |> str_remove("s$"),
      annotation = recode_values(
        annotation,
        "1to5kb" ~ "1-5 kb",
        "promoter" ~ "Promoter",
        "5UTR" ~ "5' UTR",
        "exon" ~ "Exon",
        "intron" ~ "Intron",
        "3UTR" ~ "3' UTR"
      ),
      annotation = fct_relevel(annotation, "1-5 kb", "Promoter", "5' UTR", "Exon", "Intron", "3' UTR"),
    ) |> 
    relocate(peak_id, .before = 1)

  # Extract gene biotype from org.Hs.eg.db and add to the result
  gids <- an$entrez_id |> 
    unique()
  gene_info <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = gids,
    columns = "GENETYPE",
    keytype = "ENTREZID"
  ) |> 
    as_tibble() |> 
    rename(
      entrez_id = ENTREZID,
      gene_biotype = GENETYPE
    )
 
  an |> 
    left_join(gene_info, by = c("entrez_id"))
}

# Annotatr results are very redundant, with the same
# feature (e.g. exon) appearing multiple times, from
# multiple transcripts. This function collapses redundancy.
# Primarly for counting annotations.
# Only canonical transcript are used.
reduce_annotatr <- function(an, trans) {
  tr <- trans |> 
    filter(transcript_is_canonical == 1) |> 
    select(transcript_id)

  an |> 
    mutate(transcript_id = str_remove(transcript_id, "\\..*$")) |> 
    semi_join(tr, by = "transcript_id")
}


annotate_peaks_rmsk <- function(bed, rmsk) {
  if(!("peak_id" %in% colnames(bed)))
    bed <- bed |> mutate(peak_id = str_glue("{chr}:{start}-{end}"))
  
  gr <- GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
  reps <- GenomicRanges::makeGRangesFromDataFrame(rmsk, keep.extra.columns = TRUE)

  an <- annotatr::annotate_regions(gr, reps) |> 
    as_tibble() |> 
    clean_names()  |> 
    rename(
      chr = seqnames,
      rep_id = annot_rep_id,
      rep_chr = annot_seqnames,
      rep_start = annot_start,
      rep_end = annot_end,
      rep_width = annot_width,
      rep_strand = annot_strand,
      rep_name = annot_rep_name,
      rep_class = annot_rep_class,
      rep_family = annot_rep_family,
      rep_fam = annot_rep_fam,
      rep_category = annot_rep_category
    ) |> 
    relocate(peak_id, .before = 1) |> 
    relocate(rep_id, .before = rep_chr)
}

######################################################################

build_human_annotations <- function() {
  annotations <- annotatr::build_annotations(genome = 'hg38', annotations = 'hg38_basicgenes')
  GenomeInfoDb::seqlevels(annotations) <- GenomeInfoDb::seqlevels(annotations) |> str_remove("^chr")
  return(annotations)
}

######################################################################

# annotatr leaves a lot of transcripts without gene names
annotate_missing_genes <- function(annot, transcripts) {
  an <- annot |> 
    mutate(row_id = row_number()) |> 
    mutate(transcript_id = str_remove(transcript_id, "\\..*$")) |> 
    select(-gene_biotype)

  good_data <- an |> filter(!is.na(gene_symbol))
  bad_data <- an |> filter(is.na(gene_symbol))

  added <- bad_data |>   
    select(-c(gene_symbol)) |> 
    left_join(transcripts |> select(transcript_id, gene_symbol), by = "transcript_id")

  bind_rows(good_data, added) |> 
    left_join(transcripts |> select(transcript_id, gene_biotype), by = "transcript_id") |> 
    arrange(row_id) |> 
    select(-row_id)
}

######################################################################

combine_genes_annotatr <- function(annot) {
  annot |> 
    select(peak_id:cell_line, annotation, entrez_id, gene_symbol, gene_biotype) |> 
    distinct() |> 
    filter(annotation != "promoter")
}

######################################################################

combine_hits_rbp <- function(annot) {
  annot |> 
    group_split(rbp_symbol, cell_line) |> 
    map(function(an) {
      ns <- an |> 
        count(id) |> 
        arrange(desc(n))
      lst <- map(ns$id, function(selid) {
        an |> 
          filter(id == selid) |> 
          GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
      })
      Reduce(function(x, y) IRanges::subsetByOverlaps(x, y), lst) |> 
        as_tibble()
    }) |> 
    list_rbind() |> 
    rename(chr = seqnames) |> 
    relocate(peak_id, .before = 1)
}


######################################################################

get_repeat_masker <- function() {
  asm <- "hg38"

  file_path <- str_glue("https://hgdownload.soe.ucsc.edu/goldenPath/{asm}/database/rmsk.txt.gz")
  cols <- c(
    "bin","swScore","milliDiv","milliDel","milliIns",
    "genoName","genoStart","genoEnd","genoLeft",
    "strand","repName","repClass","repFamily",
    "repStart","repEnd","repLeft","id"
  )
  r <- read_tsv(file_path, col_names = cols, show_col_types = FALSE)

  r |> select(chr = genoName, start = genoStart, end = genoEnd, strand,
           rep_name = repName, rep_class = repClass, rep_family = repFamily) |> 
    mutate(
      rep_id = sprintf("R%07d", row_number()),
      chr = str_remove(chr, "chr"),
      across(c(start, end), as.integer)
    ) |> 
    unite(rep_fam, c(rep_class, rep_family), sep = ":", remove = FALSE) |> 
    mutate(
      rep_category = str_remove(rep_class, "\\?"),
      rep_category = if_else(
        rep_category %in% c("Unknown", "snRNA", "RC", "tRNA", "rRNA", "srpRNA", "scRNA", "RNA"),
        "Other",
        rep_category
      )
    )
}


