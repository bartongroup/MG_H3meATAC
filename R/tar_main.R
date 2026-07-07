targets_main <- function() {

  get_annotations <- tar_plan(
    config_annot = list(ensembl = config$ensembl, species = config$species),
    mart = new_mart(config_annot),
    chromosomes = biomart_fetch_chromosomes(mart),
    #structures = biomart_fetch_structures(mart, gene_ids),
    gns = biomart_fetch_genes(mart),
    trans = biomart_fetch_transcripts(mart),
    terms = get_functional_terms(config_annot, gns),
    fterms = prepare_terms_fenr(terms),
    rmsk = get_repeat_masker()
  )

  setup_experiment <- tar_plan(
    tar_target(config_file, "config.json", format = "file"),
    config = jsonlite::read_json(config_file),
    tar_target(params_file, "params.json", format = "file"),
    params = jsonlite::read_json(params_file),
    tar_target(meta_file, config$metadata_file, format = "file"),
    metadata = read_metadata(meta_file)
  )

  qc <- tar_plan(
    fastq_len = read_fastq_lengths(config, metadata),

    fscrn = parse_fscreens(config, metadata),
    qcs = parse_qcs(config, metadata, paired = TRUE),
    idxstats = parse_idxstats(config, metadata),
    fig_chrom_reads = plot_chrom_proportion(idxstats, c(1:22, "X", "Y", "MT"), per_length = FALSE),
    fig_mitochondrial_prop = plot_mitochondrial_proportion(idxstats, metadata),
    
    fig_fscreen = plot_fscreen_map(fscrn),
    fig_read_qual = plot_qualities(qcs),
    fig_read_qual_clust = plot_cluster_qualities(qcs),

    bowtie_sum_human = read_bowtie_summaries(config, metadata),
    bowtie_sum_dm = read_bowtie_summaries(config, metadata, bam_dir = "bam_dm"),
    bowtie_sum_merged = merge_two_bowtie_summaries(bowtie_sum_human, bowtie_sum_dm),
    fig_bowtie_sum_human = plot_bowtie_summaries(bowtie_sum_human, metadata),
    fig_bowtie_sum_dm = plot_bowtie_summaries(bowtie_sum_dm, metadata),
    fig_spikein_ratios = plot_spikein_ratios(bowtie_sum_merged, metadata),

    picard_dup = read_picard_duplication(config, metadata),
    fig_picard_dup = plot_picard_duplication(picard_dup, metadata),

    spikein_norms = read_spikein_norms(config, metadata),
    fig_spikein_norms = plot_spikein_norms(spikein_norms, metadata) 
  )

  read_data <- tar_plan(
    chrom_sizes = readr::read_tsv(file.path(config$data_path, config$dir$genome, "chrom_sizes.txt"), col_names = c("chr", "chrom_size"), show_col_types = FALSE),
    annotations = build_human_annotations()
  )

  map_peaks <- tar_map(
    values = tibble::tibble(TYPE = c("narrow")),
    names = TYPE,

    tar_plan(
      peaks = read_peak_files(config, metadata, format = TYPE),
      fig_npeaks = plot_number_of_peaks(peaks, metadata),
      fig_peak_width = plot_peak_width_distribution(peaks),
      fig_peak_coverage = plot_peak_coverage(peaks, chrom_sizes, metadata),

      # Annotate all peaks
      peaks_annotatr_raw = annotate_peaks_annotatr(peaks, annotations),
      peaks_annotatr = reduce_annotatr(peaks_annotatr_raw, trans),
      peaks_chipseeker = annotate_peaks_chipseeker(peaks),
      peaks_chipseeker_red = peaks_chipseeker |> dplyr::filter(!(annotation %in% c("Downstream", "Distal Intergenic"))),
      peaks_rmsk = annotate_peaks_rmsk(peaks, rmsk),

      # Random peaks
      prand = make_random_regions(peaks, chrom_sizes, n_regions = 100000),
      prand_annotatr = annotate_peaks_annotatr(prand, annotations) |> reduce_annotatr(trans),
      prand_chipseeker = annotate_peaks_chipseeker(prand),
      prand_chipseeker_red = prand_chipseeker |> dplyr::filter(!(annotation %in% c("Downstream", "Distal Intergenic"))),
      prand_rmsk = annotate_peaks_rmsk(prand, rmsk),

      fig_annotatr = plot_peak_annotatios(peaks_annotatr),
      fig_chipseeker = plot_peak_annotatios(peaks_chipseeker),

      # These are proper figures
      fig_groups_annotatr = plot_group_peak_annotations(peaks_annotatr, rand_an = prand_annotatr, palette = 2),
      fig_groups_chipseeker = plot_group_peak_annotations(peaks_chipseeker, rand_an = prand_chipseeker),
      fig_groups_chipseeker_red = plot_group_peak_annotations(peaks_chipseeker_red, rand_an = prand_chipseeker_red),
      fig_groups_rmsk = plot_group_peak_annotations(peaks_rmsk, annot_var = "rep_category", rand_an = prand_rmsk, palette = 3),
      
      fig_groups_tss = plot_peaks_tss(peaks_chipseeker, rand_an = prand_chipseeker, max_dist = 3000, ncol = 5),

      # Prepare for DiffBind
      meta_db = make_diffbind_metadata(metadata, config, type = TYPE),
      sav_meta_db = make_diffbind_metadata_for_cluster(metadata, config, type = TYPE),
      # These steps are very slow and require a lot of memory, so we run it on
      # the cluster and cache the result. After creating
      # "diffbind_metadata_narrow.csv" copy it to the cluster and run the script
      # "diffbind_count_cluster.R". Copy *.rds files back here.
      # db = read_and_count_diffbind(meta_db, min_overlap = 2),
      # db_spike = spikein_normalization(dbc),
    )
  )

  # Map over normalization and peak type
  map_diffbind <- tar_map(
    values = tidyr::expand_grid(
      NORM = c("libs"),
      TYPE = c("narrow")
    ) |> dplyr::mutate(
      NAME = stringr::str_glue("{NORM}_{TYPE}")
    ),
    names = NAME,

    tar_plan(
      db = read_cache(TYPE, NORM),

      fig_clust = plot_db_clustering(db, metadata),
      fig_peaks_chr11 = plot_db_biding_region(db, metadata, make_region("11", 40e6, 60e6)),
    ),

    # Map over two models
    tar_map(
      values = tibble::tibble(MODEL = c("full", "groups")),
      names = MODEL,

      tar_plan(
        # Perform DiffBind analysis
        dbc = de_model(db, MODEL),
        # Extract DE results in a useful format
        dec = extract_de_results(dbc),

        fig_volcano = plot_volcano(dec, fdr_limit = params$fdr_limit, logfc_limit = params$logfc_limit),
        fig_pdist = plot_pdist(dec),
        fig_up_down = plot_up_down(dec, fdr_limit = params$fdr_limit, logfc_limit = params$logfc_limit),

        # Annotate DE results
        dec_annotatr = annotate_peaks_annotatr(dec, annotations) |> reduce_annotatr(trans),
        dec_chipseeker = annotate_peaks_chipseeker(dec),
        dec_rmsk = annotate_peaks_rmsk(dec, rmsk),

        dec_annotatr_sig = dec_annotatr |> dplyr::filter(fdr < params$fdr_limit & abs(log_fc) > params$logfc_limit),
        dec_chipseeker_sig = dec_chipseeker |> dplyr::filter(fdr < params$fdr_limit & abs(log_fc) > params$logfc_limit),
        dec_rmsk_sig = dec_rmsk |> dplyr::filter(fdr < params$fdr_limit & abs(log_fc) > params$logfc_limit),

        fig_dec_annotatr = plot_group_peak_annotations(dec_annotatr_sig, group_var = "contrast_direction", palette = 2, n_top = NULL),
        fig_dec_chipseeker = plot_group_peak_annotations(dec_chipseeker_sig, group_var = "contrast_direction", n_top = NULL),
        fig_dec_chipseeker_err = plot_group_peak_annotations(dec_chipseeker_sig, group_var = "contrast_direction", n_top = NULL, stacked = FALSE),
        fig_dec_rmsk_da = plot_peak_rmsk(dec_rmsk_sig, fdr_limit = params$fdr_limit),
        fig_dec_rmsk_all = plot_peak_rmsk(dec_rmsk),

        fig_dec_rmsk = plot_group_peak_annotations(dec_rmsk_sig, group_var = "contrast_direction", annot_var = "rep_category", n_top = NULL, palette = 3),

        # Test statistical sigificance of proportion results
        tab_dec_chipseeker = tabulate_annotations(dec_chipseeker_sig, group_var = "contrast_direction", annot_var = "annotation"),
        test_dec_chipseeker = annot_group_test(tab_dec_chipseeker),
        fig_test_dec_chipseeker = plot_annot_group_test(test_dec_chipseeker),
        tab_dec_rmsk = tabulate_annotations(dec_rmsk_sig, group_var = "contrast_direction", annot_var = "rep_category"),
        test_dec_rmsk = annot_group_test(tab_dec_rmsk),
        fig_test_dec_rmsk = plot_annot_group_test(test_dec_rmsk),

        sav_dec = write_dec(dec_chipseeker, paste0("diffbind_chipseeker_", MODEL, "_", NAME)),

        deg = prepare_de_genes(dec_chipseeker),
        fig_deg_biotype = plot_gene_biotype_distribution(deg, fdr_limit = params$fdr_limit),
        gse = fgsea_all_terms(deg, fterms, rank_expr = "-log_fc * log10(p_value)"),
        gse_tab = tabulate_gsea(gse, gns)
      )
    )

  )

  # This is done only for ChIPseeker, differential binding in pair-wise comparison significant genes,
  # lib-size normalisation, narrow peaks
  search_for_motifs <- tar_plan(
    dec_for_motifs = dec_chipseeker_sig_groups_libs_narrow,

    tar_target(genome_file, file.path(config$data_path, config$dir$genome, "Homo_sapiens.canonical.fa"), format = "file"),
    db_ranges = dec_for_motifs |> dplyr::distinct(peak_id, chr, start, end, strand),
    db_seqs = fetch_range_sequences(genome_file, db_ranges),
    tf_motifs = get_jaspar_matrix(),
    n_tf_motifs = tf_motifs |> length(),
    tf_res = search_sequences_for_motifs(db_seqs, tf_motifs),
    tf_res_sig = tf_res |> dplyr::filter(motif_fdr < 0.01),
    dec_res = dec_for_motifs |> dplyr::left_join(tf_res_sig |> dplyr::mutate(peak_id = as.character(peak_id)), by = "peak_id", relationship = "many-to-many"),

    motif_stats = get_motif_stats(dec_res, dec_for_motifs),
    motif_stats_promoters = get_motif_stats(dec_res |> dplyr::filter(annotation == "Promoter"), dec_for_motifs |> dplyr::filter(annotation == "Promoter"))
  )



 

  for_report <- tar_plan(
    fig_all_biotype = plot_all_genes_biotype(),

    consensus_peaks = db_libs_narrow |> extract_consensus(),
    tbl_rmsk_peaks =  annotate_peaks_rmsk(consensus_peaks, rmsk),
    # Do it manually, cache:
    tbl_rmsk_stats = rmsk_stats(rmsk, consensus_peaks, n_perm = 50, n_top = 20),
    #tbl_rmsk_stats = readr::read_rds("cache/tbl_rmsk_stats.rds"),

    sav_rmsk_peaks = tbl_rmsk_peaks |> readr::write_csv("tab/rmask_peaks.csv"),

    tbl_annot_fisher_example = annotation_fisher_example(test_dec_chipseeker_groups_libs_narrow),

    db_none_narrow = no_normalization(db_libs_narrow),
    #fig_norm_compare = plot_all_total_bindings(db_none_narrow, db_libs_narrow, db_spike_narrow, metadata),

    fig_de_abu_example_1 = plot_peak_bindings(dbc_full_libs_narrow, "5:70822898-70823402", metadata),
    #fig_de_peak_example_1 = plot_region_ggbio_peak_id(dec_full_libs_narrow |> dplyr::filter(contrast == "48h-iPS"), "5:70822898-70823402", config, metadata, peaks_narrow, peaks_narrow, normfac |> dplyr::mutate(norm = libs_norm)),

    normfac = tibble::tibble(
      sample = db_libs_narrow$samples$SampleID,
      lib_size = db_libs_narrow$norm$DESeq2$lib.sizes,
      libs_norm = db_libs_narrow$norm$DESeq2$norm.facs
    ),
    #fig_compare_spikein_norm = plot_compare_spikein_normalisation(normfac, spikein_norms),
    png_peaks_narrow_ctrl = plot_peak_distribution_sample(peaks_narrow, "iPS_Ctrl_1", chrom_sizes, bar_size = 1) |> gs("peaks_narrow_dist_ctrl", width = 8, height = 8),
    png_peaks_narrow_set = plot_peak_distribution_sample(peaks_narrow, "iPS_SET_1", chrom_sizes, bar_size = 1) |> gs("peaks_narrow_dist_set", width = 8, height = 8),
    png_peaks_narrow_48 = plot_peak_distribution_sample(peaks_narrow, "48h_SET_1", chrom_sizes, bar_size = 1) |> gs("peaks_narrow_dist_48", width = 8, height = 8),
    png_peaks_narrow_11 = plot_peak_distribution_chrom(peaks_narrow, "11", metadata) |> gs("peaks_narrow_11", width = 8, height = 8)
  )

  for_manuscript <- tar_plan(
    mn_fdr_limit = 0.01,
    mn_logfc_limit = 0.58,
    mn_ctr = c("48h_Ctrl-iPS_Ctrl", "48h_SETDB1-48h_Ctrl", "iPS_SETDB1-iPS_Ctrl"),

    png_volcano = mn_plot_volcano(dec_chipseeker_groups_libs_narrow |> filter(contrast %in% mn_ctr), fdr_limit = mn_fdr_limit, logfc_limit = mn_logfc_limit) |> gp("volcano", width = 7, height = 3),

    png_groups_chipseeker = plot_group_peak_annotations(peaks_chipseeker_narrow, rand_an = prand_chipseeker_narrow) |> 
      gp("chipseeker_annotations", width = 4, height = 4),
    png_chipseeker_annotations = mn_plot_group_peak_annotations(dec_chipseeker_sig_groups_libs_narrow, group_var = "contrast_direction", n_top = NULL) |>
      gp("chipseeker_da", width = 4, height = 4),
    png_chipseeker_rmsk = plot_group_peak_annotations(dec_rmsk_sig_groups_libs_narrow, group_var = "contrast_direction", annot_var = "rep_category", n_top = NULL, palette = 3) |>
      gp("rmsk_da", width = 4, height = 4)
    
  )

  for_shiny <- tar_plan(
    #sav_shiny = prepare_for_shiny(dbc_full_libs_narrow, dbc_full_spike_narrow, dbc_groups_libs_narrow, dbc_groups_spike_narrow,
    #                              metadata, terms, gns)
  )


  get_environment <- tar_plan(
    tar_target(version_file, file.path(config$data_path, "environment.txt"), format = "file"),
    software_versions = get_software_versions(
      version_file,
      selection = c("fastqc", "fastq-screen", "fastp", "ribodetector", "multiqc", "bowtie2", "samtools", "bedtools", "snakemake", "drmaa", "picard", "macs3")
    )
  )


  c(
    get_environment,
    get_annotations,
    setup_experiment,
    qc,
    read_data,
    map_peaks,
    map_diffbind,
    search_for_motifs,
    for_report
    #for_manuscript,
    #for_shiny
  )
  
}