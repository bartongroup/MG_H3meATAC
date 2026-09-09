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
    metadata = read_metadata(meta_file),
    tar_target(genome_file, file.path(config$data_path, config$dir$genome, "Homo_sapiens.canonical.fa"), format = "file"),
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
      median_width = median(peaks$end - peaks$start + 1),

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
      # If your machine has enough memory (>100 GB), uncomment the line below.
      # db = read_and_count_diffbind(meta_db, min_overlap = 2)
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

        png_volcano = plot_volcano(dec, fdr_limit = params$fdr_limit, logfc_limit = params$logfc_limit) |> gs(paste0("volcano", "_", MODEL, "_", NORM, "_", TYPE), 8, 6),
        fig_pdist = plot_pdist(dec),
        fig_up_down = plot_up_down(dec, fdr_limit = params$fdr_limit, logfc_limit = params$logfc_limit),

        # Annotate DE results
        dec_annotatr = annotate_peaks_annotatr(dec, annotations) |> reduce_annotatr(trans),
        dec_chipseeker = annotate_peaks_chipseeker(dec),
        dec_rmsk = annotate_peaks_rmsk(dec, rmsk),

        # Full peak data, annotated and not, used for testing
        dec_annotatr_all = add_unannotated_peaks(dec, dec_annotatr, annot_var = "annotation"),
        dec_chipseeker_all = add_unannotated_peaks(dec, dec_chipseeker, annot_var = "annotation"),
        dec_rmsk_all = add_unannotated_peaks(dec, dec_rmsk, annot_var = "rep_category"),
        dec_rmsk_all_ltr = add_unannotated_peaks(dec, dec_rmsk |> dplyr::filter(rep_category == "LTR"), annot_var = "rep_family"),

        dec_annotatr_sig = dec_annotatr |> dplyr::filter(fdr < params$fdr_limit & abs(log_fc) > params$logfc_limit),
        dec_annotatr_nsig = dec_annotatr |> dplyr::filter(fdr >= params$fdr_limit & abs(log_fc) <= params$logfc_limit),
        dec_chipseeker_sig = dec_chipseeker |> dplyr::filter(fdr < params$fdr_limit & abs(log_fc) > params$logfc_limit),
        dec_chipseeker_nsig = dec_chipseeker |> dplyr::filter(fdr >= params$fdr_limit & abs(log_fc) <= params$logfc_limit),
        dec_rmsk_sig = dec_rmsk |> dplyr::filter(fdr < params$fdr_limit & abs(log_fc) > params$logfc_limit),
        dec_rmsk_nsig = dec_rmsk |> dplyr::filter(fdr >= params$fdr_limit & abs(log_fc) <= params$logfc_limit),

        fig_dec_annotatr = plot_group_peak_annotations(dec_annotatr_sig, rand_an = dec_annotatr_nsig, group_var = "contrast_direction", palette = 2, n_top = NULL),
        fig_dec_chipseeker = plot_group_peak_annotations(dec_chipseeker_sig, rand_an = dec_chipseeker_nsig, group_var = "contrast_direction", n_top = NULL),
        fig_dec_chipseeker_err = plot_group_peak_annotations(dec_chipseeker_sig, rand_an = dec_chipseeker_nsig, group_var = "contrast_direction", n_top = NULL, stacked = FALSE),
        fig_dec_rmsk_da = plot_peak_rmsk(dec_rmsk_sig, fdr_limit = params$fdr_limit),
        fig_dec_rmsk_all = plot_peak_rmsk(dec_rmsk),

        fig_dec_rmsk = plot_group_peak_annotations(dec_rmsk_sig, rand_an = dec_rmsk_nsig, group_var = "contrast_direction", annot_var = "rep_category", n_top = NULL, palette = 3),

        # Test statistical sigificance of proportion results
        tab_dec_chipseeker = tabulate_annotations(dec_chipseeker_sig, group_var = "contrast_direction", annot_var = "annotation"),
        test_dec_chipseeker = annot_group_test(tab_dec_chipseeker),
        fig_test_dec_chipseeker = plot_annot_group_test(test_dec_chipseeker),
        tab_dec_rmsk = tabulate_annotations(dec_rmsk_sig, group_var = "contrast_direction", annot_var = "rep_category"),
        test_dec_rmsk = annot_group_test(tab_dec_rmsk),
        fig_test_dec_rmsk = plot_annot_group_test(test_dec_rmsk),

        # Test for DE annotation enrichment
        test_de_chipseeker = annot_de_test(dec_chipseeker_all, fdr_limit = params$fdr_limit, logfc_limit = params$logfc_limit, annot_var = "annotation"),
        test_de_annotatr = annot_de_test(dec_annotatr_all, fdr_limit = params$fdr_limit, logfc_limit = params$logfc_limit, annot_var = "annotation"),
        test_de_rmsk = annot_de_test(dec_rmsk_all, fdr_limit = params$fdr_limit, logfc_limit = params$logfc_limit, annot_var = "rep_category"),
        test_de_rmsk_ltr = annot_de_test(dec_rmsk_all_ltr, fdr_limit = params$fdr_limit, logfc_limit = params$logfc_limit, annot_var = "rep_family"),
        fig_test_de_chipseeker = mn_plot_annot_group_test(test_de_chipseeker),
        fig_test_de_annotatr = mn_plot_annot_group_test(test_de_annotatr),
        fig_test_de_rmsk = mn_plot_annot_group_test(test_de_rmsk, max_or = 2.9),
        fig_test_de_rmsk_ltr = mn_plot_annot_group_test(test_de_rmsk_ltr, max_or = 4.2),

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
  motif_analysis <- tar_plan(
    dec_for_motifs = dec_chipseeker_sig_groups_libs_narrow |> dplyr::filter(abs(log_fc) > 2),
    n_peaks_for_motifs = nrow(dec_for_motifs),
    consensus_peaks =  extract_consensus(db_libs_narrow),
    consensus_seqs = fetch_range_sequences(genome_file, consensus_peaks),

    tf_motifs = get_jaspar_matrix(),
    n_tf_motifs = tf_motifs |> length(),

    # Scan for motifs to find locations
    tf_scan = search_peaks_for_motifs(dec_for_motifs, tf_motifs, genome_file),
    tf_scan_sig = tf_scan |> dplyr::filter(motif_fdr < 0.01),
    dec_tf_scan = dec_for_motifs |> dplyr::left_join(tf_scan_sig |> dplyr::mutate(peak_id = as.character(peak_id)), by = "peak_id", relationship = "many-to-many"),
    motif_stats = get_motif_stats(dec_tf_scan, dec_for_motifs),
    motif_stats_promoters = get_motif_stats(dec_tf_scan |> dplyr::filter(annotation == "Promoter"), dec_for_motifs |> dplyr::filter(annotation == "Promoter")),

    # Motif enrichment
    tf_enr = motif_enrichment(dec_for_motifs, consensus_seqs, genome_file, tf_motifs, min_target = 10),
    tf_comb = dplyr::left_join(tf_enr, dec_tf_scan,  by = dplyr::join_by(motif_name, motif_id, contrast_direction), relationship = "many-to-many"),

    # Motif enrichment for selected contrasts, deeper
    dec_for_motifs_1 = dec_chipseeker_sig_groups_libs_narrow |> dplyr::filter(abs(log_fc) > 1 & contrast_direction %in% c("48h_SETDB1-48h_Ctrl_up", "iPS_SETDB1-iPS_Ctrl_up")),
    tf_enr_1 = motif_enrichment(dec_for_motifs_1, consensus_seqs, genome_file, tf_motifs, min_target = 10),

    fig_enr_contrasts_logp_down = plot_enr_contrasts_logp(tf_enr, ctr1 = "48h_Ctrl-iPS_Ctrl_down", ctr2 = "48h_SETDB1-iPS_SETDB1_down"),
    fig_enr_contrasts_logp_up = plot_enr_contrasts_logp(tf_enr, ctr1 = "48h_Ctrl-iPS_Ctrl_up", ctr2 = "48h_SETDB1-iPS_SETDB1_up")
  )


  cuttag_comparison <- tar_plan(
    peaks_cut = readr::read_rds("cache/peaks_cut.rds"),

    dec_cut = readr::read_rds("cache/dec_cut.rds") |> dplyr::mutate(sig = fdr < 0.01 & abs(log_fc) > 0.58),
    dec_atac = dec_chipseeker_groups_libs_narrow |> dplyr::mutate(sig = fdr < 0.01 & abs(log_fc) > 1),
    cnt_over_cut_atac = count_cut_atac_overlaps(dec_atac, dec_cut),
    cmp_cut_atac = compare_cut_atac(dec_atac, dec_cut),
    over_cut_atac_ctrl = get_cut_atac_overlaps(dec_atac, dec_cut, "48h_Ctrl-iPS_Ctrl", "48h_Ctrl-iPS_Ctrl"),
    over_cut_atac_set = get_cut_atac_overlaps(dec_atac, dec_cut, "48h_SETDB1-iPS_SETDB1", "48h_SETDB1-iPS_SETDB1"),
    fig_over_cut_atac_fcfc_ctrl = plot_cut_atac_fcfc(over_cut_atac_ctrl),
    fig_over_cut_atac_fcfc_set = plot_cut_atac_fcfc(over_cut_atac_set),

    png_peaks_cut_atac_tss = plot_peaks_cut_atac_tss(peaks_cut, peaks_chipseeker_narrow, n_top = 300) |> gs("cut_atac_peaks", 8, 12)
  )
 

  for_report <- tar_plan(
    fig_pca = plot_pca(db_libs_narrow, metadata),
    fig_pca_ips = plot_pca(db_libs_narrow, metadata |> filter(time_point == 'iPS')),
    fig_pca_48 = plot_pca(db_libs_narrow, metadata |> filter(time_point == '48h')),
    fig_all_biotype = plot_all_genes_biotype(),

    tbl_rmsk_peaks =  annotate_peaks_rmsk(consensus_peaks, rmsk),
    tbl_rmsk_stats = rmsk_stats(rmsk, consensus_peaks, group_var = "rep_fam", n_perm = 100, n_top = 20),
    #tbl_rmsk_stats = readr::read_rds("cache/tbl_rmsk_stats.rds"),

    sav_rmsk_peaks = tbl_rmsk_peaks |> readr::write_csv("tab/rmask_peaks.csv"),

    tbl_annot_fisher_example = annotation_fisher_example(test_dec_chipseeker_groups_libs_narrow),

    db_none_narrow = no_normalization(db_libs_narrow),
    #fig_norm_compare = plot_all_total_bindings(db_none_narrow, db_libs_narrow, db_spike_narrow, metadata),

    fig_peak_example_1 = plot_region_ggbio(metadata$sample, config, metadata, make_region(7, 143620700, 143621600), peaks_broad = NULL, peaks_narrow = peaks_narrow),
    fig_de_abu_example_1 = plot_peak_bindings(dbc_groups_libs_narrow, "8:37374556-37375887", metadata),
    dec_ctr = dec_groups_libs_narrow |> dplyr::filter(contrast == "48h_Ctrl-iPS_Ctrl"),
    fig_de_peak_example_1 = plot_region_ggbio_peak_id(dec_ctr, "8:37374556-37375887", config, metadata, peaks_narrow = peaks_narrow, normfac = normfac |> dplyr::mutate(norm = libs_norm)),

    smpls = metadata |> filter(treatment == "Ctrl" ) |> pull(sample) |> as.character(),
    png_cut_atac_example = plot_region_ggbio(smpls, config, metadata, make_region(18, 69935096, 69935793), peaks_broad = NULL, peaks_narrow = peaks_narrow) |> gs("cut_atac_example", width = 8, height = 8),

    normfac = tibble::tibble(
      sample = db_libs_narrow$samples$SampleID,
      lib_size = db_libs_narrow$norm$DESeq2$lib.sizes,
      libs_norm = db_libs_narrow$norm$DESeq2$norm.facs
    ),
    #fig_compare_spikein_norm = plot_compare_spikein_normalisation(normfac, spikein_norms),
    png_peaks_narrow_ctrl = plot_peak_distribution_sample(peaks_narrow, "iPS_Ctrl_1", chrom_sizes, bar_size = 1) |> gs("peaks_narrow_dist_ctrl", width = 8, height = 8),
    png_peaks_narrow_set = plot_peak_distribution_sample(peaks_narrow, "iPS_SET_1", chrom_sizes, bar_size = 1) |> gs("peaks_narrow_dist_set", width = 8, height = 8),
    png_peaks_narrow_48 = plot_peak_distribution_sample(peaks_narrow, "48h_SET_1", chrom_sizes, bar_size = 1) |> gs("peaks_narrow_dist_48", width = 8, height = 8),
    png_peaks_narrow_11 = plot_peak_distribution_chrom(peaks_narrow, "11", metadata) |> gs("peaks_narrow_11", width = 8, height = 8),

    sav_bed_consensus = write_bed(consensus_peaks, "tab/consensus_peaks.bed"),
    sav_bed_de = write_bed(dec_groups_libs_narrow, "tab/de_peaks.bed", extra_cols = c("contrast_direction", "log_fc")),

    sav_over_cut_atac_ctrl = over_cut_atac_ctrl |> readr::write_csv("tab/compare_cut_atac_ctrl.csv"),
    sav_over_cut_atac_set = over_cut_atac_set |> readr::write_csv("tab/compare_cut_atac_setdb1.csv")
  )

  for_manuscript <- tar_plan(
    mn_fdr_limit = 0.01,
    mn_logfc_limit = 1,

    pdf_pca = plot_pca(db_libs_narrow, metadata) |> gp("atac_pca", width = 5.5, height = 4),

    mn_ctr = c("48h_Ctrl-iPS_Ctrl", "48h_SETDB1-48h_Ctrl"),
    pdf_volcano_1 = mn_plot_volcano(dec_chipseeker_groups_libs_narrow |> filter(contrast == mn_ctr[1]), fdr_limit = mn_fdr_limit, logfc_limit = mn_logfc_limit) |>
      gp(paste0("atac_volcano_", mn_ctr[1]), width = 4, height = 4),
    pdf_volcano_2 = mn_plot_volcano(dec_chipseeker_groups_libs_narrow |> filter(contrast == mn_ctr[2]), fdr_limit = mn_fdr_limit, logfc_limit = mn_logfc_limit) |>
      gp(paste0("atac_volcano_", mn_ctr[2]), width = 4, height = 4),
    pdf_volcano_3 = mn_plot_volcano(dec_chipseeker_groups_libs_narrow |> filter(contrast %in% c("48h_SETDB1-48h_Ctrl", "iPS_SETDB1-iPS_Ctrl")), fdr_limit = mn_fdr_limit, logfc_limit = mn_logfc_limit, ncol = 1) |>
      gp(paste0("atac_volcano_2"), width = 4, height = 7),

    pdf_up_down = fig_up_down_groups_libs_narrow |> 
      gp("atac_up_down", width = 3, height = 4),


    mn_ctr_dir = c("48h_Ctrl-iPS_Ctrl_up", "48h_Ctrl-iPS_Ctrl_down"),

    pdf_groups_chipseeker = fig_groups_chipseeker_narrow |> gp("groups_chipseeker", width = 3.5, height = 3),

    pdf_test_chipseeker_annotations_full = mn_plot_annot_group_test(test_de_chipseeker_groups_libs_narrow, max_or = 2.5) |>
      gp("atac_chipseeker_test_full", width = 3.4, height = 4),
    pdf_test_rmsk_annotations_full = mn_plot_annot_group_test(test_de_rmsk_groups_libs_narrow, max_or = 3) |>
      gp("atac_rmsk_test_full", width = 3.3, height = 4),

    pdf_test_rmsk_annotations_ltr = mn_plot_annot_group_test(test_de_rmsk_ltr_groups_libs_narrow, max_or = 4.2) |>
      gp("atac_rmsk_test_ltr", width = 3.3, height = 4),

    ctr_for_ltr = c("48h_Ctrl-iPS_Ctrl_down", "48h_SETDB1-48h_Ctrl_up", "48h_SETDB1-iPS_SETDB1_down", "48h_SETDB1-iPS_SETDB1_up", "iPS_SETDB1-iPS_Ctrl_up"),
    pdf_rmsk_ltr = dec_rmsk_groups_libs_narrow |>
      filter(contrast_direction %in% ctr_for_ltr & fdr < mn_fdr_limit, abs(log_fc) > mn_logfc_limit, rep_category == "LTR") |>
      mutate(rep_family = str_remove(rep_family, "\\?"), contrast_direction = factor(contrast_direction, levels = ctr_for_ltr)) |> 
      plot_group_peak_annotations(group_var = "contrast_direction", annot_var = "rep_family", n_top = Inf) |> 
      gp("atac_rmsk_ltr", width = 3, height = 4),

    #pdf_test_chipseeker_annotations = mn_plot_annot_group_test(test_dec_chipseeker_groups_libs_narrow |> filter(grp %in% mn_ctr_dir), max_or = 1.3) |>
    # gp("atac_chipseeker_test", width = 2.8, height = 4),
    #pdf_test_rmsk_annotations = mn_plot_annot_group_test(test_dec_rmsk_groups_libs_narrow |> filter(grp %in% mn_ctr_dir), max_or = 1.3) |>
    #  gp("atac_rmsk_test", width = 2.8, height = 4),
    sav_test_chipseeker = test_dec_chipseeker_groups_libs_narrow |> readr::write_csv("tab/test_chipseeker.csv"),
    sav_test_rmsk = test_dec_rmsk_groups_libs_narrow |> readr::write_csv("tab/test_rmsk.csv"),

    pdf_enr_contrasts_logp_down = mn_plot_enr_contrasts_logp(tf_enr, ctr1 = "48h_Ctrl-iPS_Ctrl_down", ctr2 = "48h_SETDB1-iPS_SETDB1_down") |> 
      gp("atac_enr_contrasts_logp_down", width = 4, height = 4),
    pdf_enr_contrasts_logp_up = mn_plot_enr_contrasts_logp(tf_enr, ctr1 = "48h_Ctrl-iPS_Ctrl_up", ctr2 = "48h_SETDB1-iPS_SETDB1_up") |> 
      gp("atac_enr_contrasts_logp_up", width = 4, height = 4),
    pdf_closing_opening = mn_plot_enr_contrasts_logp(tf_enr, ctr1 = "48h_Ctrl-iPS_Ctrl_down", ctr2 = "48h_Ctrl-iPS_Ctrl_up") |> 
      gp("atac_enr_closing_opening", width = 4, height = 4),
    sav_closing_opening = mn_get_enr_contrasts_logp(tf_enr, ctr1 = "48h_Ctrl-iPS_Ctrl_down", ctr2 = "48h_Ctrl-iPS_Ctrl_up") |> 
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~signif(.x, 4))) |> 
      dplyr::select(-mrg) |> 
      readr::write_csv("fig/atac_enr_closing_opening.csv"),

    matsuda = read_matsuda_data(),
    sav_matsuda_de = save_matsuda(matsuda),
    sav_dec_atac = dec_annotatr_groups_libs_narrow |> dplyr::filter(fdr < 0.01 & abs(log_fc) > 1 & contrast == "48h_Ctrl-iPS_Ctrl" & annotation %in% c("Promoter", "1-5 kb")) |> readr::write_csv("tab/dec_atac.csv"),
    sav_matsuda_atac_1 = merge_matsuda(matsuda, dec_annotatr_groups_libs_narrow, "201B7") |> readr::write_csv("tab/matsuda_201B7.csv"),
    sav_matsuda_atac_2 = merge_matsuda(matsuda, dec_annotatr_groups_libs_narrow, "1231A3") |> readr::write_csv("tab/matsuda_1231A3.csv"),
    pdf_matsuda_atac_1 = mn_plot_matsuda_comparison(matsuda, dec_annotatr_groups_libs_narrow, "201B7") |> gp("matsuda_atac_201B7", 5, 4),
    pdf_matsuda_atac_2 = mn_plot_matsuda_comparison(matsuda, dec_annotatr_groups_libs_narrow, "1231A3") |> gp("matsuda_atac_1231A3", 5, 4)
  )

  for_shiny <- tar_plan(
    sav_shiny = prepare_for_shiny(dbc_full_libs_narrow, dbc_groups_libs_narrow, metadata, terms, gns)
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
    motif_analysis,
    cuttag_comparison,
    for_report,
    for_manuscript,
    for_shiny
  )
  
}