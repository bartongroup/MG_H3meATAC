# create a GRanges object for a single region
make_region <- function(chr, start, end) {
  gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(seqnames = as.character(chr), ranges = IRanges::IRanges(start, end)))
  genome(gr) <- "GRCh38"
  return(gr)
}

# expand a region for plotting
expand_region <- function(gr, f = 1, fl = NULL, fr = NULL) {
  w <- width(gr) / 2
  mid <- (start(gr) + end(gr)) / 2
  if(is.null(fl)) fl <- f
  if(is.null(fr)) fr <- f
  start <- as.integer(mid - fl * w)
  end <- as.integer(mid + fr * w)
  gr <- make_region(seqnames(gr), start, end)
  return(gr)
}

get_peak_region <- function(peaks, id, f = 2) {
  peaks |> 
    filter(peak_id == id) |> 
    GenomicRanges::makeGRangesFromDataFrame() |> 
    expand_region(f)
}

# set nice colours to tracks
track_colour <- function(smpl, meta) {
  idx <- meta |> 
    select(sample, group) |> 
    filter(sample == smpl) |> 
    pull(group) |> 
    as.integer()

  okabe_ito_palette[idx]
}


# name of the data file, either bw or bam
data_track_file <- function(smpl, config, meta, type = c("bam", "bw")) {
  type <- match.arg(type)
  rsmpl <- meta |> 
    distinct(sample, raw_sample) |> 
    filter(sample == smpl) |> 
    pull(raw_sample) |> 
    as.character()

  if(type == "bam") {
    file <- file.path(config$remote_path, config$dir$bam, str_glue("{rsmpl}.noMT.bam"))
  } else if (type == "bw") {
    file <- file.path(config$remote_path, config$dir$bigwig, str_glue("{rsmpl}.bw"))
  }
  stopifnot(file.exists(file))
  return(file)
}


# make data tracks from a list of sample names
make_data_tracks <- function(samples, config, meta, type, chr, genome = "GRCh38", ymax = NULL) {
  if(is.null(ymax)) {ylim=NULL} else {ylim=c(0, ymax)}
  dtracklist <- map(samples, function(smpl) {
    file <- data_track_file(smpl, config, meta, type)
    dtrack <- Gviz::DataTrack(
      range = file,
      type = "h",
      name = smpl,
      ylim = ylim,
      chromosome = chr,
      genome = genome,
      col = track_colour(smpl, meta),
      fontsize = 8
    )
  })
  names(dtracklist) <- samples
  return(dtracklist)
}

# make peak tracks from a list of sample names
make_peak_tracks <- function(samples, peaks, chrom) {
  ptracklist <- map(samples, function(smpl) {
    peak <- peaks |>
      filter(sample == smpl & chr == chrom) |> 
      GenomicRanges::makeGRangesFromDataFrame()
    tr <- Gviz::GeneRegionTrack(peak, name = "MACS", showId = FALSE, geneSymbol = FALSE, fill = "red")
  })
  names(ptracklist) <- paste0(samples, "_peaks")
  return(ptracklist)
}

# read a genomic region from a bw file
get_region_from_file <- function(file, region) {
  sel <- rtracklayer::BigWigSelection(ranges = region)
  rtracklayer::import(file, selection = sel)
}


# need global variables mart and fm
plot_region <- function(samples, config, meta, region, mart, ymax = NULL, peaks = NULL, type = "bam",
                        with.tracks = FALSE, stacking = "squish") {
  samples <- as.character(samples)
  
  options(ucscChromosomeNames = FALSE)
  chr <- as.character(seqnames(region))

  fm <- Gviz:::.getBMFeatureMap()
  
  # guess what
  # need to fudge chromosome name, or it will all go to hell
  idt <- Gviz::IdeogramTrack(chromosome = paste0("chr", chr), genome = "hg38")
  levels(idt@bandTable$chrom) <- sub("^chr", "", levels(idt@bandTable$chrom), ignore.case = TRUE)
  
  gat <- Gviz::GenomeAxisTrack()
  

  if(with.tracks) {
    grt <- tryCatch({
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      Gviz::GeneRegionTrack(txdb, chromosome = str_glue("chr{chr}"), start = start(region), end = end(region))
      #BiomartGeneRegionTrack(chromosome = chr, start = start(region), end = end(region), featureMap = fm, biomart = mart, stacking = stacking)
    }, error=function(err) {
      print("Error creating genomic tracks: skipped.")
      return(NULL)
    })
  } else {
    grt <- NULL
  }
  #if(!is.null(gr)) grt <- GeneRegionTrack(range=gr, name="Ensembl", stacking=stacking)
  
  # to find common ymax across tracks, we have to load and analyse them
  # right now it works only for BigWig files, not for BAM
  if(!is.null(ymax) && ymax == "auto") {
    if(type == "bam") stop('ymax="auto" does not work with BAM files.')
    max.track <- sapply(samples, function(smpl) {
      file <- data_track_file(smpl, config, meta, type)
      r <- get_region_from_file(file, region)
      max(r$score)
    })
    ymax <- max(max.track)
  }
  
  # read data tracks
  dtracks <- make_data_tracks(samples, config, meta, type, chr, ymax = ymax)
  
  # create MACS peak tracks
  if(!is.null(peaks)) {
    ptracks <- make_peak_tracks(samples, peaks, chr)
    n <- length(samples)
    tr <- list()
    for(i in 1:n) tr <- c(tr, dtracks[i], ptracks[i])
    if(length(dtracks) > n) tr <- c(tr, dtracks[(n+1):length(dtracks)])
    dtracks <- tr
  }
  
  #tracks <- append(list(idt, gat), dtracks)
  tracks <- append(list(gat), dtracks)
  if(!is.null(grt)) tracks <- append(tracks, list(grt))
  tryCatch({
    Gviz::plotTracks(tracks, transcriptAnnotation = "symbol", from = start(region), to = end(region), chromosome = chr)  
  }, error=function(err) {
     print("Error plotting tracks (too many tracks??)")
     if(!is.null(grt)) {
        print("  trying dense stacking")
        #grt <- GeneRegionTrack(range=gr, name="Ensembl", stacking="dense")
        grt <- BiomartGeneRegionTrack(chromosome=chr, start=start(region), end=end(region), featureMap=fm, biomart=mart, stacking="dense")
        tracks[length(tracks)] <- grt
        plotTracks(tracks, transcriptAnnotation = "symbol", from=start(region), to=end(region))
     } else {
       print("Giving up")
     }
  })
  
}


plot_region_peak_id <- function(dec, peakid, config, meta, peaks, mart, f = 2) {
    if(is.numeric(peakid)) {
    peakid <- dec$peak_id[peakid]
  }

  pk <- dec |> 
    filter(peak_id == peakid)
  stopifnot(nrow(pk) == 1)

  #sample <- pk$sample[1] |> as.character()
  samples <- meta |> filter(as.character(replicate) == "1") |> pull(sample) |> as.character()

  region <- make_region(pk$chr, pk$start, pk$end) |> 
    expand_region(f = f)

  plot_region(samples, config, meta, region, mart, peaks = peaks)
}