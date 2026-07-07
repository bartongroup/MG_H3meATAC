#!/usr/bin/env Rscript

# usage:
# diffbind_count_cluster.R diffbind_metadata.csv broad

args = commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)

in_file <- args[1]
type <- args[2]

library(DiffBind)

read_and_count_diffbind <- function(meta_db, min_overlap = 2) {
  
  ## Read all peaks & metadata 
  db <- DiffBind::dba(sampleSheet = meta_db)
  
  ## Build per-condition consensus peaksets
  peakset.cons <- DiffBind::dba.peakset(
    db,
    consensus = DBA_CONDITION,
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
    bParallel = FALSE
  )
  db
}

spikein_normalization <- function(dbc) {
  DiffBind::dba.normalize(dbc, normalize = DiffBind::DBA_NORM_LIB, spikein = TRUE)
}


meta <- read.csv(in_file)

dbc <- read_and_count_diffbind(meta, min_overlap = 2)
saveRDS(dbc, file = paste0("dbc_", type, ".rds"))

dbs <- spikein_normalization(dbc)
saveRDS(dbs, file = paste0("dbc_", type, "_spike.rds"))
