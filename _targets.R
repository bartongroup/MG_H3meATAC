suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
  library(qs2)
})

required_packages <- read.delim("R/packages.txt", header = FALSE, col.names = "name")$name
tar_option_set(packages = required_packages, format = "qs")
options(tidyverse.quiet = TRUE, dplyr.summarise.inform = FALSE)

# Create dirs if necessary
for (d in c("tab", "fig", "doc", "doc/images", "cache")) if (!dir.exists(d)) dir.create(d)

# For interactive session only
if (interactive()) sapply(required_packages, library, character.only = TRUE)

# Source all functions from .R files
files_R <- list.files(c("R", "targets"), pattern = "*.R$", full.names = TRUE)
sr_ <- sapply(files_R, source)

sesinfo <- list(
  tar_target(session_info, sessionInfo())
)

# prevent other packages from stealing from tidyverse
filter <- dplyr::filter
rename <- dplyr::rename
select <- dplyr::select
slice <- dplyr::slice
count <- dplyr::count


# Add session info
sessinfo <- list(
  tar_force(session_info, sessioninfo::session_info(), force = TRUE)
)

c(
  sessinfo,
  targets_main()
)

