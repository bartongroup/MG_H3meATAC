library(dexdash)
library(readr)

# datatable needs this option to display Inf correctly

options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))
options(dplyr.summarise.inform = FALSE)

# Source all R files

files_R <- list.files(c("R", "modules"), pattern = "*.R$", full.names = TRUE)
sr_ <- sapply(files_R, source)


# Read config json file

config <- jsonlite::read_json("shiny_config.json")

# Read all data, make it available to all modules

shd <- sh_read_data(c("dexset", "features", "fterms"), config, with_progress = FALSE)

dexdash::run_app(shd$dexset, shd$features, shd$fterms, title = config$title,
                 x_variable = "group", colour_variable = "replicate",
                 value_variable = "abundance", volma_fdr_limit = 0.01, volma_logfc_limit = 1)

