
required_packages <- c(
  "shiny", "factoextra", "proxy", "ggrepel", "ggforce", "progressr",
  "shinyjs", "shinydashboard", "DT", "data.table", "tidyverse",
  "furrr", "dbscan", "remotes", "promises", "future"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cran.rstudio.com/")
  }
}


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", dependencies = TRUE)
}


bioconductor_packages <- c("GenomicRanges", "plyranges")


for (pkg in bioconductor_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

suppressPackageStartupMessages({
  lapply(c(required_packages, bioconductor_packages), function(pkg) {
    library(pkg, character.only = TRUE)
  })
})


