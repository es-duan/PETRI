# scripts/install_cran_pkgs.R

# Specify the CRAN mirror
local({r <- getOption("repos")
       r["CRAN"] <- "https://cloud.r-project.org"
       options(repos=r)
})

# Install funkyheatmap if it's not already installed
if (!requireNamespace("funkyheatmap", quietly = TRUE)) {
  install.packages("funkyheatmap")
}
