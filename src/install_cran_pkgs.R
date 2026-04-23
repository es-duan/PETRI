# scripts/install_cran_pkgs.R

local({r <- getOption("repos")
       r["CRAN"] <- "https://cloud.r-project.org"
       options(repos=r)
})

# Force R to use the active Conda environment's library path
# .libPaths()[1] will point to the .snakemake/conda/.../lib/R/library folder
conda_lib_path <- .libPaths()[1]

# Install funkyheatmap ONLY into the isolated Conda path
if (!requireNamespace("funkyheatmap", quietly = TRUE)) {
  install.packages("funkyheatmap", lib = conda_lib_path)
}
