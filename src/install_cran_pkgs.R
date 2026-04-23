# scripts/install_cran_pkgs.R

# Specify the CRAN mirror
local({r <- getOption("repos")
       r["CRAN"] <- "https://cloud.r-project.org"
       options(repos=r)
})

# Install packages
packages_to_install <- c("funkyheatmap", "deSolve")

for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
