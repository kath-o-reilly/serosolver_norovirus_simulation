#!/usr/bin/env Rscript
# make_description.R
# Kamil Slowikowski
# 2019-05-15

# Get installed packages.
ip <- as.data.frame(installed.packages()[ , c(1, 3:4)])
ip <- ip[is.na(ip$Priority), 1:2, drop = FALSE]
packages <- sort(unique(as.character(ip$Package)))

is_github <- function(x) {
  exists("GithubRepo", packageDescription(x))
}

ix <- sapply(packages, is_github)
github_packages <- packages[ix]
cran_packages <- packages[!ix]

add_github <- function(x) {
  pd <- packageDescription(x)
  if (exists("GithubRepo", pd)) {
    x <- sprintf("github::%s/%s", pd$GithubUsername, pd$GithubRepo)
  }
  return(x)
}
github_packages <- unlist(lapply(github_packages, add_github))

out_file <- "DESCRIPTION"

# Template for the DESCRIPTION file.
temp <- "Package: MyPackages
Version: 0.0.0.9000
Title: My Packages
Description: My Packages
Encoding: UTF-8
Depends:
    R (>= 3.0.0)
biocViews:"
cat(temp, file = out_file)

cat("
Remotes:
    ",
    paste(github_packages, collapse = ",\n    "),
    file = out_file,
    append = TRUE,
    sep = ""
)

cat("
Imports:
    ",
    paste(packages, collapse = ",\n    "),
    file = out_file,
    append = TRUE,
    sep = ""
)