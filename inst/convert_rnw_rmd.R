library(MCMCglmm)
library(MASS)
library(lattice)
set.seed(1)
setwd("/Users/jhadfiel/Work/MCMCglmm/bookdown")

bookdown::render_book("index.Rmd", "bookdown::gitbook")
# generates book

setwd("/Users/jhadfiel/Work/MCMCglmm")
pkgdown::build_site()
# builds package documentation

dir.create("docs/course-notes", recursive = TRUE, showWarnings = FALSE)
fs::dir_copy("bookdown/_book", "docs/course-notes", overwrite = TRUE)
# moves book to doc







