rm(list=ls())
#source("/Users/jhadfiel/Work/MCMCglmm/inst/buildbook.R")
build_site<-FALSE
rm.cache<-FALSE

setwd("/Users/jhadfiel/Work/MCMCglmm/bookdown")

if(rm.cache){
 system("rm -r /Users/jhadfiel/Work/MCMCglmm/bookdown/_bookdown_files")
}

bookdown::render_book("1.intro.Rmd", "bookdown::gitbook")
# generates book

if(build_site){
setwd("/Users/jhadfiel/Work/MCMCglmm")
pkgdown::build_site()
# builds package documentation

dir.create("docs/course-notes", recursive = TRUE, showWarnings = FALSE)
fs::dir_copy("bookdown/_book", "docs/course-notes", overwrite = TRUE)
# moves book to doc
}
