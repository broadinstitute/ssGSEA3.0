## Renders a polished, standalone PDF copy of ssgsea2-vs-ssgsea3-comparison.Rmd
## (Arial via xelatex, styled headings -- matching the lab's other report
## PDFs). Not part of the package build/R CMD check -- the package's own
## vignette output is HTML (ssgsea3-vignette-style.css), which has no LaTeX
## dependency. Run from within vignettes/:
##   Rscript render-pdf.R
devtools::load_all("..", quiet = TRUE)
rmarkdown::render(
  "ssgsea2-vs-ssgsea3-comparison.Rmd",
  output_format = rmarkdown::pdf_document(
    latex_engine = "xelatex",
    toc = TRUE,
    number_sections = FALSE,
    includes = rmarkdown::includes(in_header = "ssgsea3-vignette-pdf-header.tex")
  ),
  output_file = "ssgsea2-vs-ssgsea3-comparison.pdf"
)
