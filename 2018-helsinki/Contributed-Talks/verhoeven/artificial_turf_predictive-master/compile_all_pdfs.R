rm(list = ls())

# Make sure that working directory is set to the location of this R script!
# i.e.
#setwd("C:/artificial_turf_predictive")

# needed R-packages
require(rmarkdown)
require(knitr)
require(data.table)
require(ggplot2)
require(cowplot)
require(rstan)
require(forecast)


# compile all .rmd files in the required order

rmarkdown::render("wrangle_data.Rmd", 
                  output_format = "pdf_document", envir = new.env())

dir.create("FITS", showWarnings = FALSE)

rmarkdown::render("simulate_dynamic.Rmd", 
                  output_format = "pdf_document", envir = new.env())

rmarkdown::render("fit_models.Rmd", 
                  output_format = "pdf_document", envir = new.env())

rmarkdown::render("analyse_results.Rmd", 
                  output_format = "pdf_document", envir = new.env())

rmarkdown::render("artificial_turf_predictive.Rmd", 
                  output_format = "pdf_document", envir = new.env())

