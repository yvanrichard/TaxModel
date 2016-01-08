ex_dirs = DD-example size-example PPMR-example
plots = $(addsuffix /Plots.pdf, $(ex_dirs))

$(plots): %.pdf: %.Rmd
	cd $(@D); Rscript -e 'require(knitr); knit($(<F))'
