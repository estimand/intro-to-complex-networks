
all: network_analysis.pdf

%.pdf: %.R
	@Rscript -e "rmarkdown::render('$<')"

