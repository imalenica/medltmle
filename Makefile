md:
	Rscript -e "rmarkdown::render('README.Rmd')"

site:
	Rscript -e "pkgdown::build_site()"

check:
	Rscript -e "devtools::check()"

test:
	Rscript -e "devtools::test()"

doc:
	Rscript -e "devtools::document()"

cov:
	Rscript -e "covr::package_coverage(type = 'all', combine_types = FALSE, line_exclusions = list('R/plots.R', 'R/utils.R'))"
