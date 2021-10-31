.PHONY: doc
MANDOC := man/*.Rd
RCODES := R/*.R
all: README.md doc

README.md: README.Rmd
	Rscript -e "devtools::build_readme()"

# cannot stop running
# $(MANDOC): $(RCODES)
	# Rscript -e "devtools::document()"

.PHONY: doc
doc: $(RCODES)
	Rscript -e "devtools::document()"

