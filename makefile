.PHONY: doc smmuty
MANDOC := man/*.Rd
RCODES := R/*.R
SMMUTYCODES := inst/smmuty/smmuty/*.py inst/smmuty/bin/*.py

all: README.md doc smmuty

README.md: README.Rmd
	Rscript -e "devtools::build_readme()"

# cannot stop running
# $(MANDOC): $(RCODES)
	# Rscript -e "devtools::document()"

doc: $(RCODES)
	Rscript -e "devtools::document()"

smmuty: $(SMMUTYCODES)
	cd inst; \
  pip install -e smmuty
