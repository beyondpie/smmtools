.PHONY: doc smmuty install
MANDOC := man/*.Rd
RCODES := R/*.R
CPPCODES := src/*.cpp
SMMUTYCODES := inst/smmuty/smmuty/*.py inst/smmuty/bin/*.py

all: README.md doc smmuty

README.md: README.Rmd
	Rscript -e "devtools::build_readme()"

# cannot stop running
# $(MANDOC): $(RCODES)
	# Rscript -e "devtools::document()"

doc: $(RCODES)
	-rm src/*.o
	-rm src/*.so
	Rscript -e "devtools::document(pkg = '.')"

smmuty: $(SMMUTYCODES)
	cd inst; \
  python3 -m pip install -e smmuty

install:
	-rm src/*.o
	-rm src/*.so
	Rscript -e 'devtools::install_github("beyondpie/smmtools")'
