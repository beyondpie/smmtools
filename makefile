.PHONY: doc smmuty install
MANDOC := man/*.Rd
RCODES := R/*.R
CPPCODES := src/*.cpp
SMMUTYCODES := inst/smmuty/smmuty/*.py inst/smmuty/bin/*.py

all: README.md doc

README.md: README.Rmd
	Rscript -e "devtools::build_readme()"

# cannot stop running
# $(MANDOC): $(RCODES)
	# Rscript -e "devtools::document()"

doc: $(RCODES)
	-rm src/*.o
	-rm src/*.so
	Rscript -e "devtools::document(pkg = '.')"

install_local_smmuty: $(SMMUTYCODES)
	cd inst; \
  python3 -m pip install -e smmuty

upload_smmuty:
	cd inst; \
	python setup.py sdist bdist_wheel; \
	twine upload dist/*

install:
	-rm src/*.o
	-rm src/*.so
	Rscript -e 'devtools::install_github("beyondpie/smmtools")'

install_local:
	Rscript -e 'install.packages(".", repos = NULL, type = "source")'

