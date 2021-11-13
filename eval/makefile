.ONESHELL:
.PHONY: smmuty all clean

smmuty:
	cd ../inst; pip install -e $@ ; cd -

conda_env := $(HOME)/miniconda
R := $(conda_env)/bin/R
Rscript := $(conda_env)/bin/R
python := $(conda_env)/bin/python

script_dir := $(HOME)/git-recipes/10xmm/script
samples := QY_1287 QY_1288
genome := 10mm

raw_barcode_fnms := $(foreach s, $(samples), out/$(s)_barcode)
fragment_fnms := $(foreach s, $(samples), out/$(s)_fragment.tsv.gz)
QC_fnms := $(foreach s, $(samples), out/$(s)_QC.csv)
raw_h5_file := $(foreach s, $(samples), out/$(s)_raw.h5)
sumfragment_file := $(foreach s, $(samples), out/$(s)_sumFragment.csv)
filtered_barcode_fnms := $(foreach s, $(samples), out/$(s)_barcodes_filtered.txt)
bmat_fnms := $(foreach s, $(samples), out/$(s)_)

all: $(raw_barcode_fnms) $(fragment_fnms) 

$(raw_barcode_fnms) : out/%_barcode : bam/%.bam
	${python} <<EOF
	from smmuty.snaptools import getBarcodesFromBam as f
	f("$<", "$@")
	EOF

$(fragment_fnms) : out/%_fragment.tsv.gz : bam/%.bam out/%_barcode
	${python} << EOF
	from smmuty.snaptools import SnapToolsBamTo10xFragmentBed as f
	f("bam/$*.bam", "out/$*_fragment.tsv", "out/$*_QC.tsv", "out/$*_barcode")
	EOF
	bedtools sort -i out/$*_fragment.tsv | bgzip < > $@
	tabix --preset=bed $@

$(raw_h5_file) : out/%_raw.h5 : out/%_fragment.tsv.gz
	$(Rscript) $(script_dir)/getTabixToH5.R --tabixFile $< --rawH5File $@ \
      --genome $(genome) --sampleName $* --nChunk 3

$(sumfragment_file) : out/%_sumFragment.csv : out/%_raw.h5
	$(Rscript) $(script_dir)/getSumFragment.R --rawH5File $< --sumFragFile $@ \
      --genome $(genome) --sampleName $*

$(filtered_barcode_fnms) : out/%_barcodes_filtered.txt : out/%_sumFragment.csv
	$(Rscript) $(script_dir)/filterBarcodesFromQC.R --qcf $< --outf $@ \
      --TSSE 10 --nFrag 1000 --sampleName $*

$()

clean:
	-rm $(foreach s, $(samples), out/$(s)_fragment.tsv)
