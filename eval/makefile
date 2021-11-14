## make sure to run the file under conda base environment.

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
bmat_fnms := $(foreach s, $(samples), out/$(s)_bmat.h5)
scrublet_fnms := $(foreach s, $(samples), out/$(s)_scrublet.csv)
barcodes_after_doublet := $(foreach s, $(samples), out/$(s)_barcodes_after_doublet.csv)
harmony_matrix_fnm := harmony_matrix.h5
RB_cluster_fnm := RB_cluster.csv

EMPTY :=
SPACE := $(EMPTY) $(EMPTY)
COMMA := ,
bmatfnms_str := $(subst $(SPACE), $(COMMA), $(foreach s, $(bmat_fnms), $(notdir $s)))
barcodefnms_str := $(subst $(SPACE), $(COMMA), $(foreach s, $(barcodes_after_doublet), $(notdir $s)))
samples_str := $(subst $(SPACE), $(COMMA), $(samples))


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
	bedtools sort -i out/$*_fragment.tsv | bgzip  > $@
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

$(bmat_fnms) : out/%_bmat.h5 : out/%_barcodes_filtered.txt out/%_raw.h5
	$(Rscript) $(script_dir)/getTileMatrix.R --rawH5File out/$*_raw.h5 \
             --outdir out --genome $(genome) --sampleName $* \
             --barcodesFile out/$*_filtered.txt --tileSize 5000 \
             --outfnm $*_bmat.h5

$(scrublet_fnms) : out/%_scrublet.csv : out/%_bmat.h5
	$(Rscript) $(script_dir)/getDoubletScore.R --tileMatrixFile $< \
             --outScrubletFile $@ --path_to_python $(python) \
             --zthres 1.65 --expDoubletRate 0.06 --minCount 3 --minCells 5 \
             --minConvPctl 85 --nPC 30

$(barcodes_after_doublet) : out/%_barcodes_after_doublet.csv : out/%_scrublet.csv
	$(R) --vanilla -q <<EOF
	doubletScores <- read.table(file = $<, header = TRUE, sep = ",", comment.char = "")
	out <- doubletScores$barcode[doubleScores$scrubletScore <= doubletScores$thresGMM]
	write.table(x = out, file = $@, quote = FALSE, row.names = F, col.names = F, sep = ",")
	EOF

$(harmony_matrix_fnm) : $(bmat_fnms) $(barcodes_after_doublet)
	$(Rscript) $(script_dir)/getHarmonyMatrix.R --sampleNames samples_str \
             --tileMatrixDir out --barcodeDir out \
             --tileMatrixFNames $(bmatfnms_str) --barcodeFNames $(barcodefnms_str)\
             --nLandmark 10000 --nSampleRegrn 2000 --nPCDiffusionMaps 50 \
             --nPCHarmony 30 --saveHarmonyFile $@

$(RB_cluster_fnm) : $(harmony_matrix_fnm)
	$(Rscript) $(script_dir)/getclustering.R --dimReductMatrixFile $< \
             --Knn 50 --path_to_python $(python) --reso 1.0 --partitionType RB \
             --outClusterFile $@


clean:
	-rm $(foreach s, $(samples), out/$(s)_fragment.tsv)
