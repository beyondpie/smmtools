## QY_1287, QY_1288

import smmuty
import os
from typing import Dict

from smmuty.snaptools import SnapToolsBamTo10xFragmentBed, getBarcodesFromBam

QY_1287_bam_file = "./bam/QY_1287.bam"
# barcode_cov = getBarcodesFromBam(input_bam= QY_1287_bam_file)
# with open("./out/barcode_cov.tsv", "w") as f:
#     lines = [f"{k}\t{v}" for k, v in barcode_cov.items()]
#     f.write("\n".join(lines))

with open("./out/barcode_cov.tsv", "r") as f:
    lines = [l.strip() for l in f.readlines()]
    barcodes:Dict[str, float] = {l.split("\t")[0]: float(l.split("\t")[1]) for l in lines}

SnapToolsBamTo10xFragmentBed(bam_file = QY_1287_bam_file,
                             outf = "./out/QY_1287_fragment.tsv",
                             qc_file = "./out/QY_1287_QC.csv",
                             barcode_with_cov= barcodes,
                             sep = ",",
                             step = 1000)
