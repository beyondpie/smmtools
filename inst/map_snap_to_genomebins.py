### Get Sparse Matrix from genome bin.

from typing import List, Dict, Tuple
import h5py
import collections
import argparse


def readGenomeSizeFromTxt(fname) -> Dict[str, int]:
    """
    Read genome information.

    Args:
    -----
        fname: a txt file contains genome information

    Returns:
    -----
        A dictionary contains SQ as key and SL as value, otherwise None
    """
    # first check if fname is a bam file
    res = dict()
    with open(fname) as fin:
        for line in fin:
            elems = line.split()
            chr_name = elems[0]
            chr_len = int(elems[1])
            res[chr_name] = chr_len
    return res


def getBinsFromGenomeSize(
    genome_dict: Dict[str, int], bin_size: int
) -> Dict[Tuple[str, int, int], int]:
    """Create a dictionary contains all bins of the same size across the genome

    Attributes:
        binSize: bin size (i.e. 5000)
        genomeDict: a dictionary contains chromosome sizes
    Return:
        A dictionary contains all bins and its index (start from 1)
    """
    bin_dict = collections.OrderedDict()
    i = 1
    for _chrom in genome_dict:
        for _start in range(1, genome_dict[_chrom], bin_size):
            _end = min(_start + bin_size - 1, genome_dict[_chrom])
            _binId = (_chrom, _start, _end)
            bin_dict[_binId] = i
            i = i + 1
    return bin_dict


def getBarcodesFromSnap(fname: str) -> List[str]:
    """Read barcodes from a snap file

    Attributes:
        fname - a snap-format file

    Return:
        a dictionary contains all barcode in the snap file
    """
    with h5py.File(fname, "r") as f:
        barcodes: List[str] = [item.decode() for item in f["BD/name"]]
    return barcodes


def snap_bmat(
    snap_file,
    bin_size: int,
    genome_fname: str,
    barcodes_file: str,
    bmat_outf: str,
    col_outf: str,
):
    """
    Pre-processing to create a snap file from a bam that
    contains alignments or a bed file that contains fragments.

    Args:
    --------
    snap_file:
        a snap format file.
    bin_size_list:

    Return:
    - file 1:
      - each row: barcode,[genome_bin_id:count,]
      - genome_bin_id start from one.
    - file 2, bed file to record the the column content.
      - each row: chr, start, end
    """
    f = h5py.File(snap_file, "a", libver="earliest")
    ## load genome dict
    genome_dict: Dict[str, int] = readGenomeSizeFromTxt(genome_fname)
    # extract the barcodes
    # barcode_dict = getBarcodesFromSnap(snap_file=snap_file)
    barcode_dict: Dict[str, int] = collections.OrderedDict()
    with open(barcodes_file, "r") as f:
        lines = f.readlines()
    for i, l in enumerate(lines):
        barcode = l.strip()
        barcode_dict[barcode] = i

    bin_dict: Dict[Tuple[str, int, int], int] = getBinsFromGenomeSize(
        genome_dict=genome_dict, bin_size=bin_size
    )

    idxList = collections.defaultdict(list)
    # barcode index list
    idyList = collections.defaultdict(list)
    # bin index list
    countList = collections.defaultdict(list)
    # number of count

    lines: List[str] = []
    for barcode_id, barcode_encode in enumerate(f["BD"]["name"]):
        l: List[str] = []
        barcode: str = barcode_encode.decode()
        if barcode not in barcode_dict:
            continue
        l.append(barcode)
        _chroms = f["FM"]["fragChrom"][
            (f["FM"]["barcodePos"][barcode_id] - 1) : (
                f["FM"]["barcodePos"][barcode_id]
                + f["FM"]["barcodeLen"][barcode_id]
                - 1
            )
        ]
        _chroms = [item.decode() for item in _chroms]
        _start = f["FM"]["fragStart"][
            (f["FM"]["barcodePos"][barcode_id] - 1) : (
                f["FM"]["barcodePos"][barcode_id]
                + f["FM"]["barcodeLen"][barcode_id]
                - 1
            )
        ]
        _len = f["FM"]["fragLen"][
            (f["FM"]["barcodePos"][barcode_id] - 1) : (
                f["FM"]["barcodePos"][barcode_id]
                + f["FM"]["barcodeLen"][barcode_id]
                - 1
            )
        ]
        bins: Dict[Tuple[str, int, int], int] = collections.defaultdict(lambda: 0)
        for bin_chr, bin_start, bin_end in list(zip(_chroms, _start, _start + _len)):
            for bin_pos in set(
                [
                    int(bin_start / bin_size) * bin_size + 1,
                    int(bin_end / bin_size) * bin_size + 1,
                ]
            ):
                bins[(bin_chr, bin_pos, bin_pos + bin_size - 1)] += 1

        for key in bins:
            if key in bin_dict:
                l.append(f"{bin_dict[key]:bins[key]}")
        lines.append(",".join(l))
    ## output result
    if len(lines) > 0:
        ## output genome bins
        with open(col_outf, "w") as f:
            f.write("\n".join([f"{chrom}\t{spos}\t{epos}" for chrom, spos, epos in bin_dict.keys()]))
        ## output bmat
        with open(bmat_outf, "w") as f:
            f.write("\n".join(lines))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Map Snapfile to genomeBins for the given barcodes."
    )
    parser.add_argument("--snapfile", type=str, required=True)
    parser.add_argument("--binsize", type=lambda x:int(float(x)), required=True)
    parser.add_argument("--genomefname", type=str, required=True)
    parser.add_argument("--barcodesfile", type=str, required=True)
    parser.add_argument("--bmatoutf", type=str, required=True)
    parser.add_argument("--coloutf", type=str, required=True)
    args = parser.parse_args()
    snap_bmat(
        snap_file=args.snapfile,
        bin_size=args.binsize,
        genome_fname=args.genomefname,
        barcodes_file=args.barcodesfile,
        bmat_outf=args.bmatoutf,
        col_outf=args.coloutf,
    )
