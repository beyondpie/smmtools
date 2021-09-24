from typing import List, Dict, Tuple
import collections

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

