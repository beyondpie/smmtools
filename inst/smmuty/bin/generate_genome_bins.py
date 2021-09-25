import argparse
from typing import Dict, List, Tuple
from smmuty import readGenomeSizeFromTxt, getBinsFromGenomeSize

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog = "Generate Genome Bins Bed File."
    )
    parser.add_argument("--binsize", type = lambda x:int(float(x)), required= True)
    parser.add_argument("--genomef", type = str, required= True)
    parser.add_argument("--outf", type = str, required = True)
    args = parser.parse_args()
    
    genome_dict :Dict[str, int] = readGenomeSizeFromTxt(fname = args.genomef)
    bin_dict: Dict[Tuple[str, int, int], int] = getBinsFromGenomeSize(
        genome_dict = genome_dict, bin_size = args.binsize
    )
    lines: List[str] = []
    for key, val in bin_dict.items():
        l = "\t".join([key[0], str(key[1]), str(key[2]), str(val)])
        lines.append(l)

    with open(args.outf, "w") as f:
        f.writelines("\n".join(lines))
