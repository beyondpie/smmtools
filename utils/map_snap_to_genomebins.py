### Get Sparse Matrix from genome bin.
### Output:
### - file 1:
###   - each row: barcode,[genome_bin_id,]
### - file 2, bed file
###   - each row: chr, start, end

from typing import List
import h5py
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
    res = dict();
    with open(fname) as fin:
        for line in fin:
            elems = line.split();
            chr_name = elems[0];
            chr_len = int(elems[1]);
            res[chr_name] = chr_len;
    return res

def getBinsFromGenomeSize(genome_dict, bin_size):
    """Create a dictionary contains all bins of the same size across the genome
    
    Attributes:
        binSize: bin size (i.e. 5000)
    
        genomeDict: a dictionary contains chromosome sizes
    
    Return:
        A dictionary contains all bins and its index (start from 1)
    """
    bin_dict = collections.OrderedDict();
    i = 1;
    for _chrom in genome_dict:
    	for _start in range(1, genome_dict[_chrom], bin_size):
            _end = min(_start + bin_size - 1, genome_dict[_chrom]);
            _binId = (_chrom , _start, _end);
            bin_dict[_binId] = i;
            i = i +1;
    return bin_dict;

def getBarcodesFromSnap(fname: str) -> List[str]:
    """Read barcodes from a snap file
    
    Attributes:
        fname - a snap-format file

    Return:
        a dictionary contains all barcode in the snap file
    """
    with h5py.File(fname, 'r') as f:
        barcodes: List[str] =  [item.decode() for item in f["BD/name"]]
    return barcodes


def snap_bmat(snap_file,
              bin_size_list: List[int],
              genome_fname:str):
    """
    Pre-processing to create a snap file from a bam that
    contains alignments or a bed file that contains fragments.

    Args:
    --------
    snap_file: 
        a snap format file.
    bin_size_list: 
    """
    # check if snap_file is a snap-format file
    file_format = snaptools.utilities.checkFileFormat(snap_file);
    # create the bin list
    f = h5py.File(snap_file, "a", libver='earliest');
    ## load genome dict
    genome_dict: Dict[str, int] = readGenomeSizeFromTxt(genome_fname)
    # extract the barcodes
    barcode_dict = getBarcodesFromSnap(snap_file);
    bin_dict_list = collections.defaultdict(dict);

    for bin_size in bin_size_list:
        bin_dict = snaptools.utilities.getBinsFromGenomeSize(genome_dict, bin_size);
        bin_dict_list[bin_size] = bin_dict;

    num_barcode = len(barcode_dict);
    if verbose:
        print("===== reading the barcodes and bins ======");    
        print(("@AM\tnBinSize:%d"%len(list(bin_dict_list.keys()))));
        print("@AM\tbinSizeList: %s" % str(list(bin_dict_list.keys())));    
        for bin_size in list(bin_dict_list.keys()):
            print(("@AM\tbinSize:%d\tnBin:%d"%(bin_size, len(bin_dict_list[bin_size]))));
        
    idxList   = collections.defaultdict(list);        # barcode index list
    idyList   = collections.defaultdict(list);        # bin index list
    countList = collections.defaultdict(list);        # number of count

    barcode_id = 0
    for barcode in f["BD"]["name"]:  
        _chroms = f["FM"]["fragChrom"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)];
        _chroms = [item.decode() for item in _chroms];
        _start = f["FM"]["fragStart"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)]
        _len = f["FM"]["fragLen"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)]
        frag_list_uniq = list(zip(_chroms, _start, _start + _len));    

        for bin_size in bin_dict_list:
            bin_dict = bin_dict_list[bin_size];
            bins = collections.defaultdict(lambda : 0);
            for item in frag_list_uniq:
                bin_chr = item[0];
                for bin_pos in set([int(item[1]/bin_size) * bin_size + 1, int(item[2]/bin_size) * bin_size + 1]):
                    bins[(bin_chr, bin_pos, bin_pos + bin_size - 1)] += 1;
        
            for key in bins:
                if key in bin_dict and barcode.decode() in barcode_dict:
                    idyList[bin_size].append(bin_dict[key]);
                    countList[bin_size].append(bins[key]);
                    idxList[bin_size].append(barcode_dict[barcode.decode()].id);
        
        barcode_id += 1;
        del bin_dict, bins, frag_list_uniq;
    
    dt = h5py.special_dtype(vlen=bytes)    
    f.create_dataset("AM/nBinSize", data=len(bin_dict_list),  dtype="uint32");
    f.create_dataset("AM/binSizeList", data=list(bin_dict_list.keys()),  dtype="uint32");
    
    for bin_size in bin_dict_list:
        f.create_dataset("AM/"+str(bin_size)+"/binChrom",data=[np.string_(key[0]) for key in bin_dict_list[bin_size]], dtype=h5py.special_dtype(vlen=bytes), compression="gzip", compression_opts=9);
        f.create_dataset("AM/"+str(bin_size)+"/binStart",data=[key[1] for key in bin_dict_list[bin_size]], dtype="uint32", compression="gzip", compression_opts=9);
        f.create_dataset("AM/"+str(bin_size)+"/idx", data=idxList[bin_size], dtype="uint32", compression="gzip", compression_opts=9);
        f.create_dataset("AM/"+str(bin_size)+"/idy", data=idyList[bin_size], dtype="uint32", compression="gzip", compression_opts=9);
        f.create_dataset("AM/"+str(bin_size)+"/count", data=countList[bin_size], dtype="uint8", compression="gzip", compression_opts=9);    
    f.close()  
    return 0
