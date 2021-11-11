import sys
import os
import tempfile
import time
import subprocess
import collections
from typing import Dict, Tuple, List
import pysam

class QC(object):
    """A quality control object that has the following attributes:

    Attributes:
        barcode: name of barcodes
        total: total number of sequenced fragments.
        mapped: number of mappable fragments.
        chrM: number of fragments mapped to chrM.
        paired: number of fragments are paired.
        single: number of fragments are from single read.
        proper_paired: number of paired reads are properly paired.
        usable: number of usable fragments.
        uniq: number of unique fragments.
        isize: average insert size distribution.
    """
    fields = ["barcode", "total", "mapped", "single", "secondary","paired",
                       "proper_paired","proper_flen", "usable", "uniq", "chrM"]

    def __init__(self) -> None:
        """Return a qc object"""
        self.barcode:str = ""
        self.total = 0
        self.mapped = 0
        self.single = 0
        self.secondary = 0
        self.paired = 0
        self.proper_paired = 0
        self.proper_flen = 0
        self.usable = 0
        self.uniq = 0
        self.chrM = 0
    def to_str(self, sep = "\t") -> str:
        return sep.join([str(getattr(self, i)) for i in self.fields])
        
        
        

class Fragment(object):
    """A fragment object that has the following attributes:

    Attributes:
        chrom: chromsome name
        start: start position
        end: end position
        mapq: mapping quality
        is_proper_pair: whether properly paired
        is_single: whether it is a single read
        is_secondary: whether it is a secondary alignment read
        flen: fragment length
    """
    def __init__(self, qname, chrom, pos, flen, mapq, is_single, is_secondary, is_proper_pair):
        """Return a qc object"""
        self.qname:str = qname
        self.chrom:str = chrom
        self.pos:int = pos
        self.flen:int = flen
        self.mapq:float = mapq
        self.is_single: bool = is_single
        self.is_secondary: bool = is_secondary
        self.is_proper_pair:bool = is_proper_pair


def getBarcodesFromBam(input_bam: str) -> Dict[str, float]:
    """Identify unique barcodes from the bam file
    
    Args:
        input_bam: a bam file

    Returns:
        A dictionary contains all barcodes, their coverages.
    """
    
    barcode_dict = collections.defaultdict(lambda :0.0)
    samfile = pysam.AlignmentFile(input_bam, "rb")
    for _read in samfile:
        barcode = _read.qname.split(":")[0].upper()
        # approximate counting, a read is half fragment
        if barcode not in barcode_dict:
            barcode_dict[barcode] = 0.5
        else:
            barcode_dict[barcode] += 0.5
    samfile.close()
    return barcode_dict


def group_reads_by_barcode_bam(input_bam):
    """ Group reads based on the barcodes from SnapTools.
    
    Args:
        input_bam: a bam file

    Returns:
        Generator that contains reads sharing the same barcode
    """
    read_group_list = []; 
    pre_barcode = "";
    samfile = pysam.AlignmentFile(input_bam, "rb");
    for cur_read in samfile:
        cur_barcode = cur_read.qname.split(":")[0];
        if cur_barcode == pre_barcode:
            read_group_list.append(cur_read)
        else:
            if pre_barcode != "":
                # return read group
                yield (x for x in read_group_list)
            read_group_list = [cur_read] # add the first read
            pre_barcode = cur_barcode
    # reads from the last barcode
    yield (x for x in read_group_list)

def pairReadsByName(read_list):
    """ Pair reads based on read names from SnapTools.
    
    Args:
        read_list: a list of reads that share the same barcode

    Returns:
        Generator contains read pairs from the same fragment
        and a boolen variable indicates whether it is supplementary alignment
    """
    # pair up 
    for read1 in read_list:
        # read until read1 is not a supplementary alignment
        while(read1.is_supplementary):
            yield (read1, None, False, True)
            try:
                #print "Warning: skipping ", read1.qname;
                read1 = next(read_list);
            except:
            	break        
        try:
            read2 = next(read_list);
        except:
        	break
        while(read2.is_supplementary):
            yield (read2, None, False, True)
            try:
                #print "Warning: skipping ", read2.qname;
                read2 = next(read_list);
            except:
            	break
        if(read1.qname != read2.qname):
            while (read1.qname != read2.qname):
                yield(read1, None, False, False);
                read1 = read2;
                try:
                    read2 = next(read_list);
                    while(read2.is_supplementary):
                        try:
                            #print "Warning: skipping ", read2.qname;
                            read2 = next(read_list);
                        except:
                        	break
                except:
                    break;
        yield (read1, read2, True, False)

def readToFragment(read1, is_secondary):
    """ convert a single read to fragment
    
    Args:
        read1: a single read

    Returns:
        Generator contains a fragment object
    """
    try:
        read1.qname;
    except ValueError as e:
        sys.exit('readto_fragment: can not get read1 name!');   

    barcode = read1.qname.split(":")[0];
    mapq = read1.mapq;
    try:
        chrom1 = read1.reference_name;
        start1 = read1.reference_start;
        strand1 = "-" if read1.is_reverse else "+";    
        # it is possible that flen1 is None  
        flen1 = read1.reference_length if read1.reference_length != None else 0
        end1   = start1 + flen1;
        start = min(start1, end1);
        end = max(start1, end1);
        #(qname, chrom, pos, flen, mapq, is_single, is_secondary, is_proper_pair):
        return Fragment(read1.qname, chrom1, start, abs(start - end), mapq, True, is_secondary, False);
    except ValueError as e:
        return Fragment(read1.qname, None, None, None, mapq, True, is_secondary, False);

def readPairToFragment(read1, read2, is_secondary):
    """ convert read pairs to fragments
    
    Args:
        read1: R1 read

        read2: R2 read
    Returns:
        Generator that contains a fragment object
    """
    try:
        read1.qname;
        read2.qname;
        if read1.qname != read2.qname:
            sys.exit('read_pair_to_fragment: read1 and read2 name does not match!');  
    except ValueError as e:
        sys.exit('read_pair_to_fragment: can not get read1 or read2 name!');   
    barcode = read1.qname.split(":")[0];
    mapq = min(read1.mapq, read2.mapq);
    try:
        chrom1 = read1.reference_name;
        start1 = read1.reference_start;
        strand1 = "-" if read1.is_reverse else "+";    
        chrom2 = read2.reference_name;
        start2 = read2.reference_start;
        strand2 = "-" if read1.is_reverse else "+";    
        # it is possible that flen1 is None  
        flen1 = read1.reference_length if read1.reference_length != None else 0
        flen2 = read2.reference_length if read2.reference_length != None else 0
        end1   = start1 + flen1;
        end2   = start2 + flen2;
        start = min(start1, end1, start2, end2);
        end = max(start1, end1, start2, end2);
        return Fragment(read1.qname, chrom1, start, abs(start - end), mapq, False, is_secondary, read1.is_proper_pair);
    except ValueError:
        return Fragment(read1.qname, None, None, None, mapq, False, is_secondary, False);


def SnapToolsBamTo10xFragmentBed(bam_file: str, outf:str,
                                 qc_file: str,
                                 min_mapq: int = 30,
                                 min_flen: int = 0,
                                 max_flen: int = 1000,
                                 min_cov: int = 100,
                                 num_threads: int = 1,
                                 verbose: bool = True) -> None:
    """Get 10x Fragment bed file format from bam file.

    Bam file is handled by the SnapTools' logic.
    NOTE: need ~bedtools~ in the env.
    """
    outdir = os.path.dirname(outf)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(os.path.dirname(qc_file)):
        os.makedirs(os.path.dirname(qc_file))
    ## make sure bam files are sorted by barcode name
    ftmp = tempfile.NamedTemporaryFile(delete = False, dir = outdir)
    ftmp.close()
    pysam.sort("-n", "-@", str(num_threads), "-o", ftmp.name, bam_file)

    check_point:int = 0
    start_time = time.time()

    barcode_cov = getBarcodesFromBam(input_bam=bam_file)
    barcodes = {k: v for k, v in barcode_cov.items() if v >= min_cov}

    input_bam = pysam.AlignmentFile(ftmp.name, "rb")
    output_bed = open(outf, "wa")
    output_qc = open(qc_file, "wa")
    output_qc.write(f"{'\t'.join(QC.fields)}\n")
    for read_group in group_reads_by_barcode_bam(input_bam=bam_file):
        ##  reads from read_group come from the same barcode
        frag_list: List[Tuple[str, int, int]] = []
        qc = QC()
        for (read1, read2, is_paired, is_secondary) in pairReadsByName(read_list=read_group):
            if is_paired:
                frag = readPairToFragment(read1, read2, is_secondary)
            # supplementary alignments or unmated reads
            else:
                frag = readToFragment(read1, is_secondary);            
            # extract the barcode
            barcode = frag.qname.split(":")[0].upper()
            qc.barcode = barcode
            # filter by min_cov
            if not barcode in barcodes:
                break
            # only for printing the progress
            check_point += 1;
            if verbose and check_point%100000 == 0:
                print(f"{check_point} tags, {time.time() - start_time} seconds.")
            # total number of sequencing fragments (exclude supplementary alignments)
            if frag.is_secondary == False:
                qc.total += 1
            ## 1. Filter non-unique mapped fragments
            if frag.mapq < min_mapq:
                continue
            if frag.is_single:
                if frag.is_secondary:
                    qc.secondary += 1
                else:
                    qc.single += 1
                    ## not keep the single reads
                    continue
            else:
                ## pared-end reads
                qc.paired += 1
                if frag.is_proper_pair:
                    qc.proper_paired += 1
                else:
                    continue
            # 2. check fragment size
            if frag.flen > min_flen and frag.flen < max_flen:
                qc.proper_flen += 1
            else:
                continue
            # 3. combine single and paired as fragments
            frag_list.append((frag.chrom, frag.pos, frag.pos+frag.flen))
        if(len(frag_list) < 1):
            continue
        frag_counter = collections.Counter(frag_list)
        qc.usable = len(frag_list)
        qc.uniq = len(frag_counter)
        qc.chrM = sum([ e[0] == "chrM" for e in frag_list])
        n_uniq_chrM = sum([ e[0] == "chrM" for e in frag_counter.keys()])
        if (n_uniq_chrM == qc.uniq):
            continue
        ## output the fragments
        frag_content: List[str] = [f"{k[0]}\t{k[1]}\t{k[2]}\t{v}" for k, v in frag_counter.items()]
        output_bed.write("\n".join(frag_content))
        output_qc.write(f"{qc.to_str()}\n")
    input_bam.close() 
    subprocess.check_call(["rm", ftmp.name])
    output_bed.close()
    output_qc.close()
    return None
