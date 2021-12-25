import pybedtools
import collections
import tempfile
import snaptools
from snaptools.snap import getBarcodesFromSnap
from snaptools.snap import dump_read

gene_file = "../data/gencode.vM16.geneUp2k.binarystr.bed"
gene_list= set([str(item.name) for item in pybedtools.BedTool(gene_file)])
gene_dict = collections.OrderedDict(list(zip(gene_list, list(range(1, len(gene_list) + 1)))));    

sample = "CEMBA210429_18B"
snap_file = f"snap/{sample}.snap"
buffer_size = 1000
tmp_folder = None


# extract the barcodes
barcode_dict = getBarcodesFromSnap(snap_file)

fout_frag = tempfile.NamedTemporaryFile(delete=False, dir=tmp_folder)
dump_read(snap_file, fout_frag.name, buffer_size, None, tmp_folder, True)

# in parallel find the overlap cell and peaks
frag_bt = pybedtools.BedTool(fout_frag.name)
peak_bt = pybedtools.BedTool(gene_file)

# count for frequency
cell_peak_arr = collections.defaultdict(list)
items = frag_bt.intersect(peak_bt, wa = True, wb = True)

for item in frag_bt.intersect(peak_bt, wa=True, wb=True):
    key = str(item.fields[7])
    if key in gene_dict:
        idy = gene_dict[key]
        barcode = item.name.split(":")[0]
        idx = barcode_dict[barcode].id
        cell_peak_arr[idx].append(idy)

IDX_LIST = []
IDY_LIST = []
VAL_LIST = []
for barcode_id in cell_peak_arr:
    d = collections.Counter(cell_peak_arr[barcode_id])
    IDX_LIST += [barcode_id] * len(d)
    for peak_id in d:
        IDY_LIST.append(peak_id)
        VAL_LIST.append(d[peak_id])

