import pybedtools
import collections

gene_file = "../data/gencode.vM16.geneUp2k.bed"
gene_list= set([str(item.name) for item in pybedtools.BedTool(gene_file)])

gene_dict = collections.OrderedDict(list(zip(gene_list, list(range(1, len(gene_list) + 1)))));    

