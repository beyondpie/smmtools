"""
Graph-based clustering: Leiden algorithm.

Ref: Kai's and Yang's codes
"""

from time import perf_counter as pc
import leidenalg as la
import igraph as ig
from scipy.io import mmread
from scipy import sparse
from time import perf_counter as pc
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse import save_npz
import numpy as np
import itertools
from scipy.cluster.hierarchy import linkage, leaves_list, cophenet
from scipy.spatial.distance import squareform

import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig, imshow, set_cmap
import fastcluster as fc

plt.switch_backend("agg")


def leiden(knn, reso: float = 1.0, seed: int = None, opt: str = "RB"):
    """Run and count to peak

    We validate this function with Yang's script.

    Args:
        knn (sparse matrix): scipy.sparse.csc_matrix format when using
          R reticult::r_to_py.
        reso (float): resolution for clustering, (0,1],
        seed (int): used in clustering, default is None (the same as in igraph),
        opt (str): "RB" or "CPM", default is RB
          RB short for RBConfigurationVertexPartition
          CPM short for CPMVertexPartition

    Returns:
        integer vector, cluster idx (start from 0) for the nodes in order.
        the index are ordered based on the sizes of the clusters. (0 -> max size)
    """
    start_time = pc()
    vcount = max(knn.shape)
    sources, targets = knn.nonzero()
    edgelist = list(zip(sources.tolist(), targets.tolist()))
    g = ig.Graph(vcount, edgelist)
    if opt == "RB":
        partition = la.find_partition(
            graph=g,
            partition_type=la.RBConfigurationVertexPartition,
            resolution_parameter=reso,
            seed=seed,
        )
    else:
        partition = la.find_partition(
            graph=g,
            partition_type=la.CPMVertexPartition,
            resolution_parameter=reso,
            seed=seed,
        )
    part_membership = partition.membership
    end_time = pc()
    print(
        f"leiden uses (secs) with partition type {opt}: {round(end_time - start_time, 3)}"
    )
    return part_membership


## From Yang's script
def cal_connectivity(P):
    """calculate connectivity matrix"""
    # print("=== calculate connectivity matrix ===")
    # connectivity_mat = lil_matrix((len(P), len(P)))
    # connectivity_mat = csr_matrix((len(P), len(P)))
    connectivity_mat = np.zeros((len(P), len(P)))
    classN = max(P)
    for cls in range(classN + 1):
        xidx = [i for i, x in enumerate(P) if x == cls]
        iterables = [xidx, xidx]
        for t in itertools.product(*iterables):
            connectivity_mat[t[0], t[1]] = 1
    """connectivity_mat = csr_matrix(connectivity_mat)"""
    return connectivity_mat


def plotC(prefix, C):
    """
    Plot reordered consensus matrix.

    :param C: Reordered consensus matrix.
    :type C: numpy.ndarray`
    :param rank: Factorization rank.
    :type rank: `int`
    """
    fig = plt.figure(figsize=(5, 5), dpi=100)
    imshow(C)
    set_cmap("RdBu_r")
    fileN = [prefix, "C", "png"]
    fileN = ".".join(fileN)
    fig.savefig(fileN)


def reorder(C):
    """
    Reorder consensus matrix.

    :param C: Consensus matrix.
    :type C: `numpy.ndarray`
    """
    Y = 1 - C
    Z = linkage(squareform(Y), method="average")
    ivl = leaves_list(Z)
    ivl = ivl[::-1]
    return C[:, ivl][ivl, :]


def plot_CDF(prefix, C, u1, u2, num_bins=100):
    counts, bin_edges = np.histogram(C, bins=num_bins, density=True)
    cdf = np.cumsum(counts)
    fig = plt.figure(dpi=100)
    plt.plot(bin_edges[1:], cdf)
    plt.xlabel("Consensus index value")
    plt.ylabel("CDF")
    plt.axvline(x=u1, color="grey", linestyle="--")
    plt.axvline(x=u2, color="grey", linestyle="--")
    fileN = [prefix, "cdf", "png"]
    fileN = ".".join(fileN)
    fig.savefig(fileN)
    outBinEdges = ".".join([prefix, "cdf.txt"])
    with open(outBinEdges, "w") as fo:
        fo.write("\t".join(str(i) for i in cdf) + "\n")


def cumfreq(a, numbins=100, defaultreallimits=None):
    # docstring omitted
    h, l, b, e = histogram(a, numbins, defaultreallimits)
    cumhist = np.cumsum(h * 1, axis=0)
    return cumhist, l, b, e


def cal_cophenetic(C):
    """calculate cophenetic correlation coefficient"""
    print("=== calculate cophenetic correlation coefficient ===")
    X = C
    """Z = linkage(X)"""
    Z = fc.linkage_vector(X)  # Clustering
    orign_dists = fc.pdist(X)  # Matrix of original distances between observations
    cophe_dists = cophenet(Z)  # Matrix of cophenetic distances between observations
    corr_coef = np.corrcoef(orign_dists, cophe_dists)[0, 1]
    return corr_coef


def cal_dispersion(C):
    """calculate dispersion coefficient"""
    print("=== calculate dispersion coefficient ===")
    start_t = pc()
    n = C.shape[1]
    corr_disp = np.sum(4 * np.square(np.concatenate(C - 1 / 2))) / (np.square(n))
    end_t = pc()
    print("Used (secs): ", end_t - start_t)
    return corr_disp


def cal_PAC(C, u1, u2):
    """calculate PAC (proportion of ambiguous clustering)"""
    print("=== calculate PAC (proportion of ambiguous clustering) ===")
    start_t = pc()
    totalN = C.shape[0] * C.shape[0]
    u1_fraction = (C.ravel() < u1).sum() / totalN
    u2_fraction = (C.ravel() < u2).sum() / totalN
    PAC = u2_fraction - u1_fraction
    end_t = pc()
    print("Used (secs): ", end_t - start_t)
    return PAC


def cal_stab(x):
    """calculate stability for every cell"""
    s = np.sum(abs(x - 0.5)) / (0.5 * x.shape[0])
    return s


def run(
    knn, reso: float, N: int, outf: str, u1: float = 0.05, u2: float = 0.95
) -> None:
    """Run and count to peak

    Args:
        knn (sparse Matrix)
        reso (float): resolution
        N (int): iteration of N time
        out (str): output file name prefix, incluses the directory
        u1 (float): 0.05, left interval cutoff of CDF
        u2 (float): 0.95, right interval cutoff of CDF

    Side effects:
        - saved consensus matrix into '.'.join([outf, "consensus", "npz"])
        - output the figure of consensus matrix
        - ouput file of stat into '.'.join([outf, "stat.txt"])
        - output file of stab into '.'.join([outf, "stab.txt"])

    Returns:
        None
    """
    start_time = pc()
    # inf = args.input
    # reso = args.resolution
    # N = args.N
    # sampleN = args.sample
    # u1 = args.u1
    # u2 = args.u2
    # outf = args.output

    # knn = mmread(inf)
    # knn = knn.tocsr()
    dimN = knn.shape[0]

    vcount = max(knn.shape)
    sources, targets = knn.nonzero()
    edgelist = list(zip(sources.tolist(), targets.tolist()))
    g = ig.Graph(vcount, edgelist)

    # generate consensus matrix
    dims = knn.shape[0]
    if dims > sampleN:
        dims = sampleN
        import random

        random.seed(2022)
        idxy = random.sample(range(dimN), sampleN)
        idxy = sorted(idxy)
    consensus = csr_matrix((dims, dims))

    print("=== calculate connectivity matrix ===")
    for seed in range(N):
        start_t = pc()
        partition = la.find_partition(
            g, la.RBConfigurationVertexPartition, resolution_parameter=reso, seed=seed
        )
        part_membership = partition.membership
        part_membership = np.array(part_membership)
        if len(part_membership) > sampleN:
            part_membership = part_membership[idxy]  # downsample (10,000 observations)
        outs = cal_connectivity(part_membership)
        consensus += outs
        end_t = pc()
        print("seed ", seed, " used (secs): ", end_t - start_t)
    consensus /= N

    # save consensus matrix
    consensus_sp = sparse.csr_matrix(consensus)
    outfname = ".".join([outf, "consensus", "npz"])
    save_npz(outfname, consensus_sp)

    # plotting
    order_consensus = reorder(consensus)
    plotC(outf, order_consensus)
    plot_CDF(outf, consensus, u1, u2)

    # cal measurement
    # o_cophcor = cal_cophenetic(consensus)
    o_disp = cal_dispersion(consensus)
    o_PAC = cal_PAC(consensus, u1, u2)
    print("=== calculate stability for every cell ===")
    o_stab = np.apply_along_axis(cal_stab, 0, consensus)

    # write stat
    out_list = [outf, reso, o_disp, o_PAC]
    outstat = ".".join([outf, "stat.txt"])
    with open(outstat, "w") as fo:
        fo.write("\t".join(str(i) for i in out_list) + "\n")

    outStab = ".".join([outf, "stab.txt"])
    with open(outStab, "w") as fo:
        fo.write("\t".join(str(i) for i in o_stab) + "\n")

    end_time = pc()
    print("Total used (secs): ", end_time - start_time)
