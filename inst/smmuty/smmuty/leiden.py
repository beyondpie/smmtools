"""
Graph-based clustering: Leiden algorithm.

Ref: Kai's and Yang's codes
"""

from time import perf_counter as pc
import leidenalg as la
import igraph as ig

def leiden(knn, reso: float = 1.0, seed: int = None, opt: str = "RB"):
    """Run and count to peak

    Args:
        knn (sparse matrix): scipy.sparse.csc_matrix format when using
          R reticult::r_to_py.
        reso (float): resolution for clustering, (0,1],
        seed (int): used in clustering, default is None (the same as in igraph),
        opt (str): "RB" or "CPM", default is RB
          RB short for RBConfigurationVertexPartition
          CPM short for CPMVertexPartition
    
    Returns:
        two-column matrix: [node, cluster], all the index start from zero
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
