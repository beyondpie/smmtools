"""
Graph-based clustering: Leiden algorithm.

Ref: Kai's and Yang's codes
"""

from time import perf_counter as pc
import leidenalg as la
import igraph as ig

def leiden(knn, reso: float = 1.0, seed: int = 10, opt: str = "RB"):
    """Run and count to peak"""
    start_time = pc()
    """ init input files """
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
