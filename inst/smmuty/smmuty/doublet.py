"""Use Scrublet to remove doublets within samples

Ref: yang Li and Kai Zhang's codes
"""
import numpy as np
import scrublet as scr
import numpy as np
from sklearn.mixture import BayesianGaussianMixture


def detectDoublet(
    counts_matrix,
    expected_doublet_rate: float = 0.06,
    min_counts: int = 3,
    min_cells: int = 3,
    min_gene_variability_pctl: float = 85,
    n_prin_comps: int = 30,
):
    """
    counts_matrix: could be sparse matrix
    """
    scrub = scr.Scrublet(
        counts_matrix,
        expected_doublet_rate=expected_doublet_rate,
        sim_doublet_ratio=3,
        n_neighbors=25,
    )
    doublet_scores, _ = scrub.scrub_doublets(
        min_counts=min_counts,
        min_cells=min_cells,
        min_gene_variability_pctl=min_gene_variability_pctl,
        mean_center=True,
        normalize_variance=True,
        # n_prin_comps=min(30, counts_matrix.get_shape()[0] // 10),
        n_prin_comps=n_prin_comps,
    )

    # Fit a Gaussian mixture model
    X = scrub.doublet_scores_sim_
    X = np.array([X]).T
    gmm = BayesianGaussianMixture(n_components=2, max_iter=1000, random_state=2394).fit(
        X
    )
    i = np.argmax(gmm.means_)

    probs_sim = gmm.predict_proba(X)[:, i]
    vals = X[np.argwhere(probs_sim > 0.5)].flatten()
    if vals.size == 0:
        threshold = np.amax(X.flatten())
    else:
        threshold = min(vals)

    X = np.array([doublet_scores]).T
    probs = gmm.predict_proba(X)[:, i].tolist()
    return (
        threshold,
        scrub.threshold_,
        probs,
        doublet_scores,
        scrub.doublet_scores_sim_,
    )
