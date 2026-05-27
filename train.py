from __future__ import annotations

import json
import random
from dataclasses import dataclass, asdict
from typing import Optional, Dict
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.metrics import silhouette_score as sklearn_silhouette_score
from sklearn.metrics.cluster import adjusted_rand_score

import anndata
import scanpy as sc

#HYPERPARAMETERS
@dataclass
class AutoResearchConfig:
    # --QC--
    min_genes: int = 1
    max_genes: int = 100000
    min_counts: int = 10
    max_counts: int = 100000
    max_pct_mt: float = 20.0

    # --Normalization--
    target_sum: int = 1000

    # --Highly Variable Genes--
    n_top_genes: int = 2000

    # --Scaling--
    do_scale: bool = True
    max_value: Optional[float] = 10.0

    # --PCA--
    n_pcs: int = 30
    svd_solver: str = "arpack"

    # --Neighborhood Graph--
    n_neighbors: int = 50
    metric: str = "euclidean"

    # --Clustering--
    resolution: float = 0.6

    #--Misc--
    random_seed: int = 42

#UTILS
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)

KNOWN_MARKERS = {
    "T_cell": ["CD3D", "TRBC1", "IL7R"],
    "B_cell": ["MS4A1", "CD79A"],
    "Macrophage": ["LST1", "CTSS", "FCER1G"],
    "Epithelial": ["EPCAM", "KRT18", "KRT19"],
}

#==================================
#           TRAINING PIPELINE
#==================================

#DATA LOADING
def load_data(file_path: str) -> anndata.AnnData:
    adata = sc.read_h5ad(file_path)

    adata.var_names = adata.var['feature_name']
    adata.var_names_make_unique()
    return adata

#QC
def compute_qc_metrics(adata: anndata.AnnData):
    adata.var['mt'] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    return adata

def apply_qc(adata: anndata.AnnData, config: AutoResearchConfig):
    adata = adata[
        (adata.obs.n_genes_by_counts >= config.min_genes) &
        (adata.obs.n_genes_by_counts <= config.max_genes) &
        (adata.obs.total_counts >= config.min_counts) &
        (adata.obs.total_counts <= config.max_counts) &
        (adata.obs.pct_counts_mt <= config.max_pct_mt)
    ].copy()
    return adata

#NORMALIZATION AND LOG TRANSFORMATION
def normalize_data(adata: anndata.AnnData, config: AutoResearchConfig):
    sc.pp.normalize_total(adata, target_sum=config.target_sum)
    sc.pp.log1p(adata)
    return adata

#HIGHLY VARIABLE GENES
def select_hvg(adata: anndata.AnnData, config: AutoResearchConfig):
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=config.n_top_genes,
        flavor="seurat",
        batch_key=None
    )
    adata = adata[:, adata.var.highly_variable].copy()
    return adata

#SCALING
def scale_data(adata: anndata.AnnData, config: AutoResearchConfig):
    if config.do_scale:
        sc.pp.scale(adata, max_value=config.max_value)
    return adata

#PCA
def run_pca(adata: anndata.AnnData, config: AutoResearchConfig):
    sc.tl.pca(
        adata, 
        n_comps=config.n_pcs, 
        svd_solver=config.svd_solver,
        random_state=config.random_seed)
    return adata

#NEIGHBORHOOD GRAPH
def compute_neighbors(adata: anndata.AnnData, config: AutoResearchConfig):
    sc.pp.neighbors(
        adata, 
        n_neighbors=config.n_neighbors, 
        n_pcs=config.n_pcs,
        metric=config.metric,
        random_state=config.random_seed)
    return adata

#CLUSTERING
def cluster_cells(adata: anndata.AnnData, config: AutoResearchConfig):
    sc.tl.leiden(
        adata, 
        resolution=config.resolution, 
        key_added="cluster",
        random_state=config.random_seed,
        flavor="igraph",
        n_iterations=2,
        directed=False)
    return adata


#==============================
#           EVALUATION
#==============================
#MARKER COHERENCE SCORE
def compute_marker_genes(adata: anndata.AnnData):
    sc.tl.rank_genes_groups(adata, groupby="cluster", method="wilcoxon", use_raw=True)
    return adata

def marker_coherence_score(adata):
    if "rank_genes_groups" not in adata.uns:
        return 0.0

    marker_df = pd.DataFrame(adata.uns["rank_genes_groups"]["names"])

    cluster_scores = []

    for cluster in marker_df.columns:
        top_genes = set(marker_df[cluster].head(20).astype(str))

        best_match = 0.0

        for celltype, markers in KNOWN_MARKERS.items():
            overlap = len(top_genes & set(markers))
            score = overlap / len(markers)

            best_match = max(best_match, score)

        cluster_scores.append(best_match)

    return float(np.mean(cluster_scores))

#STABILITY SCORE
def run_subsample_clustering(adata: anndata.AnnData, config: AutoResearchConfig):
    n = adata.n_obs
    subset_indices = np.random.choice(n, size=int(0.8 * n), replace=False)
    subset_adata = adata[subset_indices].copy()
    subset_adata = compute_qc_metrics(subset_adata)
    subset_adata = apply_qc(subset_adata, config)
    subset_adata = normalize_data(subset_adata, config)
    subset_adata = select_hvg(subset_adata, config)
    subset_adata = scale_data(subset_adata, config)
    subset_adata = run_pca(subset_adata, config)
    subset_adata = compute_neighbors(subset_adata, config)
    subset_adata = cluster_cells(subset_adata, config)
    return subset_adata

def cluster_stability_score(adata: anndata.AnnData, adata_raw: anndata.AnnData, config: AutoResearchConfig) -> float:
    scores = []
    base = adata.obs["cluster"].astype(str)
    for _ in range(3):
        sub = run_subsample_clustering(adata_raw, config)
        overlap = np.intersect1d(adata.obs_names, sub.obs_names)
        if len(overlap) < 10:
            continue

        base_labels = base[overlap]
        sub_labels = sub.obs["cluster"].loc[overlap].astype(str)

        ari = adjusted_rand_score(base_labels, sub_labels)
        scores.append(ari)
    
    if not scores:
        return 0.0
    return float(np.mean(scores))

#SILHOUETTE SCORE
def compute_silhouette_score(adata: anndata.AnnData) -> float:
    X = adata.obsm["X_pca"]
    labels = adata.obs["cluster"].astype(str)

    if len(np.unique(labels)) < 2:
        return -1.0
    
    score = sklearn_silhouette_score(X, labels)
    score = (score + 1) / 2
    return float(score)


#CLUSTER BALANCE SCORE
def cluster_balance_score(adata: anndata.AnnData) -> float:
    counts = adata.obs["cluster"].value_counts().values
    proportions = counts / counts.sum()
    entropy = -np.sum(proportions * np.log(proportions + 1e-8))
    max_entropy = np.log(len(counts))
    if max_entropy == 0:
        return 0.0
    return float(entropy / max_entropy)

#OBJECTIVE SCORE
def compute_objective(scores:Dict[str, float]) -> float:
    objective = (
        0.4 * scores["marker_coherence"] +
        0.3 * scores["cluster_stability"] +
        0.2 * scores["silhouette"] +
        0.1 * scores["cluster_balance"]
    )
    return float(objective)

#VISUALIZATION
def save_umap(adata):
    sc.tl.umap(adata)
    sc.pl.umap(adata, color="cluster", show=False)
    plt.savefig(OUTPUT_DIR / "umap.png", dpi=200)
    plt.close()

#RESULTS LOGGING
def save_scores(scores: Dict[str, float]):
    with open(OUTPUT_DIR / "scores.json", "w") as f:
        json.dump(scores, f, indent=2)

def append_results(config: AutoResearchConfig, scores: Dict[str, float]):
    row = {
        **asdict(config),
        **scores,
    }

    df = pd.DataFrame([row])
    path = OUTPUT_DIR / "results.csv"
    if path.exists():
        df.to_csv(path, mode="a", header=False, index=False)
    else:
        df.to_csv(path, index=False)

def write_summary(scores: Dict[str, float], cfg: AutoResearchConfig):
    lines = []
    lines.append("=== AutoResearch Configuration ===")
    lines.append(json.dumps(asdict(cfg), indent=2))

    lines.append("\n=== Evaluation Scores ===")
    lines.append(json.dumps(scores, indent=2))

    with open(OUTPUT_DIR / "summary.txt", "w") as f:
        f.write("\n".join(lines))

#MAIN TRAINING FUNCTION
def main():
    config = AutoResearchConfig()
    adata_raw = load_data("data/lims_lung_celltype_demo.h5ad")
    adata = adata_raw.copy()
    adata = compute_qc_metrics(adata)
    adata = apply_qc(adata, config)
    adata = normalize_data(adata, config)
    adata.raw = adata.copy()
    adata = select_hvg(adata, config)
    adata = scale_data(adata, config)
    adata = run_pca(adata, config)
    adata = compute_neighbors(adata, config)
    adata = cluster_cells(adata, config)

    adata = compute_marker_genes(adata)

    scores = {
        "marker_coherence": marker_coherence_score(adata),
        "cluster_stability": cluster_stability_score(adata, adata_raw, config),
        "silhouette": compute_silhouette_score(adata),
        "cluster_balance": cluster_balance_score(adata),
    }
    scores["objective"] = compute_objective(scores)

    save_scores(scores)
    append_results(config, scores)
    write_summary(scores, config)

    print(json.dumps(scores, indent=2))

#ENTRY POINT
if __name__ == "__main__":
    main()