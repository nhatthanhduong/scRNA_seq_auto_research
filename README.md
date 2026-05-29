# AutoResearch: scRNA-seq Hyperparameter Tuning

An autonomous hyperparameter search framework for single-cell RNA-seq (scRNA-seq) clustering pipelines. This project implements a loop that iteratively tunes preprocessing and clustering parameters to maximize a composite objective score, automatically committing successful configurations and reverting regressions.

## Overview

The pipeline processes 10x scRNA-seq data through a standard Scanpy workflow — QC, normalization, highly variable gene selection, scaling, PCA, neighborhood graph construction, and Leiden clustering — then evaluates the clustering quality using four complementary metrics. The experiment loop operates autonomously: tune parameters, commit, run, evaluate, keep or discard, and repeat indefinitely.

**Objective function** (range ≈[-0.2, 1.0], higher is better):

```
objective = 0.4 × marker_coherence + 0.3 × cluster_stability
           + 0.2 × silhouette            + 0.1 × cluster_balance
```

## Project Structure

```
├── train.py              # Training pipeline + AutoResearchConfig (the only file you edit)
├── prepare.py            # Dataset download and preparation from cellxgene-census
├── program.md            # Experiment loop instructions and protocol
├── visualization.py      # Objective improvement plotting over runs
├── pyproject.toml        # Project configuration and dependencies
├── .gitignore            # Git ignore rules
├── data/
│   └── lims_lung_celltype_demo.h5ad   # Input AnnData dataset
├── output/
│   ├── results.csv       # Experiment results (hyperparameters + scores)
│   ├── umap.png          # UMAP visualization of the best clustering
│   └── objective_improvement.png  # Objective score progression plot
└── README.md             # This file
```

## Setup

### Prerequisites

- Python ≥3.10, <3.13
- [uv](https://docs.astral.sh/uv/) package manager

### Installation

```bash
# Clone the repository
git clone <repository-url>
cd lims-single-cell

# Create branch for experiments
git checkout -b autoresearch/may26

# Install dependencies
uv sync

# Prepare dataset (if not already present)
uv run prepare.py
```

The `prepare.py` script downloads a human lung scRNA-seq dataset from CZ CELLxGENE Census, filtering for normal and COVID-19 samples and sampling up to 500 cells per cell type. Resulting cell types include T cells, B cells, natural killer cells, macrophages, monocytes, dendritic cells, epithelial cells, endothelial cells, and fibroblasts.

## Pipeline

The pipeline follows a standard Scanpy workflow:

| Step | Function | Description |
|------|----------|-------------|
| **Load** | `load_data()` | Read the h5ad file, set gene names as `var_names` |
| **QC** | `compute_qc_metrics()` + `apply_qc()` | Calculate QC metrics, filter cells by gene/count/mito thresholds |
| **Normalize** | `normalize_data()` | Total-count normalize to `target_sum`, then log1p transform |
| **HVG** | `select_hvg()` | Select highly variable genes (Seurat flavor) |
| **Scale** | `scale_data()` | Z-score scaling with optional `max_value` clipping |
| **PCA** | `run_pca()` | Dimensionality reduction to `n_pcs` components |
| **Neighbors** | `compute_neighbors()` | Compute k-nearest neighbor graph |
| **Cluster** | `cluster_cells()` | Leiden community detection |
| **Evaluate** | Various | Compute marker coherence, cluster stability, silhouette, and balance scores |

## Hyperparameter Tuning

Edit `AutoResearchConfig` in `train.py` (lines 20-52) to tune parameters. Only this dataclass should be modified — do not change pipeline logic, `prepare.py`, or `pyproject.toml`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_genes` | int | 100 | Minimum genes per cell |
| `max_genes` | int | 100000 | Maximum genes per cell |
| `min_counts` | int | 10 | Minimum UMI counts per cell |
| `max_counts` | int | 100000 | Maximum UMI counts per cell |
| `max_pct_mt` | float | 10.0 | Max % mitochondrial reads |
| `target_sum` | int | 1000 | Normalization target sum |
| `n_top_genes` | int | 4000 | Number of highly variable genes |
| `n_pcs` | int | 60 | Number of PCA components |
| `svd_solver` | str | `"randomized"` | SVD solver (`"arpack"` / `"randomized"`) |
| `n_neighbors` | int | 45 | Number of neighbors for graph |
| `metric` | str | `"euclidean"` | Distance metric (`"cosine"` / `"euclidean"`) |
| `resolution` | float | 0.8 | Leiden clustering resolution |
| `random_seed` | int | 42 | Random seed for reproducibility |

**Do not touch**: `flavor`, `batch_key`, `do_scale`, `max_value`

## Experiment Loop

The experiment runs an autonomous optimization loop:

1. **Check state** — note the current git branch/commit
2. **Tune** — modify 1-3 `AutoResearchConfig` parameter defaults
3. **Commit** — `git add -A && git commit -m "exp_N: try <params>=<values>"`
4. **Run** — `uv run train.py > run.log 2>&1`
5. **Evaluate** — parse the JSON objective score from `run.log`
6. **Handle crashes** — if the run fails, fix typos/imports or log "crash" and move on
7. **Record** — append a tab-separated row to `results.tsv` (do not commit this file)
8. **Keep or discard** — if objective improved, keep the commit; otherwise `git reset --hard HEAD~1`
9. **Repeat** — never stop

```bash
# Quick reference
git checkout -b autoresearch/may26
git add -A && git commit -m "exp_001: baseline defaults"
uv run train.py > run.log 2>&1
python -c "import sys,json; d=json.load(open('run.log')); print(d.get('objective','CRASH'))"
tail -n 50 run.log  # on crash
git reset --hard HEAD~1  # to revert
echo -e "exp_002\t0.535\tresolution=0.8,n_neighbors=25" >> results.tsv
```

## Evaluation Metrics

### Marker Coherence (weight: 0.4)
Evaluates whether cluster marker genes match known cell-type markers (T cell, B cell, macrophage, epithelial). For each cluster, the top 20 marker genes are compared against a curated reference; the best overlap score is averaged across clusters.

### Cluster Stability (weight: 0.3)
Measures reproducibility by subsampling 80% of the data, re-running the full pipeline, and computing the Adjusted Rand Index (ARI) between the original and subsample cluster assignments. Averaged over 3 subsample iterations.

### Silhouette Score (weight: 0.2)
Computes the silhouette coefficient on PCA embeddings using the cluster labels. Scores are normalized from [-1, 1] to [0, 1] range via `(score + 1) / 2`.

### Cluster Balance (weight: 0.1)
Measures how evenly cells are distributed across clusters using normalized entropy. A perfectly balanced clustering (all clusters equal size) scores 1.0; a single-cluster result scores 0.0.

## Results

Results are stored in `output/results.csv` with each row containing the full `AutoResearchConfig` parameters and corresponding evaluation scores. The objective improvement chart (`output/objective_improvement.png`) visualizes best-so-far progress across experiment runs.

### Best Known Configuration

Based on experiments run so far, the best objective score (~0.569) was achieved with:

```
min_genes=100, max_genes=100000, min_counts=10, max_counts=100000,
max_pct_mt=10.0, target_sum=1000, n_top_genes=4000, n_pcs=60,
svd_solver=arpack, n_neighbors=45, metric=euclidean, resolution=0.8,
random_seed=42
```

## Output Files

| File | Description |
|------|-------------|
| `output/results.csv` | CSV log of all experiment runs with parameters and scores |
| `output/umap.png` | UMAP projection colored by cluster |
| `output/objective_improvement.png` | Improvement trace of best objective over runs |
| `output/scores.json` | Latest evaluation scores (JSON) |
| `output/summary.txt` | Human-readable config and scores summary |
| `run.log` | Full stdout/stderr of the latest run |

## License

This project is for research and educational purposes.