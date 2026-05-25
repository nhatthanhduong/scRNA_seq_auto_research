# Implementation Plan

## Overview

Create a `program.md` document that describes an autonomous hyperparameter search protocol (following Andrej Karpathy's autoresearch format) for tuning single-cell RNA-seq clustering pipeline hyperparameters to maximize a composite objective score.

The document instructs an AI coding agent to iteratively modify the `AutoResearchConfig` dataclass in `train.py` within defined search spaces, run the pipeline, evaluate results via the composite objective score (0.4×marker_coherence + 0.3×cluster_stability + 0.2×silhouette + 0.1×cluster_balance), and log outcomes to `output/results.csv`. The agent operates autonomously in a loop — trying hyperparameters, keeping improvements, and discarding regressions — until manually stopped.

**Important context**: The agent is ONLY permitted to change hyperparameter values in the `AutoResearchConfig` class, never the pipeline logic or evaluation code.

## Types

Define the hyperparameter search space as typed constraints with min/max bounds and valid categorical choices.

### Hyperparameter Search Space

| Parameter | Type | Range / Choices | Default |
|-----------|------|-----------------|---------|
| `min_genes` | int | [100, 500] | 200 |
| `max_genes` | int | [4000, 8000] | 6000 |
| `min_counts` | int | [200, 1000] | 500 |
| `max_counts` | int | [30000, 80000] | 50000 |
| `max_pct_mt` | float | [5.0, 20.0] | 10 |
| `target_sum` | float | [1e4, 1e5] | 10000 |
| `n_top_genes` | int | [1000, 5000] | 2000 |
| `n_pcs` | int | [20, 90] | 30 |
| `svd_solver` | str | {"arpack", "randomized"} | "arpack" |
| `n_neighbors` | int | [15, 50] | 15 |
| `metric` | str | {"cosine", "euclidean"} | "cosine" |
| `resolution` | float | [0.2, 2.0] | 1.0 |
| `random_seed` | int | any int | 42 |

### Objective Score Definition

```python
objective = (
    0.4 * marker_coherence +
    0.3 * cluster_stability +
    0.2 * silhouette +
    0.1 * cluster_balance
)
```

Higher is better. Score range is approximately [-0.2, 1.0] based on the components.

### Results TSV Schema (output/results.csv)

Current format in train.py is CSV, not TSV. The `program.md` should instruct the agent to read results from the existing `output/scores.json` file (written by `save_scores()`) and log to `output/results.csv` (already handled by `append_results()`).

## Files

Create one new file: `program.md` at the repository root.

**No existing files are modified.** The agent is only allowed to change hyperparameter *values* in `train.py`'s `AutoResearchConfig`, not its logic or structure. The `program.md` instructs the agent on how to do this.

### program.md Structure (following Karpathy autoresearch format)

The document will have these sections:

1. **Title** - Autoresearch: scRNA-seq Hyperparameter Tuning
2. **Setup** - Steps to: create run tag, verify data exists, read in-scope files, initialize results
3. **Experimentation rules** - What the agent CAN do (change hyperparameter values in AutoResearchConfig) and CANNOT do (modify prepare.py, change pipeline logic, install packages)
4. **Search space** - Table of all tunable parameters with ranges
5. **First run** - Establish baseline (which will crash due to the subset_adata bug — note it and move on)
6. **Output format** - What the script prints (scores from main())
7. **Logging results** - How to record to results.csv
8. **The experiment loop** - Autonomous loop: tune → run → evaluate → keep/discard → repeat
9. **Stopping** - Run for 200 iterations
10. **Known bugs** - Document the subset_adata bug

## Functions

No new functions. The `program.md` is a documentation/instruction file, not executable code.

The document may reference these existing functions from train.py (read-only context for the agent):
- `main()` — entry point, runs pipeline and prints scores
- `compute_objective(scores)` — computes the composite objective from sub-scores
- `save_scores(scores)` — writes `output/scores.json`
- `append_results(config, scores)` — appends to `output/results.csv`
- `cluster_stability_score(...)` — subsampling-based stability (takes ~3× longer than main pipeline)
- `marker_coherence_score(...)` — marker gene overlap
- `compute_silhouette_score(...)` — silhouette on PCA space
- `cluster_balance_score(...)` — entropy of cluster sizes

## Classes

No new classes. Referenced class (read-only):
- `AutoResearchConfig` in `train.py` (lines 20-52) — a `@dataclass` with all tunable hyperparameters as fields. The agent changes values in this class directly in the source code.

## Dependencies

No new dependencies. The existing `pyproject.toml` includes all required packages: anndata, cellxgene-census, matplotlib, numpy, pandas, scanpy[leiden], scikit-learn, seaborn.

## Testing

No unit tests. Validation is performed by running `uv run train.py` and checking:
1. Exit code (0 = success, non-zero = crash)
2. Grepping `objective` from the JSON output to extract the score
3. Recording the result in results.csv

If the script crashes, the agent should note the hyperparameter choice, log the crash, and move on.

## Implementation Order

1. **Write `program.md`** — Create the full instruction document for the AI agent following the structure described above
2. **Verify structure** — Ensure all sections from Karpathy's format are present and adapted correctly
3. **Cross-check search space** — Verify each hyperparameter in the search space matches the fields in `AutoResearchConfig`
4. **Final review** — Ensure all rules (CAN/CANNOT, crash handling, keep/discard logic) are clearly stated