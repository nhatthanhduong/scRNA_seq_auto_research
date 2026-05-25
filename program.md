# Autoresearch: scRNA-seq Hyperparameter Tuning

Maximize composite objective (higher better, range ≈[-0.2, 1.0]):

```
objective = 0.4×marker_coherence + 0.3×cluster_stability + 0.2×silhouette + 0.1×cluster_balance
```

Pipeline: QC → Normalize+log1p → HVG → Scale → PCA → Neighbors → Leiden → Eval

## Setup

1. **Tag**: `exp_001`, `exp_002`, ...
2. **Data**: verify `data/lims_lung_celltype_demo.h5ad` exists, else `uv run prepare.py`
3. **Read**: `train.py` (only file you edit, only hyperparameter values), `prepare.py` (do not touch), `pyproject.toml` (no new packages)

## Hyperparameter Search Space

Edit `AutoResearchConfig` in `train.py` (lines 20–52). Only these 13 params are tunable:

| Param | Type | Range | Default |
|-------|------|-------|---------|
| `min_genes` | int | [100, 500] | 200 |
| `max_genes` | int | [4000, 8000] | 6000 |
| `min_counts` | int | [200, 1000] | 500 |
| `max_counts` | int | [30000, 80000] | 50000 |
| `max_pct_mt` | float | [5.0, 20.0] | 10.0 |
| `target_sum` | int | [1e4, 1e5] | 10000 |
| `n_top_genes` | int | [1000, 5000] | 2000 |
| `n_pcs` | int | [20, 90] | 30 |
| `svd_solver` | str | `"arpack"` or `"randomized"` | `"arpack"` |
| `n_neighbors` | int | [15, 50] | 15 |
| `metric` | str | `"cosine"` or `"euclidean"` | `"cosine"` |
| `resolution` | float | [0.2, 2.0] | 1.0 |
| `random_seed` | int | any | 42 |

**Do not change**: `flavor` (`"seurat"`), `batch_key` (`None`), `do_scale` (`True`), `max_value` (`10.0`)

## First Run (Baseline)

```
uv run train.py
```

## Experiment Loop

Repeat until 200 runs. Each iteration:

1. **Tune** — pick 1–3 params, set new values within ranges above
2. **Edit** — change defaults in `AutoResearchConfig` in `train.py`
3. **Run** — `uv run train.py`
4. **Read score** — parse JSON stdout or `output/scores.json`
5. **Keep/Discard/Commit**:
   - Score > best → **Keep** (new best), **git commit** with message `exp_<tag>: new best objective=<score> with <params>=<vals>, ...`
   - Score ≤ best → **Discard** (revert `train.py` to previous best config)
   - Crash → **Discard** (revert), log sentinel row
6. **Log** — auto-logged on success; manual CSV append on crash
7. **Repeat** — increment tag

```
Tune → Edit → Run → Evaluate →
  if better: Keep + git commit
  if worse:  Discard (revert)
  if crash:  Discard, log -999
→ Log → Repeat
```

## Stopping

200 iterations total. If >10 crashes in a row, stop early.

## Quick Reference

```bash
uv run train.py                                            # run pipeline
uv run train.py 2>&1 | python -c "import sys,json; print(json.load(sys.stdin)['objective'])"  # parse score
cat output/scores.json                                     # latest scores
cat output/results.csv                                     # experiment history
git add -A && git commit -m "exp_N: new best objective=..." # commit new best