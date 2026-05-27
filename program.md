# Autoresearch: scRNA-seq Hyperparameter Tuning

Maximize composite objective (higher better, range ≈[-0.2, 1.0]):

```
objective = 0.4×marker_coherence + 0.3×cluster_stability + 0.2×silhouette + 0.1×cluster_balance
```

Pipeline: QC → Normalize+log1p → HVG → Scale → PCA → Neighbors → Leiden → Eval

## Setup

- Run on branch `autoresearch/may26`
- Verify `data/lims_lung_celltype_demo.h5ad` exists; if not `uv run prepare.py`
- Read `train.py` — only file you edit. Only change `AutoResearchConfig` values.
- DO NOT modify `prepare.py`, `pyproject.toml`, or pipeline logic.

## Search Space

Edit `AutoResearchConfig` (train.py lines 20–52). Tunable params:

| Param | Type | Range | Default |
|-------|------|-------|---------|
| `min_genes` | int | any | 1 |
| `max_genes` | int | any | 10 |
| `min_counts` | int | any | 100 |
| `max_counts` | int | any | 1000 |
| `max_pct_mt` | float | [1, 100] | 1 |
| `target_sum` | int | any | 1000 |
| `n_top_genes` | int | any | 100 |
| `n_pcs` | int | any | 10 |
| `svd_solver` | str | `"arpack"` / `"randomized"` | `"arpack"` |
| `n_neighbors` | int | any | 20 |
| `metric` | str | `"cosine"` / `"euclidean"` | `"cosine"` |
| `resolution` | float | [0.1, 10.0] | 0.1 |
| `random_seed` | int | any | 42 |

**Do not touch**: `flavor`, `batch_key`, `do_scale`, `max_value`

## First Run

```
uv run train.py 2>&1 | tee run.log
```

## Experiment Loop

LOOP FOREVER:

1. **Look at git state** — note the current branch/commit
2. **Tune** — pick 1–3 params, directly hack the defaults in `AutoResearchConfig`
3. **Commit** — `git add -A && git commit -m "exp_N: try <params>=<vals>"` (commit everything)
4. **Run** — `uv run train.py > run.log 2>&1` (redirect everything, NO tee, NO output flood)
5. **Read results** — grep the scores from the JSON in `run.log`. If `objective` found, extract it
6. **If empty (crash)** — `tail -n 50 run.log` to see stack trace. If it's a dumb bug (typo, missing import) fix it and re-run. If fundamentally broken, log "crash" in the tsv and move on
7. **Record** — append row to `results.tsv` (tab-separated). DO NOT commit `results.tsv`, leave it untracked
8. **Keep/Discard**:
   - If `objective` improved (higher) → advance, keep the commit
   - If `objective` equal or worse → `git reset --hard HEAD~1` (discard the commit, revert code)
   - If crash → `git reset --hard HEAD~1` (discard the commit)
9. **NEVER STOP** (go to step 1)

```
Look at git state → Tune → Commit → Run → Read score →
  if improved: Keep commit
  if worse/crash: git reset --hard HEAD~1
→ Record in tsv → Repeat (forever)
```

## Rules

- **NEVOR STOP**: Once the loop begins, do NOT pause to ask the human. Do NOT ask "should I continue?". The human might be asleep. You are autonomous. If out of ideas, think harder. Loop until manually interrupted.
- **Timeout**: Each run should take ~5 min. If >10 min, kill it, treat as failure, revert.
- **TSV, not CSV**: Append tab-separated rows to `results.tsv`. Do not commit this file.
- **Commit every try**: git commit BEFORE running. This way you can always `git reset --hard HEAD~1` to revert.
- **Crashes**: Dumb bug (typo, missing import)? Fix it, re-run. Fundamental broken idea? Log "crash", revert, move on.
- **No scripts**: Do NOT write automation scripts. Hack `train.py` directly.


## Quick Reference

```bash
# Create branch
git checkout -b autoresearch/may26

# Initial commit of working state
git add -A && git commit -m "exp_001: baseline defaults"

# Run
uv run train.py > run.log 2>&1

# Parse score
python -c "import sys,json; d=json.load(open('run.log')); print(d.get('objective','CRASH'))"

# See crash
tail -n 50 run.log

# If worse/crash, revert
git reset --hard HEAD~1

# Record results (tab-separated, do NOT commit)
echo -e "exp_002\t0.535\tresolution=0.8,n_neighbors=25" >> results.tsv