# Autoresearch: scRNA-seq Hyperparameter Tuning

## Objective

Maximize the composite objective score of a single-cell RNA-seq clustering pipeline by iteratively tuning hyperparameters.

**Objective function** (higher is better, approximate range [-0.2, 1.0]):

```
objective = 0.4 × marker_coherence + 0.3 × cluster_stability + 0.2 × silhouette + 0.1 × cluster_balance
```

The pipeline uses **scanpy** for single-cell RNA-seq analysis with the following stages:
- QC filtering → Normalization + log1p → HVG selection → Scaling → PCA → Neighbor graph → Leiden clustering → Evaluation

---

## Setup

### 1. Create a run tag
Give each experiment a unique tag (e.g. `exp_001`, `exp_002`) to track iterations.

### 2. Verify data exists
Check that `data/lims_lung_celltype_demo.h5ad` is present. If not, run `uv run prepare.py` once to download it.

### 3. Read in-scope files for context
Read these files to understand the pipeline before making changes:

| File | Purpose |
|------|---------|
| `train.py` | The main pipeline. Contains `AutoResearchConfig` with all tunable hyperparameters. **This is the only file you may edit. You may ONLY change hyperparameter values.** |
| `prepare.py` | Data preparation script. **DO NOT MODIFY.** Read-only context. |
| `pyproject.toml` | Project dependencies. **Do not modify or install additional packages.** |

### 4. Initialize results tracking
The first time you run the pipeline (see First Run below), the output files will be created automatically:
- `output/scores.json` — JSON of the latest scores (overwritten each run)
- `output/results.csv` — Cumulative CSV log (appended each run)
- `output/summary.txt` — Human-readable summary (overwritten each run)

Read `output/results.csv` after each run to review the history.

---

## Experimentation Rules

### ✅ You CAN:
- Change **hyperparameter values** in the `AutoResearchConfig` dataclass in `train.py` (lines 20–52)
- Read any file in the repository for context
- Run `uv run train.py` to execute the pipeline
- Read `output/scores.json` to parse the objective score
- Read `output/results.csv` to review experiment history
- Keep the best configuration seen so far

### ❌ You CANNOT:
- **Modify pipeline logic, evaluation functions, or any code outside of `AutoResearchConfig` values**
- **Modify `prepare.py`** under any circumstances
- **Install additional Python packages** or modify `pyproject.toml`
- **Change the output format, file paths, or logging behaviour**
- **Fix known bugs** (see Known Bugs section)
- **Use git branches, tags, or version control** — simple CSV logging is sufficient
- **Modify `test.py` if it exists** — it is not part of this experiment

---

## Search Space

All tunable hyperparameters are in `AutoResearchConfig` (lines 20–52 of `train.py`). The agent may set any parameter to any value within the ranges below.

| # | Parameter | Type | Range | Default | Notes |
|---|-----------|------|-------|---------|-------|
| 1 | `min_genes` | `int` | [100, 500] | 200 | Minimum genes per cell (QC filter) |
| 2 | `max_genes` | `int` | [4000, 8000] | 6000 | Maximum genes per cell (QC filter) |
| 3 | `min_counts` | `int` | [200, 1000] | 500 | Minimum UMI counts per cell (QC filter) |
| 4 | `max_counts` | `int` | [30000, 80000] | 50000 | Maximum UMI counts per cell (QC filter) |
| 5 | `max_pct_mt` | `float` | [5.0, 20.0] | 10.0 | Maximum % mitochondrial counts (QC filter) |
| 6 | `target_sum` | `int` | [1e4, 1e5] | 10000 | Target sum for library-size normalization |
| 7 | `n_top_genes` | `int` | [1000, 5000] | 2000 | Number of highly variable genes to select |
| 8 | `n_pcs` | `int` | [20, 90] | 30 | Number of principal components |
| 9 | `svd_solver` | `str` | `"arpack"` or `"randomized"` | `"arpack"` | SVD solver for PCA |
| 10 | `n_neighbors` | `int` | [15, 50] | 15 | Number of neighbours for graph construction |
| 11 | `metric` | `str` | `"cosine"` or `"euclidean"` | `"cosine"` | Distance metric for neighbour graph |
| 12 | `resolution` | `float` | [0.2, 2.0] | 1.0 | Resolution parameter for Leiden clustering |
| 13 | `random_seed` | `int` | any integer | 42 | Random seed for reproducibility |

**Parameters in `AutoResearchConfig` that are NOT in the search space (do not change):**
- `flavor` — must remain `"seurat_v3"` (affects HVG selection order)
- `batch_key` — keep as `None` (single-sample dataset)
- `do_scale` — keep as `True` (enables scaling)
- `max_value` — keep as `10.0` (scaling clip value)

---

## First Run (Baseline)

Run the pipeline with **default hyperparameters** to establish a baseline:

```bash
uv run train.py
```

**Expect the first run to crash.** This is due to a known bug (see Known Bugs section). Do not attempt to fix it.

**Crash handling procedure:**
1. The script will exit with a non-zero code and a Python traceback
2. Note the hyperparameter configuration that caused the crash
3. Log the crash in `output/results.csv` by adding a row with the hyperparameter values and `objective = -999` (or a sentinel value). You can append to results.csv manually:
   ```bash
   # Example manual logging of a crash:
   echo "200,6000,500,50000,10,10000,2000,30,arpack,15,cosine,1.0,42,0,0,0,0,-999" >> output/results.csv
   ```
   (The column order matches the config fields followed by scores.)
4. Move on to the experiment loop with a new configuration

If the pipeline runs successfully despite the bug, record the scores and proceed to the experiment loop.

---

## Output Format

When `train.py` runs successfully, it prints the scores as JSON to stdout. Example:

```json
{
  "marker_coherence": 0.45,
  "cluster_stability": 0.32,
  "silhouette": 0.61,
  "cluster_balance": 0.78,
  "objective": 0.492
}
```

The objective is computed inside `train.py` by `compute_objective()` and printed as part of the score dictionary. You can extract it by parsing the JSON output.

To capture the objective score programmatically:

```bash
uv run train.py 2>&1 | python -c "import sys,json; print(json.load(sys.stdin)['objective'])"
```

---

## Logging Results

Two output files are managed automatically by `train.py`:

### `output/scores.json`
Overwritten every run. Contains only the latest scores. Read this file after each run to get the objective value.

### `output/results.csv`
Appended every run by `append_results()`. Contains one row per experiment with all hyperparameter values and all score components.

**CSV columns (in order):**
`min_genes,max_genes,min_counts,max_counts,max_pct_mt,target_sum,n_top_genes,n_pcs,svd_solver,n_neighbors,metric,resolution,random_seed,marker_coherence,cluster_stability,silhouette,cluster_balance,objective`

You should periodically read `output/results.csv` to track your best objective so far and determine which configurations to explore next.

---

## The Experiment Loop

Run autonomously in a loop until stopped. Each iteration follows these steps:

### Step 1: Tune
Select one or more hyperparameters to modify. Choose a new value within the defined search space.

**Search strategy suggestions:**
- **Grid/random search:** Sample random values within ranges to explore broadly
- **Greedy refinement:** Take your best configuration so far and nudge one parameter at a time (e.g. increase/decrease by 10–20%)
- **Bayesian-style:** Focus on parameters that most affect the objective based on history
- **Local search near good runs:** If resolution=1.2 gave the best score, try 1.1, 1.3, 1.4

Do not change ALL parameters at once — change 1–3 parameters per iteration to isolate effects.

### Step 2: Edit
Open `train.py` and change the default values in `AutoResearchConfig` (lines 20–52) to your chosen values.

**Example:** To set `resolution=1.5` and `n_neighbors=30`:
```python
    resolution: float = 1.5
    n_neighbors: int = 30
```

**You must edit the source file directly.** No command-line arguments or environment variables are supported.

### Step 3: Run
Execute the pipeline:
```bash
uv run train.py
```

### Step 4: Evaluate
- **If the script crashes:** Log the configuration with objective = -999 (or sentinel value), then go to Step 7
- **If the script succeeds:** Parse the JSON output to extract all scores, especially the `objective` value

### Step 5: Keep or Discard

| Condition | Action |
|-----------|--------|
| Objective is higher than the previous best | **Keep** — this becomes your new best configuration |
| Objective is lower than the previous best | **Discard** — revert to the previous best configuration in `AutoResearchConfig` |
| Script crashed | **Discard** — revert to the previous best configuration |
| First successful run (no prior baseline) | **Keep** — this is your baseline |

**Reversion means:** edit `train.py` to restore the `AutoResearchConfig` values back to the previous best configuration.

### Step 6: Log
The pipeline automatically appends results to `output/results.csv` on successful runs. For crashes, you need to manually append a row (see First Run section).

### Step 7: Repeat
Increment your run tag and go to Step 1.

**Summary of the loop:**
```
Tune → Edit train.py → Run pipeline → Parse objective → 
  if better: Keep (new best)
  if worse:  Discard (revert to previous best)
  if crash:  Discard (revert), log sentinel
→ Log → Repeat
```

---

## Stopping Condition

Run for **200 iterations** total (including the baseline attempt). After 200 iterations, stop and report the best configuration found along with its objective score.

If the pipeline consistently crashes for every configuration tried (e.g. more than 10 crashes in a row), you may stop early and report the issue.

---

## Known Bugs

### `subset_adata` reference in `main()` (lines 294–299)
There is a variable naming issue in the `main()` function where `subset_adata` is referenced but was never defined in that scope. This may cause the pipeline to crash during the HVG selection / normalization branching block.

**You must NOT fix this bug.** The purpose of this experiment is hyperparameter tuning, not code debugging. Acknowledge the crash, log it, revert to the previous working configuration, and move on.

If this bug causes the pipeline to crash on ALL runs, note this in your final report and terminate the experiment early.

---

## Quick Reference

```bash
# Run the pipeline
uv run train.py

# Parse objective from stdout
uv run train.py 2>&1 | python -c "import sys,json; print(json.load(sys.stdin)['objective'])"

# Read latest scores
cat output/scores.json

# Read experiment history
cat output/results.csv

# Count experiments so far
wc -l output/results.csv