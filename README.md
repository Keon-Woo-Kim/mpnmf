# mpnmf (Metaprogram NMF)

Metaprogram discovery via non-negative matrix factorization for single-cell RNA-seq data.


## Overview

`mpnmf` is a Python implementation of the metaprogram discovery method described in Gavish et al. (2023, *Nature*). It identifies recurrent transcriptional programs — "metaprograms" — across samples in single-cell RNA-seq data through three steps:

1. **Per-sample NMF** (`run`) — factorizes each sample's expression matrix across a range of ranks.
2. **Program refinement** (`refine`) — retains programs that are intra-sample reproducible and inter-sample recurrent while removing intra-sample redundancy.
3. **Program clustering** (`cluster`) — iteratively merges filtered programs sharing gene overlap into metaprograms.


## Differences from the original method

### Deterministic NMF initialization

The original R implementation runs NMF with random initialization and averages results across multiple runs. For each sample, `mpnmf` fits NMF once per rank using NNDSVDa initialization (Boutsidis & Gallopoulos, 2008), which is deterministic and produces identical output on repeated runs. This yields reproducible outputs without the need for consensus averaging, reducing computational cost substantially.

### Preprocessing modes

Before NMF, `mpnmf` applies gene selection, centering, and optional scaling. Two modes are offered, with defaults matching common analytical conventions:

- **HEG mode** (`mode='heg'`, `scale=False` by default): top genes by mean expression, centered per gene, clipped to non-negative. Matches the original method.
- **HVG mode** (`mode='hvg'`, `scale=True` by default): highly variable genes (dispersion-based) filtered further by minimum mean expression, centered and divided by gene standard deviation, clipped to non-negative. Full z-normalization amplifies lowly expressed but strongly variable genes, enabling detection of rare or trace signals.

The default scaling behavior in each mode can be overridden via the scale argument.


## Installation

**Requirements:** Python ≥ 3.9.

We recommend installing `mpnmf` in a dedicated conda environment to avoid dependency conflicts with other single-cell tools:

```bash
conda create -n mpnmf
conda activate mpnmf
pip install mpnmf
```

Or install into an existing environment:

```bash
pip install mpnmf
```


## Usage

```python
import scanpy as sc
import mpnmf

adata = sc.read_h5ad("your_data.h5ad")        # anndata should be log-normalized

sample_key  = "batch"
sample_list = adata.obs[sample_key].unique().tolist()
krange      = range(7, 13)

# NMF run: HVG mode
nmf_run     = mpnmf.run(adata, krange=krange, sample_key=sample_key, sample_list=sample_list, mode="hvg", n_top_genes=7000, scale=True, title="test")
# NMF run: HEG mode
nmf_run     = mpnmf.run(adata, krange=krange, sample_key=sample_key, sample_list=sample_list, mode="heg", n_top_genes=7000, scale=False, title="test")

# Program refinement: intra-sample reproducibility, inter-sample recurrence, intra-sample non-redundancy
nmf_refined = mpnmf.refine(nmf_run, thres_intra=0.7, thres_inter=0.2, thres_redun=0.2, title="test")

# Metaprogram clustering: iteratively merge programs into metaprograms
nmf_df      = mpnmf.cluster(nmf_refined, thres_overlap=0.3, min_overlap=5, title="test")
```


## Input requirements

- `adata.X` is log-normalized expression (not raw counts, not z-scored).
- `adata.var_names` contains unique gene symbols.
- `adata.obs` contains a column identifying the sample of each cell.
- Each sample has enough cells to factorize at `max(krange)` (rule of thumb: ≥ 50).


## APIs

### `mpnmf.run(adata, krange, sample_key, sample_list, ...)`

Runs NMF per sample across a range of ranks.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `adata` | — | Log-normalized AnnData object. |
| `krange` | — | Iterable of NMF ranks to try (e.g., `range(4, 10)`). |
| `sample_key` | — | Column in `adata.obs` used to split cells by sample. |
| `sample_list` | — | List of sample values to run NMF on. |
| `n_genes` | `50` | Number of top genes retained per program. |
| `max_iter` | `5000` | Max NMF iterations per fit. |
| `mode` | `'hvg'` | Gene selection: `'hvg'` (dispersion-based) or `'heg'` (mean expression). |
| `n_top_genes` | `7000` | Number of genes kept after selection. |
| `min_exp_pct` | `0.2` | In HVG mode, drop bottom fraction by mean expression. |
| `scale` | `'auto'` | Whether to divide by gene std after centering. `'auto'` = `True` for HVG, `False` for HEG. Centering is always applied regardless. |
| `title` | `None` | Prefix for output files; defaults to `"mpnmf"`. |
| `savepath` | `None` | Output directory; defaults to `./mpnmf/`. |

**Returns:** `nmf_run` dict, keyed by sample → rank → `{W, H, rank}`.

### `mpnmf.refine(nmf_run, ...)`

Filters programs through three sequential criteria: intra-sample reproducibility, inter-sample recurrence, and intra-sample non-redundancy.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nmf_run` | — | Output of `mpnmf.run`. |
| `samples` | `None` | Subset of samples to use; defaults to all. |
| `krange` | `None` | Subset of ranks to use; defaults to all. |
| `n_genes` | `None` | Program length; inferred from `nmf_run` if not given. |
| `thres_intra` | `0.7` | Min fraction of top genes shared with another rank in the same sample. |
| `thres_inter` | `0.2` | Min fraction of top genes shared with the best-matching program in another sample. |
| `thres_redun` | `0.2` | Max allowed overlap with programs already kept in the same sample. |
| `title` | `None` | Prefix for output files; defaults to `"mpnmf"`. |
| `savepath` | `None` | Output directory; defaults to `./mpnmf/`. |

**Returns:** `nmf_refined` dict, keyed by program name → `{genes, scores}`.

### `mpnmf.cluster(nmf_refined, ...)`

Iteratively merges refined programs into metaprograms by gene overlap.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nmf_refined` | — | Output of `mpnmf.refine`. |
| `n_genes` | `50` | Expected program length; all programs must match. |
| `thres_overlap` | `0.3` | Min fraction of shared genes to merge two programs. |
| `min_overlap` | `5` | Min number of qualifying partners for a program to seed a new metaprogram. |
| `title` | `None` | Prefix for output files; defaults to `"mpnmf"`. |
| `savepath` | `None` | Output directory; defaults to `./mpnmf/`. |

**Returns:** `nmf_df`, a DataFrame of genes × metaprograms.


## Output files

| Function | File | Content |
|----------|------|---------|
| `run` | `{prefix}_run.pkl` | `nmf_run` dict |
| `refine` | `{prefix}_refined.pkl` | `nmf_refined` dict |
| `cluster` | `{prefix}_clustered.pkl` | `MP_dict` (genes + scores + freq per MP) |
| `cluster` | `{prefix}.csv` | Gene × MP table |

`{prefix}` = `title` if given, else `"mpnmf"`.


## Citation

If you use `mpnmf` in your research, please cite the original paper:

> Gavish, A., Tyler, M., Greenwald, A.C., et al. Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours. *Nature* 618, 598–606 (2023).