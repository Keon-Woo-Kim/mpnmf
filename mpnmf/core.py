import os
import pickle
import time
import warnings
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import ImplicitModificationWarning
from scipy.sparse import issparse
from sklearn.decomposition import NMF


# =============================================================================
# Internal helpers
# =============================================================================

def _resolve_savepath(savepath, title):
    # Default: {cwd}/mpnmf/ ; else user-specified directory.
    if savepath is None:
        savepath = os.path.join(os.getcwd(), 'mpnmf')
    os.makedirs(savepath, exist_ok=True)
    # Prefix for output filenames (defaults to 'mpnmf' if no title).
    prefix = title if title is not None else 'mpnmf'
    return savepath, prefix


def _run_single(adata, krange, sample_key, sample, n_genes, max_iter):
    # Runs NMF on a single sample across the full range of ranks (k).

    gene_names = adata.var_names
    # Subset expression matrix to current sample.
    temp = adata[adata.obs[sample_key] == sample].X
    if temp.shape[0] == 0:
        raise ValueError(f"Sample '{sample}' has no cells.")
    if temp.shape[0] < max(krange):
        raise ValueError(f"Sample '{sample}' has {temp.shape[0]} cells, fewer than max(krange)={max(krange)}.")
    result = {}
    start = time.time()

    # Iterate over ranks; fit NMF independently at each k.
    for k in krange:
        # Deterministic NMF (nndsvda init) — reproducible across runs.
        model = NMF(n_components=k, init='nndsvda', max_iter=max_iter, solver='cd', beta_loss='frobenius', random_state=0)
        W = model.fit_transform(temp)   # W: (n_cells, k) cell loadings
        H = model.components_            # H: (k, n_genes) gene loadings

        # For each program (row of H), rank genes by loading and keep top n_genes.
        rank = {}
        for i, program in enumerate(H):
            top_idx = np.argsort(program)[-n_genes:][::-1]
            rank[f'{sample}-k{k}-program{i+1}'] = gene_names[top_idx].tolist()

        # Store W, H, and ranked gene lists for this rank.
        result[k] = {'W': W, 'H': H, 'rank': rank}

    elapsed = (time.time() - start) / 60
    end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    print(f"{sample} (n={temp.shape[0]}) finished! (elapsed: {elapsed:.2f} min / {end_time})")

    return result


def _overlap_df(d):
    # Build a pairwise overlap matrix: number of shared top genes between programs.
    df = pd.DataFrame(index=d.keys(), columns=d.keys(), dtype=int).fillna(0)
    for p1, v1 in d.items():
        for p2, v2 in d.items():
            df.at[p1, p2] = len(set(v1["genes"]).intersection(v2["genes"]))
    return df


def _merge_programs(p1, p2, gene_history, n_genes=50):
    # Merge two programs into one metaprogram representative.
    # Priority: (1) gene recurrence across merged partners, (2) max score.

    combined = {}
    freq = Counter(gene_history)
    # Collect best score per gene across both programs.
    for g, s in zip(p1["genes"], p1["scores"]):
        combined[g] = max(combined.get(g, 0), s)
    for g, s in zip(p2["genes"], p2["scores"]):
        combined[g] = max(combined.get(g, 0), s)

    # Pool all genes seen so far (from full history, not just the current pair).
    all_genes = set(gene_history)
    # Rank by frequency first, then by score as tiebreaker.
    sorted_genes = sorted(all_genes, key=lambda g: (freq[g], combined.get(g, 0)), reverse=True)
    top = sorted_genes[:n_genes]
    return {"genes": top, "scores": [combined.get(g, 0) for g in top], "freq": [freq[g] for g in top]}


# =============================================================================
# Main API
# =============================================================================

def run(adata, krange, sample_key, sample_list, n_genes=50, max_iter=5000,
        mode='hvg', n_top_genes=7000, min_exp_pct=0.2, scale='auto',
        title=None, savepath=None):

    start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    t0 = time.time()
    print(f"\n=== NMF run started. (starting time: {start_time}) ===\n")

    adata_nmf = adata.copy()

    # Gene selection + preprocessing
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ImplicitModificationWarning)

        # Gene selection: HVG (dispersion-based) or HEG (mean expression)
        if mode == 'hvg':
            # Keep top n_top_genes highly variable genes (scanpy's dispersion method).
            sc.pp.highly_variable_genes(adata_nmf, n_top_genes=n_top_genes, subset=True)
            # Further drop genes below min_exp_pct percentile of mean expression.
            # This removes HVGs that are highly variable but barely expressed.
            if min_exp_pct > 0.0:
                mean_expr = np.asarray(adata_nmf.X.mean(axis=0)).flatten()
                thres = np.percentile(mean_expr, min_exp_pct * 100)
                adata_nmf = adata_nmf[:, mean_expr >= thres]
        elif mode == 'heg':
            # Keep top n_top_genes by mean expression.
            mean_expr = np.asarray(adata_nmf.X.mean(axis=0)).flatten()
            top_idx = np.sort(np.argsort(mean_expr)[-n_top_genes:])
            adata_nmf = adata_nmf[:, top_idx]
        else:
            raise ValueError(f"mode must be 'hvg' or 'heg', got '{mode}'")

        # Materialize subset to avoid view-assignment issues downstream.
        adata_nmf = adata_nmf.copy() if adata_nmf.is_view else adata_nmf

        # Centering of gene expression + optionally scale by std
        # 'auto': scale in HVG mode (amplifies rare but strongly variable genes), skip in HEG mode.
        X = adata_nmf.X.toarray() if issparse(adata_nmf.X) else adata_nmf.X.copy()
        X = X - X.mean(axis=0)
        if scale == 'auto':
            scale_flag = (mode == 'hvg')
        else:
            scale_flag = bool(scale)
        if scale_flag:
            stds = X.std(axis=0, ddof=0)
            stds[stds == 0] = 1
            X = X / stds

        # Clip negative values (NMF requires non-negative input)
        adata_nmf.X = np.maximum(X, 0)

    # Run NMF per sample
    nmf_run = {}
    for sample in sample_list:
        nmf_run[sample] = _run_single(adata_nmf, krange, sample_key, sample, n_genes, max_iter)

    # Save + report
    savepath, prefix = _resolve_savepath(savepath, title)
    save_file = os.path.join(savepath, f"{prefix}_run.pkl")
    with open(save_file, 'wb') as f:
        pickle.dump(nmf_run, f)
    print(f"\n** Saved: {save_file} **")

    end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    elapsed = (time.time() - t0) / 60
    print(f"\n=== NMF run finished! (elapsed: {elapsed:.2f} min / {end_time}) ===\n")

    return nmf_run


def refine(nmf_run, samples=None, krange=None, n_genes=None, thres_intra=0.7, thres_inter=0.2, thres_redun=0.2, title=None, savepath=None):

    start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    t0 = time.time()
    print(f"\n=== NMF refine started. (starting time: {start_time}) ===\n")

    # Defaults pulled from the input structure.
    samples = samples or list(nmf_run.keys())
    krange = krange or list(nmf_run[samples[0]].keys())
    n_genes = n_genes or len(list(nmf_run[samples[0]][krange[0]]['rank'].values())[0])

    # Stage 1: Intra-sample reproducibility
    # Keep programs whose top genes reappear across ranks in the same sample.
    filt1 = defaultdict(dict)
    for s in samples:
        for k1 in krange:
            for p, v1 in nmf_run[s][k1]['rank'].items():
                # Reproducible if any program at a DIFFERENT rank shares >= thres_intra of top genes.
                if any(len(set(v1).intersection(v2)) >= thres_intra * n_genes
                       for k2 in krange if k2 != k1
                       for v2 in nmf_run[s][k2]['rank'].values()):
                    scores = nmf_run[s][k1]['H'][int(p.split('-program')[-1]) - 1]
                    filt1[s][p] = {"genes": v1, "scores": scores.tolist()}

    # Stage 2: Inter-sample recurrence
    # Keep programs whose top genes also appear (partially) in other samples.
    filt2, filt3 = defaultdict(dict), defaultdict(dict)
    for s1 in samples:
        for p1, v1 in filt1[s1].items():
            # Max overlap with any program in every other sample.
            # Skip empty samples and handle single-sample edge case safely.
            overlap = []
            for s2 in samples:
                if s2 == s1 or not filt1[s2]:
                    continue
                overlap.append(max(len(set(v1["genes"]).intersection(v2["genes"])) for v2 in filt1[s2].values()))
            if not overlap:
                continue
            max_overlap = max(overlap)
            # Retain if recurrence passes inter-sample threshold.
            if max_overlap >= thres_inter * n_genes:
                filt2[s1][p1] = {"max_overlap": max_overlap, "genes": v1["genes"], "scores": v1["scores"]}
        # Sort programs within sample by recurrence strength (high first).
        sorted_keys = sorted(filt2[s1], key=lambda p: filt2[s1][p]["max_overlap"], reverse=True)
        filt3[s1] = {p: {"genes": filt2[s1][p]["genes"], "scores": filt2[s1][p]["scores"]} for p in sorted_keys}

    # Stage 3: Intra-sample non-redundancy
    # Within each sample, drop programs that are too similar to already-kept ones.
    # Ordering by recurrence (stage 2) ensures the most recurrent program is kept.
    filt4 = defaultdict(dict)
    for s in samples:
        for p1, v in filt3[s].items():
            if all(len(set(v["genes"]).intersection(v2["genes"])) <= thres_redun * n_genes
                   for v2 in filt4[s].values()):
                filt4[s][p1] = {"genes": v["genes"], "scores": v["scores"]}

    # Flatten per-sample dicts into a single dict of surviving programs.
    nmf_refined = {p: v for s in samples for p, v in filt4[s].items()}
    print(f"\n** After refining, {len(nmf_refined)} programs are left. **\n")

    # Save + report.
    savepath, prefix = _resolve_savepath(savepath, title)
    save_file = os.path.join(savepath, f"{prefix}_refined.pkl")
    with open(save_file, 'wb') as f:
        pickle.dump(nmf_refined, f)
    print(f"** Saved: {save_file} **")
    
    end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    elapsed = (time.time() - t0) / 60
    print(f"\n=== NMF refine finished! (elapsed: {elapsed:.2f} min / {end_time}) ===\n")

    return nmf_refined


def cluster(nmf_refined, n_genes=50, thres_overlap=0.3, min_overlap=5, title=None, savepath=None):

    start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    t0 = time.time()
    print(f"\n=== NMF cluster started. (starting time: {start_time}) ===\n")
    
    if not nmf_refined:
        raise ValueError("nmf_refined is empty; nothing to cluster.")

    # Sanity check: all programs must have the same length of top genes.
    if set([len(p["genes"]) for p in nmf_refined.values()]) != {n_genes}:
        raise ValueError(f"Each program must contain equal {n_genes} number of genes.")

    n_total = len(nmf_refined)
    MP_dict = {}
    idx = 1

    # Pairwise overlap matrix; reused and shrunk as programs are consumed.
    overlap_df = _overlap_df(nmf_refined)

    # Iterative metaprogram assembly
    # Each iteration: pick founder -> greedily merge neighbors -> record MP.
    while nmf_refined:
        # Re-rank remaining programs by how many partners they have above threshold.
        # Subtract 1 to exclude the diagonal (self-overlap = n_genes is always above threshold).
        overlap_count = overlap_df.apply(lambda row: (row >= thres_overlap * n_genes).sum() - 1, axis=1).to_dict()
        nmf_refined = dict(sorted(nmf_refined.items(), key=lambda item: overlap_count[item[0]], reverse=True))

        # Find founder program (highest partner count, above min_overlap) 
        founder = None
        for p, data in nmf_refined.items():
            if overlap_count[p] >= min_overlap:
                founder = p
                founder_data = data
                break
        # No more founder candidates => stop.
        if founder is None:
            break

        # Initialize the metaprogram with the founder.
        MP = founder_data
        added = {founder}
        del nmf_refined[founder]
        # Full gene history across all merged partners (used for frequency tally).
        gene_history = list(founder_data["genes"])

        # Merge the single most-overlapping partner until none qualifies 
        while True:
            if not nmf_refined:
                break
            # Compute overlap of every remaining program with the current MP.
            mp_genes = set(MP["genes"])
            overlaps = {p: len(mp_genes & set(v["genes"])) for p, v in nmf_refined.items()}
            max_p = max(overlaps, key=overlaps.get)
            max_v = overlaps[max_p]
            if max_v >= thres_overlap * n_genes:
                # Merge: extend history, update MP representative.
                gene_history.extend(nmf_refined[max_p]["genes"])
                MP = _merge_programs(MP, nmf_refined[max_p], gene_history=gene_history, n_genes=n_genes)
                added.add(max_p)
                del nmf_refined[max_p]
            else:
                break

        # Finalize this MP and remove its members from the pool.
        MP_name = f"MP{idx}"
        MP_dict[MP_name] = MP
        idx += 1
        print(f"{MP_name} completed with {len(added)} programs! (Remaining: {len(nmf_refined)}/{n_total})")
        overlap_df = overlap_df.drop(index=list(added), columns=list(added), errors='ignore')

    # Assemble the human-readable gene × MP table.
    nmf_df = pd.DataFrame({name: pd.Series(mp["genes"]) for name, mp in MP_dict.items()})

    # Save + report.
    savepath, prefix = _resolve_savepath(savepath, title)
    pkl_file = os.path.join(savepath, f"{prefix}_clustered.pkl")
    csv_file = os.path.join(savepath, f"{prefix}.csv")
    with open(pkl_file, 'wb') as f:
        pickle.dump(MP_dict, f)
    nmf_df.to_csv(csv_file, index=False)
    print(f"\n** Saved: {pkl_file} **")
    print(f"** Saved: {csv_file} **")

    end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    elapsed = (time.time() - t0) / 60
    print(f"\n=== NMF cluster finished! (elapsed: {elapsed:.2f} min / {end_time}) ===\n")

    return nmf_df
