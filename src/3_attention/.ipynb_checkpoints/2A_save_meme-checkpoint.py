
# jobs 2 3 4 5 6 pending

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
conda activate dnabert_aug_2025_jupyter

DNABERT2-compatible visualization and motif region extraction
-------------------------------------------------------------

This script refactors a notebook-style workflow into a clean CLI.

Per selected job, it:
  1) Loads attention parquet (2_Attention_motif.parquet) and prediction CSV (preds/2_combined.csv)
  2) Keeps only correctly predicted sequences (true_label == pred_label)
  3) Merges by (Sequence, Label); computes per-(Sequence, layer) specificity = (max - mean) / std
  4) Keeps top-3 layers by average specificity + the last layer present in the data
  5) From mask==1 rows, counts motif hits (Aho-Corasick) in pos/neg sequence sets
  6) Performs hypergeometric enrichment + BH-FDR; keeps adj_pval < alpha
  7) Merges similar motifs using a simple gap-free alignment score
  8) Extracts 2*flank bp windows centered on motif spans; filters to Label==1 and motifs with >= min_instances
  9) Saves: 3_sig_motifs.csv, 4sig_motifs_prePWM.csv, FASTA files (weblogo_label_{label}.fa), and a histogram figure

Usage examples
--------------
# Common case: read jobs from your JSON file and run job_index
python /data/private/psurana/TSProm/src/3_attention/2A_save_meme.py \
  --jobs_file /data/private/psurana/TSProm/src/files/jobs1.json \
  --job_index 0

# Fallback: collect jobs from directories (if you don't have a jobs file)
python /data/private/psurana/TSProm/src/3_attention/2A_save_meme.py \
  --base_path /path/to/runs_root \
  --data_path  /path/to/data_root \
  --length_nm 3k \
  --base_model_nms TSp_vs_nonProm,TSp_vs_genNullseqs \
  --subsets tspAll_,tspliver_,tsptestis_,tspbrain_ \
  --lr_prefix lr3e-5 \
  --refresh_jobs \
  --job_index 0

Notes
-----
- If --jobs_file exists, it will be used. Use --refresh_jobs to override and rebuild from directories.
- Motif window size and FDR alpha are configurable.
"""

import os
import re
import json
import argparse
from itertools import combinations
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

# plotting (saved to file; no interactive display)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

try:
    import ahocorasick  # python-aho-corasick
except ImportError as e:
    raise SystemExit("python-aho-corasick is required. Please install via: pip install pyahocorasick")

# ─────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────
def log(msg: str):
    print(f"[dnabert2-motif] {msg}", flush=True)


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def _try_load_jobs_any_format(jobs_file: str) -> List[Dict]:
    """
    Load jobs from a JSON file. Supports:
      - JSON array of objects (recommended)
      - JSON Lines (one object per line)
    """
    with open(jobs_file, "r") as f:
        content = f.read().strip()

    # Try as a proper JSON array
    try:
        data = json.loads(content)
        if isinstance(data, dict):  # single object -> wrap
            return [data]
        if isinstance(data, list):
            return data
    except json.JSONDecodeError:
        pass

    # Try JSON Lines
    jobs = []
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            obj = json.loads(line)
            if isinstance(obj, dict):
                jobs.append(obj)
        except json.JSONDecodeError:
            continue
    if jobs:
        return jobs

    raise ValueError(f"Could not parse jobs from {jobs_file}. "
                     "Expected a JSON array of objects or JSON Lines (one object per line).")


def collect_jobs(base_path: str,
                 base_model_nms: List[str],
                 length_nm: str,
                 subsets: List[str],
                 data_path: str,
                 lr_prefix: str = "lr3e-5") -> List[Dict]:
    """
    Walk the training directory structure and collect jobs.
    Expected structure:
      {base_path}/{base_model_nm}/{length_nm}/{subset}/{lr_dir}/checkpoint-XXXX
    Outputs to:
      {base_path}/RESULT/{lr_dir}/{base_model_nm}_{length_nm}_{subset}
    """
    jobs: List[Dict] = []
    for base_model_nm in base_model_nms:
        length_dir = os.path.join(base_path, base_model_nm, length_nm)
        if not os.path.isdir(length_dir):
            log(f"⚠️ Missing length dir: {length_dir}")
            continue

        all_subdirs = [d for d in os.listdir(length_dir)
                       if os.path.isdir(os.path.join(length_dir, d))]
        matched_subdirs = [d for d in all_subdirs if any(d.startswith(pref) for pref in subsets)]
        if not matched_subdirs:
            log(f"⚠️ No matching subset dirs under {length_dir}")
            continue

        for subset in matched_subdirs:
            search_dir = os.path.join(length_dir, subset)
            lr_dirs = [d for d in os.listdir(search_dir) if d.startswith(lr_prefix)]
            if not lr_dirs:
                log(f"⚠️ No {lr_prefix} dirs under {search_dir}")
                continue

            for lr_dir in lr_dirs:
                model_dir = os.path.join(search_dir, lr_dir)
                checkpoints = [d for d in os.listdir(model_dir) if d.startswith("checkpoint-")]
                if not checkpoints:
                    log(f"⚠️ No checkpoints in {model_dir}")
                    continue

                # Pick latest checkpoint by step number if possible
                try:
                    latest_ckpt = sorted(checkpoints, key=lambda x: int(x.split("-")[1]))[-1]
                except Exception:
                    latest_ckpt = sorted(checkpoints)[-1]  # fallback lexicographic if parsing fails

                model_path = os.path.join(model_dir, latest_ckpt)
                data_dir = os.path.join(data_path, base_model_nm, length_nm, subset)
                res_pdir = os.path.join(base_path, "RESULT", lr_dir, f"{base_model_nm}_{length_nm}_{subset}")
                ensure_dir(res_pdir)

                jobs.append({
                    "base_model_nm": base_model_nm,
                    "subset": subset,
                    "lr_dir": lr_dir,
                    "model_path": model_path,
                    "data_dir": data_dir,
                    "res_pdir": res_pdir,
                })
                log(f"✓ Job collected: {base_model_nm} | {subset} | {lr_dir} | {os.path.basename(model_path)}")
    return jobs


def save_jobs(jobs: List[Dict], jobs_file: str):
    ensure_dir(os.path.dirname(jobs_file))
    with open(jobs_file, "w") as f:
        json.dump(jobs, f, indent=2)
    log(f"✅ Saved {len(jobs)} jobs to {jobs_file}")


def load_jobs(jobs_file: str) -> List[Dict]:
    jobs = _try_load_jobs_any_format(jobs_file)
    log(f"📦 Loaded {len(jobs)} jobs from {jobs_file}")
    # light validation of required keys
    required = {"base_model_nm", "subset", "lr_dir", "model_path", "data_dir", "res_pdir"}
    missing = [i for i,j in enumerate(jobs) if not required.issubset(set(j.keys()))]
    if missing:
        log(f"⚠️ Some jobs missing expected keys at indices: {missing}")
    return jobs


# ─────────────────────────────────────────────────────────────
# Core computation
# ─────────────────────────────────────────────────────────────
def compute_specificity(group: pd.DataFrame) -> pd.Series:
    """specificity = (max(attn) - mean(attn)) / std(attn)"""
    max_val = group["attn_score"].max()
    mean_val = group["attn_score"].mean()
    std_val = group["attn_score"].std()
    specificity = (max_val - mean_val) / std_val if std_val and std_val > 0 else 0.0
    return pd.Series({"specificity": specificity})


def count_motif_hits(seqs: List[str], motifs: List[str]) -> Dict[str, int]:
    """Counts how many sequences contain each motif at least once (presence/absence per seq)."""
    motif_count = {motif: 0 for motif in motifs}
    A = ahocorasick.Automaton()
    for idx, motif in enumerate(motifs):
        A.add_word(motif, (idx, motif))
    A.make_automaton()

    for seq in seqs:
        found = set()
        for _, (_, motif) in A.iter(seq):
            found.add(motif)
        for motif in found:
            motif_count[motif] += 1
    return motif_count


def simple_gapless_score(a: str, b: str) -> int:
    """Return number of matching characters in the best gap-free alignment of a and b."""
    scores = []
    len_a, len_b = len(a), len(b)
    for i in range(-len_b + 1, len_a):
        score = 0
        for j in range(max(0, i), min(len_a, i + len_b)):
            if a[j] == b[j - i]:
                score += 1
        scores.append(score)
    return max(scores) if scores else 0


def merge_similar_motifs(tokens: List[str]) -> Dict[str, str]:
    """Greedy merge: if best gap-free match >= max(min_len-1, min_len//2), map token b -> a."""
    merge_map: Dict[str, str] = {}
    visited = set()
    for i, a in enumerate(tokens):
        if a in visited:
            continue
        for b in tokens[i+1:]:
            if b in visited:
                continue
            score = simple_gapless_score(a, b)
            min_len = min(len(a), len(b))
            req_len = max(min_len - 1, min_len // 2)
            if score >= req_len:
                merge_map[b] = a
                visited.add(b)
    return merge_map


def extract_centered_window(seq: str, start_pos: int, end_pos: int, flank: int = 12) -> str:
    """Return 2*flank window centered at (start_pos, end_pos). None if out of bounds."""
    mid = (int(start_pos) + int(end_pos)) // 2
    s = max(0, mid - flank)
    e = mid + flank
    window_seq = seq[s:e]
    return window_seq if len(window_seq) == 2 * flank else None


def plot_attention_histogram(df: pd.DataFrame, out_png: str):
    """Save a density-like histogram of attn_score split by Label using matplotlib only."""
    plt.figure(figsize=(7, 5))
    for label in sorted(df["Label"].unique()):
        vals = df.loc[df["Label"] == label, "attn_score"].values
        if len(vals) == 0:
            continue
        plt.hist(vals, bins=30, alpha=0.5, density=True, label=f"Label {label}")
    plt.title("Significant Motifs (adj_pval < threshold)")
    plt.xlabel("Attention Score")
    plt.ylabel("Density")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    log(f"📈 Saved histogram: {out_png}")

def process_job(job: Dict,
                fdr_alpha: float = 0.01,
                flank: int = 12,
                min_instances: int = 3,
                last_layer_value: int = None) -> None:
    """Run the full pipeline for a single job dict."""
    res_pdir = job["res_pdir"]
    ensure_dir(res_pdir)
    ensure_dir(os.path.join(res_pdir, "preds"))
    ensure_dir(os.path.join(res_pdir, "motifs"))

    # Load attention + predictions
    attn_path = os.path.join(res_pdir, "2_Attention_motif.parquet")
    pred_path = os.path.join(res_pdir, "preds", "2_combined.csv")

    if not os.path.isfile(attn_path):
        raise FileNotFoundError(f"Missing file: {attn_path}")
    if not os.path.isfile(pred_path):
        raise FileNotFoundError(f"Missing file: {pred_path}")

    df_attn = pd.read_parquet(attn_path)
    df_pred = pd.read_csv(pred_path)

    # Normalize column names
    df_attn = df_attn.rename(columns={"label": "Label", "sequence": "Sequence"})
    # Keep only correct predictions; drop duplicate sequence col if present
    if "sequence" in df_pred.columns and "Sequence" in df_pred.columns:
        df_pred = df_pred.drop(columns=["sequence"])
    df_pred = df_pred[df_pred["true_label"] == df_pred["pred_label"]]

    # Merge on (Sequence, Label)
    merged_df = pd.merge(df_pred, df_attn, on=["Sequence", "Label"], how="inner")
    log(f"Merged shape: {merged_df.shape} (attn x preds)")

    # Compute specificity per (Sequence, layer)
    seq_layer_spec = (
        merged_df.groupby(["Sequence", "layer"])
                 .apply(compute_specificity)
                 .reset_index()
    )

    # Average per layer
    layer_summary = (
        seq_layer_spec.groupby("layer")["specificity"]
                      .mean()
                      .reset_index()
                      .sort_values("specificity", ascending=False)
    )
    log(f"Top layers by specificity:\n{layer_summary.head(5).to_string(index=False)}")

    # Select top-3 + last layer
    top3 = layer_summary.sort_values("specificity", ascending=False).head(3)
    if last_layer_value is None:
        last_layer_value = int(merged_df["layer"].max())  # fallback: max layer present
    layers_to_keep = list(top3["layer"].values) + [last_layer_value]

    filtered_df = merged_df[(merged_df["mask"] == 1) & (merged_df["layer"].isin(layers_to_keep))].copy()
    log(f"Filtered (mask==1 & layers in top3+last): {filtered_df.shape}")

    # Prepare motif set and pos/neg groups
    motif_tokens = filtered_df["token"].astype(str).unique().tolist()
    pos_seqs = filtered_df.loc[filtered_df["Label"] == 1, "Sequence"].astype(str).unique().tolist()
    neg_seqs = filtered_df.loc[filtered_df["Label"] == 0, "Sequence"].astype(str).unique().tolist()

    N = len(pos_seqs) + len(neg_seqs)
    K = len(pos_seqs)
    log(f"Pos seqs: {K}, Neg seqs: {len(neg_seqs)}, Total: {N}, Motifs: {len(motif_tokens)}")

    # Count motif hits (presence per sequence)
    pos_hits = count_motif_hits(pos_seqs, motif_tokens)
    neg_hits = count_motif_hits(neg_seqs, motif_tokens)

    # Compute hypergeometric p-values
    pvals = []
    motif_keys = []
    n_total_dict = {}
    pos_dict = {}
    neg_dict = {}
    hit_type_dict = {}

    for motif in motif_tokens:
        k = pos_hits.get(motif, 0)
        neg = neg_hits.get(motif, 0)
        n = k + neg
        pval = hypergeom.sf(k - 1, N, n, K)  # P(X >= k)

        if k > 0 and neg == 0:
            hit_type = "positive_only"
        elif neg > 0 and k == 0:
            hit_type = "negative_only"
        elif k > 0 and neg > 0:
            hit_type = "both"
        else:
            hit_type = "none"

        pvals.append(pval)
        motif_keys.append(motif)
        n_total_dict[motif] = int(n)
        pos_dict[motif] = int(k)
        neg_dict[motif] = int(neg)
        hit_type_dict[motif] = hit_type

    # FDR correction
    adj_pvals = multipletests(pvals, alpha=fdr_alpha, method="fdr_bh")[1]
    adj_pval_dict = dict(zip(motif_keys, adj_pvals))

    motif_df = filtered_df.copy()
    motif_df["n_total"] = motif_df["token"].map(n_total_dict)
    motif_df["pos_hits"] = motif_df["token"].map(pos_dict)
    motif_df["neg_hits"] = motif_df["token"].map(neg_dict)
    motif_df["pval"] = motif_df["token"].map(dict(zip(motif_keys, pvals)))
    motif_df["adj_pval"] = motif_df["token"].map(adj_pval_dict)
    motif_df["hit_type"] = motif_df["token"].map(hit_type_dict)

    # Significant motifs
    significant_df = motif_df[motif_df["adj_pval"] < fdr_alpha].sort_values("adj_pval").copy()
    sig_csv = os.path.join(res_pdir, "3_sig_motifs.csv")
    significant_df.to_csv(sig_csv, index=False)
    log(f"💾 Wrote: {sig_csv} ({significant_df.shape})")

    # Plot histogram of attention scores for significant motifs
    if not significant_df.empty:
        hist_png = os.path.join(res_pdir, "3_sig_motifs_attn_hist.png")
        plot_attention_histogram(significant_df[["Label", "attn_score"]], hist_png)

    # Merge similar motifs (gap-free score heuristic)
    tokens = significant_df["token"].astype(str).unique().tolist()
    merge_map = merge_similar_motifs(tokens)
    significant_df["merged_token"] = significant_df["token"].map(lambda x: merge_map.get(x, x))

    # Extract centered windows around token spans
    def _row_window(row):
        return extract_centered_window(row["Sequence"], row["start_pos"], row["end_pos"], flank=flank)
    significant_df["window_seq"] = significant_df.apply(_row_window, axis=1)
    significant_df = significant_df.dropna(subset=["window_seq"]).reset_index(drop=True)

    # Keep only positives; require minimum instances per merged motif
    pos_df = significant_df[significant_df["Label"] == 1].copy()
    counts = pos_df["merged_token"].value_counts()
    valid = counts[counts >= min_instances].index.tolist()
    sig_for_pwm = pos_df[pos_df["merged_token"].isin(valid)].drop_duplicates(subset=["token"]).reset_index(drop=True)

    pre_pwm_csv = os.path.join(res_pdir, "4sig_motifs_prePWM.csv")
    sig_for_pwm.to_csv(pre_pwm_csv, index=False)
    log(f"💾 Wrote: {pre_pwm_csv} ({sig_for_pwm.shape})")

    # Export FASTA per label
    fasta_dir = os.path.join(res_pdir, "motifs")
    ensure_dir(fasta_dir)
    for label in significant_df["Label"].unique():
        subset = significant_df[significant_df["Label"] == label].reset_index(drop=True)
        fasta_lines = [f">motif_{i}\n{seq}" for i, seq in enumerate(subset["window_seq"])]
        fasta_path = os.path.join(fasta_dir, f"weblogo_label_{label}.fa")
        with open(fasta_path, "w") as f:
            f.write("\n".join(fasta_lines))
        log(f"💾 Wrote FASTA: {fasta_path} (n={len(subset)})")


# ─────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="DNABERT2 motif enrichment + window extraction pipeline")

    # Primary mode: read jobs from JSON file
    p.add_argument("--jobs_file",
                   default="/data/private/psurana/TSpDNA2/src/files/jobs.json",
                   help="Path to jobs.json. If present and --refresh_jobs is not set, this file is used.")

    # Fallback/regen mode: collect from directories (only used if --refresh_jobs OR jobs_file missing)
    p.add_argument("--base_path", help="Root path containing model runs and RESULT dir")
    p.add_argument("--data_path", help="Root path of mirrored data splits")
    p.add_argument("--length_nm", help="Length folder name (e.g., 3k)")
    p.add_argument("--base_model_nms", help="Comma-separated list (e.g., TSp_vs_nonProm,TSp_vs_genNullseqs)")
    p.add_argument("--subsets", help="Comma-separated subset prefixes (e.g., tspAll_,tspliver_,tsptestis_,tspbrain_)")
    p.add_argument("--lr_prefix", default="lr3e-5", help="LR dir prefix to match (default: lr3e-5)")

    # Control which job to run
    p.add_argument("--job_index", type=int, default=0, help="Index of job to run from jobs list (default: 0)")

    # Rebuild toggle
    p.add_argument("--refresh_jobs", action="store_true",
                   help="Re-scan directories and overwrite --jobs_file with discovered jobs")

    # Pipeline knobs
    p.add_argument("--fdr_alpha", type=float, default=0.01, help="FDR threshold (default: 0.01)")
    p.add_argument("--flank", type=int, default=12, help="Flank size for centered window (default: 12 -> 24bp window)")
    p.add_argument("--min_instances", type=int, default=3, help="Min instances per merged motif for prePWM (default: 3)")
    return p.parse_args()

def main():
    args = parse_args()

    # Decide source of jobs:
    use_jobs_file = (not args.refresh_jobs) and args.jobs_file and os.path.isfile(args.jobs_file)

    if use_jobs_file:
        jobs = load_jobs(args.jobs_file)
    else:
        # Validate required args for collecting from dirs
        needed = ["base_path", "data_path", "length_nm", "base_model_nms", "subsets"]
        missing = [n for n in needed if getattr(args, n) in (None, "")]
        if missing:
            raise SystemExit(
                "To (re)build jobs you must supply: "
                "--base_path --data_path --length_nm --base_model_nms --subsets"
            )

        base_model_nms = [x for x in args.base_model_nms.split(",") if x]
        subsets = [x for x in args.subsets.split(",") if x]

        jobs = collect_jobs(
            base_path=args.base_path,
            base_model_nms=base_model_nms,
            length_nm=args.length_nm,
            subsets=subsets,
            data_path=args.data_path,
            lr_prefix=args.lr_prefix,
        )
        if not jobs:
            raise SystemExit("No jobs found. Check your paths/patterns.")
        # Persist to jobs_file path (either provided or default)
        jobs_file = args.jobs_file or os.path.join(args.base_path, "jobs.json")
        save_jobs(jobs, jobs_file)
        log(f"Rebuilt jobs at: {jobs_file}")

    if args.job_index < 0 or args.job_index >= len(jobs):
        raise SystemExit(f"job_index {args.job_index} out of range (0..{len(jobs)-1})")

    job = jobs[args.job_index]
    log(f"▶ Running job_index={args.job_index}: {job.get('base_model_nm')} | {job.get('subset')} | {job.get('lr_dir')}")

    process_job(
        job=job,
        fdr_alpha=args.fdr_alpha,
        flank=args.flank,
        min_instances=args.min_instances,
    )
    log("✅ Done.")


if __name__ == "__main__":
    main()