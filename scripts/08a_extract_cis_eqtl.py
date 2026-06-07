#!/usr/bin/env python
"""
08a_extract_cis_eqtl.py
=======================
Helper for Task 2 (COLOC). For each lead SNP × ±500kb region, extract all
cis-eQTL associations from GTEx v10 .parquet across 13 brain tissues.

Reads:
  - data/coloc_inputs/lead_loci.tsv    (SNP, CHR, POS, region_id) — from R prep
  - /arf/scratch/ycicek/h2_paper2/data/gtex_v10_eqtl_brain/Brain_*.v10.eQTLs.signif_pairs.parquet
    NOTE: GTEx v10 parquet has all *significant* nominal pairs only.

Writes:
  - data/coloc_inputs/cis_eqtl_<TISSUE>.tsv (region_id, gene_id, variant_id,
    chr, pos, ref, alt, pval_nominal, slope, slope_se, ma_count, af)

Tissue list (13 GTEx brain regions):
  Brain_Amygdala, Brain_Anterior_cingulate_cortex_BA24,
  Brain_Caudate_basal_ganglia, Brain_Cerebellar_Hemisphere, Brain_Cerebellum,
  Brain_Cortex, Brain_Frontal_Cortex_BA9, Brain_Hippocampus,
  Brain_Hypothalamus, Brain_Nucleus_accumbens_basal_ganglia,
  Brain_Putamen_basal_ganglia, Brain_Spinal_cord_cervical_c-1,
  Brain_Substantia_nigra
"""
import os
import sys
import time
import pandas as pd
import pyarrow.parquet as pq

GTEX_DIR = "/arf/scratch/ycicek/h2_paper2/data/gtex_v10_eqtl_brain"
BRAIN_TISSUES = [
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
]
REGION_HALF_WIDTH = 500_000  # +/- 500 kb


def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def main(loci_tsv: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)
    loci = pd.read_csv(loci_tsv, sep="\t")
    required = {"region_id", "CHR", "POS_START", "POS_END"}
    missing = required - set(loci.columns)
    if missing:
        sys.exit(f"FATAL: lead_loci.tsv missing columns: {missing}")
    log(f"Loaded {len(loci)} regions from {loci_tsv}")
    # Standardize chromosome notation for GRCh38 -- GTEx v10 uses 'chr1'-'chr22'
    loci["chr_str"] = loci["CHR"].astype(str).apply(lambda x: x if x.startswith("chr") else f"chr{x}")

    for tissue in BRAIN_TISSUES:
        out_path = os.path.join(outdir, f"cis_eqtl_{tissue}.tsv")
        if os.path.exists(out_path):
            log(f"  [skip] {tissue} -> {out_path} exists")
            continue
        parq = os.path.join(GTEX_DIR, f"{tissue}.v10.eQTLs.signif_pairs.parquet")
        if not os.path.exists(parq):
            log(f"  [WARN] parquet missing for {tissue}: {parq}")
            continue
        log(f"  [load] {tissue} parquet ({os.path.getsize(parq) / 1e6:.1f} MB)")
        # Stream by row group to limit memory; reads ~5-50M rows per tissue
        df = pq.read_table(parq).to_pandas()
        log(f"    rows: {len(df):,}, cols: {list(df.columns)[:15]}")
        # GTEx v10 parquet schema in this build: only variant_id like
        # 'chr1_12345_A_G_b38'. Derive chr + variant_pos by splitting.
        if "chr" not in df.columns or "variant_pos" not in df.columns:
            if "variant_id" in df.columns:
                log(f"    [info] parsing chr/variant_pos from variant_id")
                # variant_id format: chrN_pos_ref_alt_build
                parts = df["variant_id"].astype(str).str.split("_", n=4, expand=True)
                if parts.shape[1] >= 4:
                    df["chr"] = parts[0]
                    df["variant_pos"] = pd.to_numeric(parts[1], errors="coerce").astype("Int64")
                    df["ref"] = parts[2]
                    df["alt"] = parts[3]
                else:
                    log(f"    [ERR] variant_id split unexpected: {parts.shape}")
                    continue
            else:
                log(f"    [ERR] no chr/variant_pos and no variant_id either")
                continue
        # Build per-region subset
        keep_frames = []
        for _, row in loci.iterrows():
            sub = df[(df["chr"] == row["chr_str"]) &
                     (df["variant_pos"] >= row["POS_START"]) &
                     (df["variant_pos"] <= row["POS_END"])].copy()
            if len(sub) == 0:
                continue
            sub["region_id"] = row["region_id"]
            keep_frames.append(sub)
        if not keep_frames:
            log(f"    [WARN] no rows match any region in {tissue}")
            # Still write empty header so R doesn't break
            empty = pd.DataFrame(columns=["region_id", "gene_id", "variant_id",
                                          "chr", "variant_pos", "ref", "alt",
                                          "pval_nominal", "slope", "slope_se",
                                          "ma_count", "af", "rs_id_dbSNP155_GRCh38p13"])
            empty.to_csv(out_path, sep="\t", index=False)
            continue
        out = pd.concat(keep_frames, ignore_index=True)
        # Keep only useful cols (some may not exist in this schema)
        candidate_cols = ["region_id", "gene_id", "variant_id", "chr", "variant_pos",
                          "ref", "alt", "pval_nominal", "slope", "slope_se",
                          "ma_count", "af", "rs_id_dbSNP155_GRCh38p13"]
        keep_cols = [c for c in candidate_cols if c in out.columns]
        out = out[keep_cols]
        out.to_csv(out_path, sep="\t", index=False)
        n_genes = out["gene_id"].nunique()
        log(f"    [write] {out_path}: {len(out):,} rows, {n_genes} genes")
        del df, out, keep_frames

    log("DONE — extract_cis_eqtl complete")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: 08a_extract_cis_eqtl.py <lead_loci.tsv> <outdir>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
