#!/usr/bin/env python
"""
09a_prep_magma_inputs.py
=========================
Prepare MAGMA inputs:
  1. Build SNP-P TSV for each trait (MAGMA gene-level expects "SNP P" + optional N).
  2. Build MitoCarta3.0 .gmt gene set file (Symbol -> pathway membership)
     from Human.MitoCarta3.0.xls sheet "C MitoPathways" (154 pathways)
     + sheet "A Human MitoCarta3.0" (1136 anchor genes -> single "MitoCarta3_all" set).

Outputs:
  data/magma_inputs/<trait>.snp_p.tsv (SNP, P, N)
  data/magma_inputs/mitocarta3.gmt    (pathway_name <tab> Symbol1 Symbol2 ...)
  data/magma_inputs/mitocarta3_set_annot.tsv (set_id <tab> Entrez gene IDs)
                                              [for MAGMA --set-annot]
"""
import os
import sys
import re
import gzip
import time
import pandas as pd

PROJECT = "/arf/scratch/ycicek/iter_003_mtdna_psy_mr"
OUT_DIR = os.path.join(PROJECT, "data", "magma_inputs")
os.makedirs(OUT_DIR, exist_ok=True)

# Map trait -> (file, snp_col, p_col, n_col_or_n_value)
TRAITS = {
    # Exposures
    "Longchamps2022":   dict(path=f"{PROJECT}/data/exposures_harmonized/Longchamps2022.tsv",
                              snp="SNP", p="P", n_col="N"),
    "Chong2022":        dict(path=f"{PROJECT}/data/exposures_harmonized/Chong2022.tsv",
                              snp="SNP", p="P", n_col="N"),
    "Gupta2023_raw":    dict(path=f"{PROJECT}/data/exposures_harmonized/Gupta2023_raw.tsv",
                              snp="SNP", p="P", n_col="N"),
    "Gupta2023_adj":    dict(path=f"{PROJECT}/data/exposures_harmonized/Gupta2023_adjusted.tsv",
                              snp="SNP", p="P", n_col="N"),
    # Outcomes
    "ASD_Grove2019":    dict(path=f"{PROJECT}/data/outcomes/ASD_Grove2019.gz",
                              snp="SNP", p="P", n_value=46350),  # Grove 2019: 18,381 + 27,969
    "ADHD_Demontis2023": dict(path=f"{PROJECT}/data/outcomes/ADHD_Demontis2023.meta.gz",
                              snp="SNP", p="P", n_value=225534),  # 38691 + 186843
    "CDG3_F3_Neurodev": dict(path=f"{PROJECT}/data/outcomes/CDG3_F3_Neurodev.tsv.gz",
                              snp="SNP", p="P", n_value=500000),  # CDG3 effective N
}


def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def read_gwas(meta):
    log(f"  Reading: {meta['path']}")
    # pandas 3.0 prohibits low_memory with sep=None (python engine).
    # Detect compression and pick separator manually for speed + memory safety.
    p = meta["path"]
    comp = "gzip" if p.endswith(".gz") else "infer"
    # Try tab first (most), then whitespace
    try:
        df = pd.read_csv(p, sep="\t", compression=comp, low_memory=False)
        if df.shape[1] == 1:
            raise ValueError("single col — wrong sep")
    except Exception:
        df = pd.read_csv(p, sep=r"\s+", engine="python", compression=comp)
    # Standardize SNP column
    if meta["snp"] not in df.columns:
        cand = [c for c in df.columns if c.upper() in ("SNP", "RSID", "RS_ID", "MARKER")]
        if cand:
            df[meta["snp"]] = df[cand[0]]
    if meta["p"] not in df.columns:
        cand = [c for c in df.columns if c.upper() in ("P", "PVAL", "P_VALUE", "PVALUE", "P-VALUE")]
        if cand:
            df[meta["p"]] = df[cand[0]]
    df = df[[meta["snp"], meta["p"]] + (
        [meta["n_col"]] if "n_col" in meta and meta["n_col"] in df.columns else []
    )].copy()
    df.columns = ["SNP", "P"] + (["N"] if "n_col" in meta and meta["n_col"] in df.columns else [])
    df = df.dropna(subset=["SNP", "P"])
    df = df[df["SNP"].astype(str).str.startswith("rs")]
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df.dropna(subset=["P"])
    df = df[(df["P"] > 0) & (df["P"] <= 1)]
    if "N" not in df.columns:
        df["N"] = meta.get("n_value", 100000)
    df = df.drop_duplicates(subset=["SNP"])
    return df


def main():
    log("=== Build SNP-P input files for MAGMA ===")
    for trait, meta in TRAITS.items():
        out_path = os.path.join(OUT_DIR, f"{trait}.snp_p.tsv")
        if os.path.exists(out_path):
            log(f"  [skip] {trait} -> exists")
            continue
        try:
            df = read_gwas(meta)
        except Exception as e:
            log(f"  [ERR] {trait}: {e}")
            continue
        df.to_csv(out_path, sep="\t", index=False)
        log(f"  [write] {out_path}: {len(df):,} SNPs")

    # --- Build MitoCarta gmt + set_annot ---
    mc_xls = f"{PROJECT}/data/magma_aux/Human.MitoCarta3.0.xls"
    log(f"\n=== Build MitoCarta gene-set files from {mc_xls} ===")
    anchor = pd.read_excel(mc_xls, sheet_name="A Human MitoCarta3.0")
    pw = pd.read_excel(mc_xls, sheet_name="C MitoPathways")
    log(f"  Anchor genes (sheet A): {len(anchor)}")
    log(f"  Pathways (sheet C): {len(pw)}")

    # Symbol -> Entrez map from sheet A (HumanGeneID = Entrez, some entries are
    # Ensembl IDs for mtDNA-encoded genes -> drop these for MAGMA Entrez-based
    # set-annot; gmt symbol file retains all entries via raw symbol list below.)
    anchor["_EntrezID"] = pd.to_numeric(anchor["HumanGeneID"], errors="coerce")
    valid_entrez = anchor.dropna(subset=["_EntrezID"]).copy()
    valid_entrez["_EntrezID"] = valid_entrez["_EntrezID"].astype(int)
    sym_to_entrez = dict(zip(valid_entrez["Symbol"].astype(str).str.strip(),
                              valid_entrez["_EntrezID"]))

    # Build sets
    sets = {}
    # All MitoCarta3.0 anchor (1 set; Entrez-only genes for set_annot, all symbols for gmt)
    all_mito_entrez = sorted(valid_entrez["_EntrezID"].unique())
    all_mito_symbols = sorted(anchor["Symbol"].dropna().astype(str).str.strip().unique())
    sets["MitoCarta3_all_anchor"] = dict(entrez=all_mito_entrez, symbols=all_mito_symbols)

    # 154 MitoPathways (sheet C)
    pw = pw.dropna(subset=["MitoPathway", "Genes"])
    for _, r in pw.iterrows():
        name = str(r["MitoPathway"]).strip()
        # Genes column is comma-sep symbols, sometimes with spaces / mixed delim
        raw = str(r["Genes"])
        symbols = re.split(r"[,;\s]+", raw)
        symbols = [s.strip() for s in symbols if s.strip()]
        entrez = [sym_to_entrez[s] for s in symbols if s in sym_to_entrez]
        # Sanitize name for filenames
        safe_name = re.sub(r"\W+", "_", name).strip("_")
        set_id = f"MITO_PW__{safe_name}"
        sets[set_id] = dict(entrez=sorted(set(entrez)), symbols=sorted(set(symbols)))

    log(f"  Total sets built: {len(sets)}")

    # Symbol-based .gmt (for downstream HGNC -> Ensembl conversions if needed)
    gmt_path = os.path.join(OUT_DIR, "mitocarta3.gmt")
    with open(gmt_path, "w") as fp:
        for sid, payload in sets.items():
            line = "\t".join([sid, "MitoCarta3.0"] + payload["symbols"])
            fp.write(line + "\n")
    log(f"  [write] {gmt_path}")

    # MAGMA --set-annot format: <set_id> <entrez_id1> <entrez_id2> ...
    # (one set per line, Entrez IDs separated by whitespace)
    set_annot_path = os.path.join(OUT_DIR, "mitocarta3_set_annot.tsv")
    with open(set_annot_path, "w") as fp:
        for sid, payload in sets.items():
            if not payload["entrez"]:
                continue
            line = " ".join([sid] + [str(e) for e in payload["entrez"]])
            fp.write(line + "\n")
    log(f"  [write] {set_annot_path}")

    log("DONE -- prep_magma_inputs")


if __name__ == "__main__":
    main()
