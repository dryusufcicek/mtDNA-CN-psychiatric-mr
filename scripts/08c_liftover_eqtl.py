#!/usr/bin/env python
"""
08c_liftover_eqtl.py
====================
Liftover GTEx v10 eqtl variant positions from GRCh38 to GRCh37,
adding chr_grch37 and pos_grch37 columns to cis_eqtl_<TISSUE>.tsv,
writing cis_eqtl_<TISSUE>.lifted.tsv

Reads chain file from /arf/home/ycicek/.pyliftover/hg38ToHg19.over.chain.gz
"""
import os, sys, time
import pandas as pd
from pyliftover import LiftOver

TISSUES = [
    "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus", "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra"
]

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def main(coloc_dir):
    lo = LiftOver('hg38', 'hg19')
    log("LiftOver loaded hg38ToHg19")
    for tissue in TISSUES:
        in_p = os.path.join(coloc_dir, f"cis_eqtl_{tissue}.tsv")
        out_p = os.path.join(coloc_dir, f"cis_eqtl_{tissue}.lifted.tsv")
        if not os.path.exists(in_p):
            log(f"  [skip] {tissue} input missing")
            continue
        if os.path.exists(out_p):
            log(f"  [skip] {tissue} already lifted")
            continue
        df = pd.read_csv(in_p, sep="\t")
        if len(df) == 0:
            log(f"  [skip] {tissue} empty")
            df["chr_grch37"] = pd.Series(dtype="object")
            df["pos_grch37"] = pd.Series(dtype="Int64")
            df.to_csv(out_p, sep="\t", index=False)
            continue
        # df.chr is e.g. 'chr1'
        df["chr_clean"] = df["chr"].astype(str).str.replace("^chr", "", regex=True)
        chr_g37 = []
        pos_g37 = []
        miss = 0
        for chr_n, pos in zip(df["chr_clean"], df["variant_pos"]):
            try:
                out = lo.convert_coordinate(f"chr{chr_n}", int(pos) - 1)
            except Exception:
                out = None
            if out:
                chr_g37.append(out[0][0].replace("chr", ""))
                pos_g37.append(out[0][1] + 1)  # 1-based
            else:
                chr_g37.append(None)
                pos_g37.append(None)
                miss += 1
        df["chr_grch37"] = chr_g37
        df["pos_grch37"] = pos_g37
        df = df.dropna(subset=["chr_grch37", "pos_grch37"])
        df.to_csv(out_p, sep="\t", index=False)
        log(f"  [{tissue}] lifted {len(df):,} rows, {miss} unmappable")
    log("DONE liftover")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: 08c_liftover_eqtl.py <coloc_dir>")
        sys.exit(1)
    main(sys.argv[1])
