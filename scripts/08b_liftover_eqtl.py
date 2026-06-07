#!/usr/bin/env python
# 08b_liftover_eqtl.py — add GRCh37 position (pos37) to each cis_eqtl file via
# pyliftover (hg38->hg19), so eQTL (GRCh38) can be matched to outcome GWAS
# (rsid/GRCh37) by chr:pos. Run in the v11 env (has pyliftover; chain cached).
import os, pandas as pd, pyliftover
P = "/arf/scratch/ycicek/iter_003_mtdna_psy_mr"
IN = f"{P}/data/coloc_inputs"
TISSUES = [
 "Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia",
 "Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9",
 "Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia",
 "Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra"]
lo = pyliftover.LiftOver("hg38", "hg19")
for t in TISSUES:
    fin = f"{IN}/cis_eqtl_{t}.tsv"; fout = f"{IN}/cis_eqtl_lifted_{t}.tsv"
    if not os.path.exists(fin):
        print(f"[skip-missing] {t}"); continue
    e = pd.read_csv(fin, sep="\t")
    if len(e) == 0:
        e.to_csv(fout, sep="\t", index=False); print(f"[empty] {t}"); continue
    e["chrn"] = e["chr"].astype(str).str.replace("chr", "", regex=False)
    u = e[["chrn", "variant_pos"]].drop_duplicates().reset_index(drop=True)
    p37 = []
    for c, p in zip(u["chrn"], u["variant_pos"]):
        r = lo.convert_coordinate("chr" + str(c), int(p) - 1)
        p37.append(r[0][1] + 1 if r else -1)
    u["pos37"] = p37
    e = e.merge(u, on=["chrn", "variant_pos"])
    e = e[e["pos37"] > 0].copy()
    e.to_csv(fout, sep="\t", index=False)
    print(f"[ok] {t}: {len(e)} rows, {e['gene_id'].nunique()} genes lifted")
print("DONE liftover")
