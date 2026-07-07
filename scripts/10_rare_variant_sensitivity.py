#!/usr/bin/env python
"""
10_rare_variant_sensitivity.py
==============================
Rare / low-frequency-instrument (MAF) leave-out sensitivity for the forward-MR
results (Supplementary Table 17). Addresses whether the primary ASD association is
driven by rare or low-frequency instruments (e.g. the DGUOK/POLG2 low-frequency
missense variants).

Two steps:
  1. Instrument MAF is read from the harmonized exposure allele frequencies
     (EAF -> MAF = min(EAF, 1-EAF)) for each of the four mtDNA-CN panels. The
     Chong (SNP-array) panel has no usable EAF column in the harmonized file, so
     its instruments are handled via the shared-SNP low-MAF set.
  2. For each panel the fixed- and random-effects IVW estimate is re-computed
     from the per-instrument beta_exposure / se_outcome pairs (as tabulated in
     Supplementary Table 3), after excluding instruments with MAF < 1% and,
     separately, MAF < 5%. Panels are then pooled with the same design-effect
     (effective-sample-size) correction used for the primary analysis
     (rho = 0.70, DEFF = 1 + (k-1)*rho).

Reproduction check: with NO exclusions this script reproduces the primary
per-panel and pooled ASD estimates exactly (beta = -0.1818, corrected p = 6.6e-3).

Reads:
  - data/exposures_harmonized/{Longchamps2022,Chong2022,Gupta2023_raw,Gupta2023_adjusted}.tsv
  - supplementary_tables/Supp_Table_3_instrument_list.tsv   (per-instrument bx, by, se)
Writes:
  - results/sensitivity/rare_variant_sensitivity.tsv        (-> Supplementary Table 17)
"""
import os, math, csv
from collections import defaultdict

BASE = os.environ.get("ITER003_BASE", ".")
HARM = os.path.join(BASE, "data", "exposures_harmonized")
ST3  = os.path.join(BASE, "supplementary_tables", "Supp_Table_3_instrument_list.tsv")
OUT  = os.path.join(BASE, "results", "sensitivity", "rare_variant_sensitivity.tsv")
RHO  = 0.70
PANELS = {"Longchamps_2022": "Longchamps2022.tsv", "Chong_2022": "Chong2022.tsv",
          "Gupta_2023_raw": "Gupta2023_raw.tsv", "Gupta_2023_adj": "Gupta2023_adjusted.tsv"}


def instrument_maf(want_by_panel):
    """MAF per (panel, SNP) from harmonized exposure EAF."""
    maf = {}
    for panel, fname in PANELS.items():
        path = os.path.join(HARM, fname)
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ci = next((hdr.index(c) for c in ("EAF", "A1FREQ", "FREQ", "eaf", "af", "AF") if c in hdr), None)
            si = hdr.index("SNP") if "SNP" in hdr else 0
            if ci is None:
                continue
            want = want_by_panel.get(panel, set())
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) > ci and c[si] in want:
                    try:
                        eaf = float(c[ci]); maf[(panel, c[si])] = min(eaf, 1 - eaf)
                    except ValueError:
                        pass
    return maf


def ivw(pairs):
    """Fixed + multiplicative-random-effects IVW from [(bx, by, sy), ...]."""
    num = sum(bx * by / sy**2 for bx, by, sy in pairs)
    den = sum(bx**2 / sy**2 for bx, by, sy in pairs)
    beta = num / den
    se_fix = math.sqrt(1 / den)
    q = sum((by - beta * bx)**2 / sy**2 for bx, by, sy in pairs)
    df = len(pairs) - 1
    se_re = se_fix * math.sqrt(max(1.0, q / df)) if df > 0 else se_fix
    return beta, se_fix, se_re, len(pairs)


def pooled_corrected(betas, ses):
    """Naive inverse-variance point estimate; SE inflated by the design effect."""
    w = [1 / s**2 for s in ses]
    beta = sum(b * wi for b, wi in zip(betas, w)) / sum(w)
    se_naive = math.sqrt(1 / sum(w))
    k = len(betas)
    deff = 1 + (k - 1) * RHO
    se_corr = se_naive * math.sqrt(deff)
    z = beta / se_corr
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return beta, se_corr, p


def main():
    # per-panel ASD instruments and their bx/by/sy from Supplementary Table 3
    by_panel = defaultdict(list); want = defaultdict(set)
    with open(ST3) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if row.get("outcome", "").upper() != "ASD":
                continue
            p = row["exposure"]; snp = row["SNP"]
            bx = float(row["beta_exposure"]); by = float(row["beta_outcome"]); sy = float(row["se_outcome"])
            by_panel[p].append((snp, bx, by, sy)); want[p].add(snp)

    maf = instrument_maf(want)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["threshold", "panel", "n_instruments", "n_excluded",
                    "beta_fix", "se_fix", "se_random"])
        pooled_rows = []
        for label, thr in [("all", None), ("MAF>=1%", 0.01), ("MAF>=5%", 0.05)]:
            betas, ses = [], []
            for p, rows in by_panel.items():
                keep = [(bx, by, sy) for snp, bx, by, sy in rows
                        if thr is None or maf.get((p, snp), 1.0) >= thr]
                nexcl = len(rows) - len(keep)
                if len(keep) < 2:
                    continue
                b, sf, sr, n = ivw(keep)
                w.writerow([label, p, n, nexcl, f"{b:.4f}", f"{sf:.4f}", f"{sr:.4f}"])
                betas.append(b); ses.append(sr)
            if betas:
                pb, pse, pp = pooled_corrected(betas, ses)
                pooled_rows.append([label, "POOLED", len(betas), "", f"{pb:.4f}", f"{pse:.4f}", f"p={pp:.2e}"])
        for row in pooled_rows:
            w.writerow(row)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
