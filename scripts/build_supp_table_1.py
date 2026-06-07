#!/usr/bin/env python3
"""
build_supp_table_1.py  — Supplementary Table 1: GWAS metadata
iter_003 mtDNA-CN psychiatric MR

Computes per-file per-variant N statistics from the actual exposure/outcome
GWAS summary statistics on TRUBA, then emits a TSV row per dataset.

ZERO-FABRICATION rule:
  - every numeric column comes from the file
  - if a column cannot be computed -> "N/A — author to supply" with reason
"""

import gzip
import os
import sys
import io
from pathlib import Path

EXP_DIR = "/arf/scratch/ycicek/iter_003_mtdna_psy_mr/data/exposures_harmonized"
OUT_DIR = "/arf/scratch/ycicek/iter_003_mtdna_psy_mr/data/outcomes"
PROJECT_ROOT = "/arf/scratch/ycicek/iter_003_mtdna_psy_mr"

OUTPUT_TSV = "/arf/scratch/ycicek/iter_003_mtdna_psy_mr/manuscript_draft/supplementary/Supp_Table_1_GWAS_metadata.tsv"


def open_text(path):
    """Open .gz or plain text transparently."""
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def median(seq):
    s = sorted(seq)
    n = len(s)
    if n == 0:
        return None
    if n % 2 == 1:
        return float(s[n // 2])
    return 0.5 * (s[n // 2 - 1] + s[n // 2])


def stream_columns(path, col_indices_dict, sep=None, comment_prefix=None,
                   skip_vcf_meta=False, header_must_have=None):
    """
    Stream a (possibly gzipped) summary stats file and return summary stats
    for the requested column indices using a FREQUENCY-TABLE approach
    (memory-bounded: only distinct values stored, not every record).

    For Nca/Nco/N/Neff style columns, the number of distinct values is at most
    a few thousand even across 16M variants (typically <100). This keeps RAM
    flat regardless of file size.

    Returns: dict {logical_name -> {"min": ..., "median": ..., "max": ..., "n_records": int}}
             plus key "_n_lines" = total non-header non-meta lines counted
    """
    stats = {}
    n_lines = 0
    col_index = {}
    # value -> count;   memory-bounded by distinct-value count
    freq = {name: {} for name in col_indices_dict}

    print(f"  scanning {path}", file=sys.stderr, flush=True)
    with open_text(path) as fh:
        # Find header
        header_line = None
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            if skip_vcf_meta and line.startswith("##"):
                continue
            header_line = line
            break
        if header_line is None:
            raise RuntimeError(f"No header found in {path}")

        header_clean = header_line.lstrip("#")
        if sep is None:
            header_cols = header_clean.split()
        else:
            header_cols = header_clean.split(sep)

        for logical, header_name in col_indices_dict.items():
            if header_name not in header_cols:
                col_index[logical] = None
            else:
                col_index[logical] = header_cols.index(header_name)

        if header_must_have is not None:
            for must in header_must_have:
                if must not in header_cols:
                    raise RuntimeError(
                        f"Header missing required column '{must}' in {path}: {header_cols[:20]}"
                    )

        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if skip_vcf_meta and line.startswith("##"):
                continue
            if sep is None:
                fields = line.split()
            else:
                fields = line.split(sep)
            n_lines += 1
            for logical, idx in col_index.items():
                if idx is None or idx >= len(fields):
                    continue
                raw_val = fields[idx]
                if raw_val == "NA" or raw_val == "":
                    continue
                try:
                    v = float(raw_val)
                except ValueError:
                    continue
                # quantize to integer-ish for thinned N columns; otherwise keep float
                key = v
                freq[logical][key] = freq[logical].get(key, 0) + 1

    for logical, fmap in freq.items():
        if not fmap:
            stats[logical] = None
            continue
        keys_sorted = sorted(fmap.keys())
        total = sum(fmap.values())
        # Compute exact median by walking cumulative count
        target_lo = (total + 1) // 2  # 1-indexed lower middle
        target_hi = (total // 2) + 1  # 1-indexed upper middle (handles even n)
        cum = 0
        med_lo = med_hi = None
        for k in keys_sorted:
            cum += fmap[k]
            if med_lo is None and cum >= target_lo:
                med_lo = k
            if med_hi is None and cum >= target_hi:
                med_hi = k
                break
        median_val = 0.5 * (med_lo + med_hi)
        stats[logical] = {
            "min": keys_sorted[0],
            "median": median_val,
            "max": keys_sorted[-1],
            "n_records": total,
        }
    stats["_n_lines"] = n_lines
    return stats


def fmt_int(x):
    if x is None:
        return "NA"
    return f"{int(round(x))}"


def fmt_eff_n_cc(nca, nco):
    """Effective N for case-control GWAS:  4*Nca*Nco/(Nca+Nco)."""
    if nca is None or nco is None or (nca + nco) == 0:
        return "NA"
    return f"{int(round(4.0 * nca * nco / (nca + nco)))}"


# =====================================================
# Per-file specification (label, file, columns, source citation, etc.)
# =====================================================

REL_ROOT = "data/"  # all paths reported as "data/<...>"

ROWS = []


def add_row(**kwargs):
    ROWS.append(kwargs)


# --------------- EXPOSURES (4) ---------------

def process_exposure(label, fname, citation, ancestry, headline_n, notes):
    print(f"[{label}] starting", file=sys.stderr, flush=True)
    full = os.path.join(EXP_DIR, fname)
    # Exposure harmonized files: SNP CHR BP A1 A2 BETA SE P EAF N INFO exposure
    stats = stream_columns(
        full,
        col_indices_dict={"N": "N"},
        header_must_have=["N"],
    )
    n_stats = stats["N"]
    if n_stats is None:
        n_min = n_med = n_max = None
    else:
        n_min = int(round(n_stats["min"]))
        n_med = int(round(n_stats["median"]))
        n_max = int(round(n_stats["max"]))
    total_snps = stats["_n_lines"]
    add_row(
        dataset=label,
        role="exposure",
        path=f"{REL_ROOT}exposures_harmonized/{fname}",
        citation=citation,
        ancestry=ancestry,
        headline_n=headline_n,
        nca_min="NA",
        nca_median="NA",
        nca_max="NA",
        nco_min="NA",
        nco_median="NA",
        nco_max="NA",
        n_min=fmt_int(n_min),
        n_median=fmt_int(n_med),
        n_max=fmt_int(n_max),
        n_eff=fmt_int(n_max),  # continuous exposure: effective N = total N
        total_snps=str(total_snps),
        file_format="harmonized TSV (post-pipeline; SNP CHR BP A1 A2 BETA SE P EAF N INFO)",
        notes=notes,
    )


def process_cc_pgc_daner(label, fname, citation, ancestry, headline_n, role, notes,
                          nca_col="Nca", nco_col="Nco"):
    """PGC daner format with Nca/Nco per variant + FRQ_A_X/FRQ_U_X headers."""
    print(f"[{label}] starting", file=sys.stderr, flush=True)
    full = os.path.join(OUT_DIR, fname)
    stats = stream_columns(
        full,
        col_indices_dict={"Nca": nca_col, "Nco": nco_col},
        header_must_have=[nca_col, nco_col],
    )
    nca = stats["Nca"]
    nco = stats["Nco"]
    nca_min = int(round(nca["min"])) if nca else None
    nca_med = int(round(nca["median"])) if nca else None
    nca_max = int(round(nca["max"])) if nca else None
    nco_min = int(round(nco["min"])) if nco else None
    nco_med = int(round(nco["median"])) if nco else None
    nco_max = int(round(nco["max"])) if nco else None
    eff_n = fmt_eff_n_cc(nca_max, nco_max) if nca_max and nco_max else "NA"
    add_row(
        dataset=label,
        role=role,
        path=f"{REL_ROOT}outcomes/{fname}",
        citation=citation,
        ancestry=ancestry,
        headline_n=headline_n,
        nca_min=fmt_int(nca_min),
        nca_median=fmt_int(nca_med),
        nca_max=fmt_int(nca_max),
        nco_min=fmt_int(nco_min),
        nco_median=fmt_int(nco_med),
        nco_max=fmt_int(nco_max),
        n_min="NA",
        n_median="NA",
        n_max="NA",
        n_eff=eff_n,
        total_snps=str(stats["_n_lines"]),
        file_format="PGC daner (case-control, FRQ_A_X / FRQ_U_X headers; per-variant Nca/Nco)",
        notes=notes,
    )


def process_grove_no_n(label, fname, citation, ancestry, headline_n, role, notes):
    """ASD Grove 2019 - no per-variant N, no Nca/Nco columns at all."""
    print(f"[{label}] starting", file=sys.stderr, flush=True)
    full = os.path.join(OUT_DIR, fname)
    stats = stream_columns(
        full,
        col_indices_dict={"OR": "OR"},
        header_must_have=["OR", "SE", "P"],
    )
    add_row(
        dataset=label,
        role=role,
        path=f"{REL_ROOT}outcomes/{fname}",
        citation=citation,
        ancestry=ancestry,
        headline_n=headline_n,
        nca_min="N/A — author to supply (Grove 2019 release omits per-variant Nca)",
        nca_median="N/A",
        nca_max="N/A",
        nco_min="N/A — author to supply (Grove 2019 release omits per-variant Nco)",
        nco_median="N/A",
        nco_max="N/A",
        n_min="NA",
        n_median="NA",
        n_max="NA",
        n_eff="NA",
        total_snps=str(stats["_n_lines"]),
        file_format="PGC daner-style (CHR SNP BP A1 A2 INFO OR SE P; no Nca/Nco/N)",
        notes=notes,
    )


def process_neg_ctrl(label, fname, citation, ancestry, headline_n, notes):
    """Yengo et al. BMI/height - continuous N column."""
    print(f"[{label}] starting", file=sys.stderr, flush=True)
    full = os.path.join(OUT_DIR, fname)
    stats = stream_columns(
        full,
        col_indices_dict={"N": "N"},
        header_must_have=["N"],
    )
    n_stats = stats["N"]
    n_min = int(round(n_stats["min"])) if n_stats else None
    n_med = int(round(n_stats["median"])) if n_stats else None
    n_max = int(round(n_stats["max"])) if n_stats else None
    add_row(
        dataset=label,
        role="negative control",
        path=f"{REL_ROOT}outcomes/{fname}",
        citation=citation,
        ancestry=ancestry,
        headline_n=headline_n,
        nca_min="NA",
        nca_median="NA",
        nca_max="NA",
        nco_min="NA",
        nco_median="NA",
        nco_max="NA",
        n_min=fmt_int(n_min),
        n_median=fmt_int(n_med),
        n_max=fmt_int(n_max),
        n_eff=fmt_int(n_max),  # continuous trait
        total_snps=str(stats["_n_lines"]),
        file_format="Yengo 2022 GIANT (continuous trait; per-variant N)",
        notes=notes,
    )


def process_cdg3(label, fname, citation, ancestry, headline_n, notes):
    """Grotzinger 2022 cross-disorder factor GWAS - no N column."""
    print(f"[{label}] starting", file=sys.stderr, flush=True)
    full = os.path.join(OUT_DIR, fname)
    stats = stream_columns(
        full,
        col_indices_dict={"BETA": "BETA"},
        header_must_have=["BETA", "SE", "P"],
    )
    add_row(
        dataset=label,
        role="CDG3 factor",
        path=f"{REL_ROOT}outcomes/{fname}",
        citation=citation,
        ancestry=ancestry,
        headline_n=headline_n,
        nca_min="N/A — CDG3 GWAS is genomic-SEM factor, no per-variant Nca",
        nca_median="N/A",
        nca_max="N/A",
        nco_min="N/A — CDG3 GWAS is genomic-SEM factor, no per-variant Nco",
        nco_median="N/A",
        nco_max="N/A",
        n_min="N/A — no per-variant N column in source release",
        n_median="N/A",
        n_max="N/A",
        n_eff="N/A — author to supply from source paper (Grotzinger 2022)",
        total_snps=str(stats["_n_lines"]),
        file_format="Grotzinger 2022 CDG3 (SNP CHR BP MAF A1 A2 BETA SE P Q_P; no N column)",
        notes=notes,
    )


def process_ppd_2023(label, fname, citation, ancestry, headline_n, notes):
    """PPD 2023 - TotalSampleSize column, no Nca/Nco."""
    print(f"[{label}] starting", file=sys.stderr, flush=True)
    full = os.path.join(OUT_DIR, fname)
    stats = stream_columns(
        full,
        col_indices_dict={"TotalSampleSize": "TotalSampleSize"},
        header_must_have=["TotalSampleSize"],
    )
    n_stats = stats["TotalSampleSize"]
    if n_stats is None:
        n_min = n_med = n_max = None
    else:
        n_min = int(round(n_stats["min"]))
        n_med = int(round(n_stats["median"]))
        n_max = int(round(n_stats["max"]))
    add_row(
        dataset=label,
        role="secondary outcome",
        path=f"{REL_ROOT}outcomes/{fname}",
        citation=citation,
        ancestry=ancestry,
        headline_n=headline_n,
        nca_min="N/A — author to supply (PPD 2023 release reports TotalSampleSize only)",
        nca_median="N/A",
        nca_max="N/A",
        nco_min="N/A",
        nco_median="N/A",
        nco_max="N/A",
        n_min=fmt_int(n_min),
        n_median=fmt_int(n_med),
        n_max=fmt_int(n_max),
        n_eff="N/A — author to supply (no Nca/Nco)",
        total_snps=str(stats["_n_lines"]),
        file_format="METAL meta-analysis output (CHR BP SNP A1 A2 ... BETA SE P ... TotalSampleSize)",
        notes=notes,
    )


def process_sud_2023(label, fname, citation, ancestry, headline_n, notes):
    """Hatoum SUD 2023 - SNP Chr BP A1 A2 Beta P; no SE, no N."""
    print(f"[{label}] starting", file=sys.stderr, flush=True)
    full = os.path.join(OUT_DIR, fname)
    stats = stream_columns(
        full,
        col_indices_dict={"Beta": "Beta"},
        header_must_have=["Beta", "P"],
    )
    add_row(
        dataset=label,
        role="secondary outcome (excluded from forward MR)",
        path=f"{REL_ROOT}outcomes/{fname}",
        citation=citation,
        ancestry=ancestry,
        headline_n=headline_n,
        nca_min="N/A — author to supply (Hatoum 2023 release omits per-variant Nca)",
        nca_median="N/A",
        nca_max="N/A",
        nco_min="N/A",
        nco_median="N/A",
        nco_max="N/A",
        n_min="N/A — no per-variant N column",
        n_median="N/A",
        n_max="N/A",
        n_eff="N/A — author to supply from source paper",
        total_snps=str(stats["_n_lines"]),
        file_format="custom TSV (SNP Chr BP A1 A2 Beta P; no SE, no N) — excluded from forward MR",
        notes=notes,
    )


def process_an_vcf(label, fname, citation, ancestry, headline_n, notes):
    """AN PGC VCF v1.0 - NCAS / NCON per variant."""
    print(f"[{label}] starting", file=sys.stderr, flush=True)
    full = os.path.join(OUT_DIR, fname)
    stats = stream_columns(
        full,
        col_indices_dict={"NCAS": "NCAS", "NCON": "NCON"},
        skip_vcf_meta=True,
        header_must_have=["NCAS", "NCON"],
    )
    nca = stats["NCAS"]
    nco = stats["NCON"]
    nca_min = int(round(nca["min"])) if nca else None
    nca_med = int(round(nca["median"])) if nca else None
    nca_max = int(round(nca["max"])) if nca else None
    nco_min = int(round(nco["min"])) if nco else None
    nco_med = int(round(nco["median"])) if nco else None
    nco_max = int(round(nco["max"])) if nco else None
    eff_n = fmt_eff_n_cc(nca_max, nco_max) if nca_max and nco_max else "NA"
    add_row(
        dataset=label,
        role="secondary outcome",
        path=f"{REL_ROOT}outcomes/{fname}",
        citation=citation,
        ancestry=ancestry,
        headline_n=headline_n,
        nca_min=fmt_int(nca_min),
        nca_median=fmt_int(nca_med),
        nca_max=fmt_int(nca_max),
        nco_min=fmt_int(nco_min),
        nco_median=fmt_int(nco_med),
        nco_max=fmt_int(nco_max),
        n_min="NA",
        n_median="NA",
        n_max="NA",
        n_eff=eff_n,
        total_snps=str(stats["_n_lines"]),
        file_format="PGC sumstats VCF v1.0 (CHROM POS ID REF ALT BETA SE PVAL NGT IMPINFO NEFFDIV2 NCAS NCON)",
        notes=notes,
    )


def process_ptsd_vcf(label, fname, citation, ancestry, headline_n, notes):
    """PTSD 2024 PGC VCF v1.0 - NEFF only (no Nca/Nco)."""
    print(f"[{label}] starting", file=sys.stderr, flush=True)
    full = os.path.join(OUT_DIR, fname)
    stats = stream_columns(
        full,
        col_indices_dict={"NEFF": "NEFF"},
        skip_vcf_meta=True,
        header_must_have=["NEFF"],
    )
    n_stats = stats["NEFF"]
    n_min = int(round(n_stats["min"])) if n_stats else None
    n_med = int(round(n_stats["median"])) if n_stats else None
    n_max = int(round(n_stats["max"])) if n_stats else None
    add_row(
        dataset=label,
        role="secondary outcome",
        path=f"{REL_ROOT}outcomes/{fname}",
        citation=citation,
        ancestry=ancestry,
        headline_n=headline_n,
        nca_min="N/A — VCF release reports NEFF only; per-variant Nca not stored",
        nca_median="N/A",
        nca_max="N/A",
        nco_min="N/A",
        nco_median="N/A",
        nco_max="N/A",
        n_min=fmt_int(n_min),
        n_median=fmt_int(n_med),
        n_max=fmt_int(n_max),
        n_eff=fmt_int(n_max),  # NEFF already reported
        total_snps=str(stats["_n_lines"]),
        file_format="PGC sumstats VCF v1.0 Z-score (CHROM ID POS A1 A2 FREQ NEFF Z P DIRE)",
        notes=notes,
    )


def process_ts_no_n(label, fname, citation, ancestry, headline_n, notes):
    """TS 2019 — SNP CHR BP A1 A2 INFO OR SE P, no N at all."""
    print(f"[{label}] starting", file=sys.stderr, flush=True)
    full = os.path.join(OUT_DIR, fname)
    stats = stream_columns(
        full,
        col_indices_dict={"OR": "OR"},
        header_must_have=["OR", "SE", "P"],
    )
    add_row(
        dataset=label,
        role="secondary outcome",
        path=f"{REL_ROOT}outcomes/{fname}",
        citation=citation,
        ancestry=ancestry,
        headline_n=headline_n,
        nca_min="N/A — Yu 2019 release omits per-variant Nca",
        nca_median="N/A",
        nca_max="N/A",
        nco_min="N/A — Yu 2019 release omits per-variant Nco",
        nco_median="N/A",
        nco_max="N/A",
        n_min="N/A — no per-variant N column",
        n_median="N/A",
        n_max="N/A",
        n_eff="N/A — author to supply from source paper",
        total_snps=str(stats["_n_lines"]),
        file_format="PGC daner-style (SNP CHR BP A1 A2 INFO OR SE P; no Nca/Nco/N)",
        notes=notes,
    )


# =====================================================
# Drive everything
# =====================================================

def main():
    # --- EXPOSURES (4) ---
    process_exposure(
        "Longchamps_2022",
        "Longchamps2022.tsv",
        "Longchamps RJ et al., Hum Genet 2022; PMID 34859289 [ref 30 in v10]",
        "EUR",
        "N = 465,809 (CHARGE + UKB qPCR meta)",
        "INFO >= 0.8 filter applied during harmonization; autosomes only; MHC excluded.",
    )
    process_exposure(
        "Chong_2022",
        "Chong2022.tsv",
        "Chong M et al., eLife 2022; PMID 35023831 [ref 31 in v10]",
        "EUR",
        "N = 383,476 (UKB AutoMitoC SNP-array pipeline)",
        "AutoMitoC pipeline; UKB European subset.",
    )
    process_exposure(
        "Gupta_2023_raw",
        "Gupta2023_raw.tsv",
        "Gupta R et al., Nat Genet 2023; PMID 37563329 [ref 32 in v10]",
        "EUR",
        "N = 155,998 (UKB WGS, unadjusted)",
        "WGS read-depth based; chr_bp_a1_a2 -> rsID lookup from 1KG Phase 3 EUR (9.09M variants mapped).",
    )
    process_exposure(
        "Gupta_2023_adjusted",
        "Gupta2023_adjusted.tsv",
        "Gupta R et al., Nat Genet 2023; PMID 37563329 [ref 32 in v10]",
        "EUR",
        "N = 155,998 (UKB WGS, blood-cell-composition adjusted)",
        "WGS, residualized for nucleated blood cell counts; same rsID lookup.",
    )

    # --- PRIMARY OUTCOMES (5) ---
    process_cc_pgc_daner(
        "ADHD_Demontis_2023",
        "ADHD_Demontis2023.meta.gz",
        "Demontis D et al., Nat Genet 2023; PMID 36702997 [ref 33 in v10]",
        "EUR",
        "FRQ_A 38,691 / FRQ_U 186,843",
        "primary outcome",
        "PGC ADHD 2023 meta release.",
    )
    process_grove_no_n(
        "ASD_Grove_2019",
        "ASD_Grove2019.gz",
        "Grove J et al., Nat Genet 2019; PMID 30804558 [ref 34 in v10]",
        "EUR",
        "18,381 cases / 27,969 controls (headline; per-variant not in file)",
        "primary outcome",
        "Grove 2019 raw release — per-variant Nca/Nco not stored; headline N taken from manuscript v10 Methods.",
    )
    process_cc_pgc_daner(
        "BD_Mullins2024_EUR",
        "BD_Mullins2024_EUR.gz",
        "Mullins N et al., PGC 2024 (BD) EUR no-23andMe public release [ref 35 in v10]",
        "EUR",
        "max per-variant 59,287 / 781,022",
        "primary outcome",
        "no-23andMe public release; SNP CHR BP A1 A2 ... format (non-canonical column order).",
    )
    process_cc_pgc_daner(
        "MDD_Adams2025_EUR",
        "MDD_Adams2025_EUR.gz",
        "Adams MJ et al., Cell 2025 (MDD Adams2025) EUR no-23andMe release [ref 36 in v10]",
        "EUR",
        "FRQ_A 412,305 / FRQ_U 1,588,397",
        "primary outcome",
        "no-23andMe public release; also includes Neff column.",
    )
    process_cc_pgc_daner(
        "SCZ_PGC3_UKBdedup",
        "SCZ_PGC3_UKBdedup.gz",
        "Trubetskoy V et al., Nature 2022 (PGC3-SCZ); PMID 35396580 [ref 37 in v10]",
        "EUR",
        "FRQ_A 67,323 / FRQ_U 93,456",
        "primary outcome",
        "PGC3 European UK-Biobank-deduplicated public release.",
    )

    # --- SECONDARY OUTCOMES (9 actually: ANX OCD PTSD AN BPD TS PPD SUD + ADR x2) ---
    process_cc_pgc_daner(
        "ANX_2026",
        "ANX_2026.gz",
        "Friligkou E et al., 2026 (Anxiety GWAS) [ref 21 in v10]",
        "EUR",
        "FRQ_A 122,083 / FRQ_U 729,602",
        "secondary outcome",
        "PGC anxiety 2026 release.",
    )
    process_cc_pgc_daner(
        "OCD_2025",
        "OCD_2025.gz",
        "Strom NI et al., 2025 (PGC OCD) [ref 22 in v10]",
        "EUR",
        "FRQ_A 23,493 / FRQ_U 1,114,613",
        "secondary outcome",
        "PGC OCD 2025 release.",
    )
    process_ptsd_vcf(
        "PTSD_Maihofer_2024_EUR",
        "PTSD_2024_EUR.vcf.gz",
        "Maihofer AX et al., Nat Genet 2024 (PGC-PTSD F3); PMID 38637617; DOI 10.1038/s41588-024-01707-9 [ref 23 in v10]",
        "EUR",
        "137,136 cases / 1,085,746 controls (EUR subset; VCF reports NEFF only)",
        "Z-score VCF format; Z->beta/SE transform applied during MR; per-variant Nca/Nco not in file.",
    )
    process_an_vcf(
        "AN_Watson_2019",
        "AN_2019.gz",
        "Watson HJ et al., Nat Genet 2019 (PGC-ED AN2) [ref 24 in v10]",
        "EUR",
        "16,992 cases / 55,525 controls (headline)",
        "PGC sumstats VCF v1.0; per-variant NCAS/NCON present.",
    )
    process_cc_pgc_daner(
        "BPD_2025",
        "BPD_2025.gz",
        "PGC-BPD 2025 release [ref 25 in v10]",
        "EUR",
        "FRQ_A 12,339 / FRQ_U 1,041,717",
        "secondary outcome",
        "PGC borderline personality disorder 2025 release (primary, includes UKB).",
    )
    process_ts_no_n(
        "TS_Yu_2019",
        "TS_2019.gz",
        "Yu D et al., Am J Psychiatry 2019 (PGC Tourette syndrome) [ref 26 in v10]",
        "EUR",
        "4,819 cases / 9,488 controls (headline)",
        "Yu 2019 release omits per-variant Nca/Nco.",
    )
    process_ppd_2023(
        "PPD_2023_EUR",
        "PPD_2023_EUR.tsv.gz",
        "Guintivano J et al., 2023 (postpartum depression PGC) [ref 27 in v10]",
        "EUR",
        "TotalSampleSize max ~54,475 (headline; case/control split in source paper)",
        "METAL output; per-variant TotalSampleSize present; per-variant Nca/Nco not in file.",
    )
    process_sud_2023(
        "SUD_Hatoum_2023_EUR",
        "SUD_2023_Hatoum_EUR.txt.gz",
        "Hatoum AS et al., Nat Ment Health 2023 (general addiction factor) [ref 28 in v10]",
        "EUR",
        "N approx 1,025,000 (headline; per source paper)",
        "SE column absent in release -> excluded from forward MR; covered indirectly via CDG3-F5.",
    )
    process_cc_pgc_daner(
        "ADR_PercImprov_2021",
        "ADR_PercImprovEUR_2021.gz",
        "Pain O et al., 2021 (antidepressant response percentage improvement) [ref 44 in v10]",
        "EUR",
        "N = 5,376 (FRQ_A 5,376 / FRQ_U 5,376; quantitative residual)",
        "sensitivity outcome",
        "Antidepressant percentage-improvement; quantitative trait modelled in daner-style format.",
    )
    process_cc_pgc_daner(
        "ADR_NonRemission_2021",
        "ADR_NonRemissionEUR_2021.gz",
        "Pain O et al., 2021 (antidepressant non-remission) [ref 44 in v10]",
        "EUR",
        "FRQ_A 3,303 / FRQ_U 1,852",
        "sensitivity outcome",
        "Antidepressant non-remission endpoint.",
    )

    # --- SENSITIVITY OUTCOMES (3 UKBB-dedup pairs) ---
    process_cc_pgc_daner(
        "BD_2021_noUKBB",
        "BD_2021_noUKBB.gz",
        "Mullins N et al., Nat Genet 2021 (PGC2-BD non-UKBB) [ref 59 in v10]",
        "EUR",
        "FRQ_A 40,463 / FRQ_U 313,436",
        "sensitivity outcome",
        "PGC2-BD without UKBB; for UKB-overlap sensitivity replication of BD primary.",
    )
    process_cc_pgc_daner(
        "MDD_2018_noUKBB",
        "MDD_2018_noUKBB.gz",
        "Wray NR et al., Nat Genet 2018 (PGC2-MDD non-UKBB) [ref 60 in v10]",
        "EUR",
        "FRQ_A 45,396 / FRQ_U 97,250",
        "sensitivity outcome",
        "PGC2-MDD without UKBB; for UKB-overlap sensitivity replication of MDD primary.",
    )
    process_cc_pgc_daner(
        "BPD_2025_noUKBB",
        "BPD_2025_noUKBB.gz",
        "PGC-BPD 2025 release (sub-cohort excluding UKB) [ref 25 in v10]",
        "EUR",
        "FRQ_A 12,157 / FRQ_U 1,039,897",
        "sensitivity outcome",
        "BPD2025 sub-cohort with UKB-overlapping samples removed.",
    )

    # --- CDG3 FACTORS (6) ---
    process_cdg3(
        "CDG3_pFactor",
        "CDG3_PFactor.tsv.gz",
        "Grotzinger AD et al., Nat Hum Behav 2022 (CDG3 14-disorder genomic-SEM) [ref 43 in v10]",
        "EUR",
        "N ~ 2,168,621 (general psychopathology pFactor; headline from v10)",
        "Genomic-SEM factor; no per-variant N column in release; effective N is composite from genomic-SEM weights.",
    )
    process_cdg3(
        "CDG3_F1_Compulsive",
        "CDG3_F1_Compulsive.tsv.gz",
        "Grotzinger AD et al., Nat Hum Behav 2022 (CDG3 F1 Compulsive) [ref 43 in v10]",
        "EUR",
        "N ~ 188,000 (F1 compulsive factor; headline from v10)",
        "Genomic-SEM factor; F1.",
    )
    process_cdg3(
        "CDG3_F2_SchizoBipolar",
        "CDG3_F2_SchizoBipolar.tsv.gz",
        "Grotzinger AD et al., Nat Hum Behav 2022 (CDG3 F2 Schizo-Bipolar) [ref 43 in v10]",
        "EUR",
        "N ~ 600,000 (F2 schizo-bipolar factor; headline from v10)",
        "Genomic-SEM factor; F2.",
    )
    process_cdg3(
        "CDG3_F3_Neurodev",
        "CDG3_F3_Neurodev.tsv.gz",
        "Grotzinger AD et al., Nat Hum Behav 2022 (CDG3 F3 Neurodevelopmental) [ref 43 in v10]",
        "EUR",
        "N = 84,760 (F3 neurodevelopmental factor; headline from v10)",
        "Genomic-SEM factor; F3.",
    )
    process_cdg3(
        "CDG3_F4_Internalizing",
        "CDG3_F4_Internalizing.tsv.gz",
        "Grotzinger AD et al., Nat Hum Behav 2022 (CDG3 F4 Internalizing) [ref 43 in v10]",
        "EUR",
        "N = 1,637,337 (F4 internalizing factor; headline from v10)",
        "Genomic-SEM factor; F4.",
    )
    process_cdg3(
        "CDG3_F5_SubstanceUse",
        "CDG3_F5_SubstanceUse.tsv.gz",
        "Grotzinger AD et al., Nat Hum Behav 2022 (CDG3 F5 Substance/externalizing) [ref 43 in v10]",
        "EUR",
        "N ~ 920,000 (F5 substance-use factor; headline from v10)",
        "Genomic-SEM factor; F5.",
    )

    # --- NEGATIVE CONTROLS (2) ---
    process_neg_ctrl(
        "NegCtrl_BMI_Yengo_2022",
        "NegCtrl_BMI.gz",
        "Yengo L et al., Nature 2022 (GIANT BMI) [ref 45 in v10]",
        "EUR",
        "max per-variant N ~ 598,895",
        "Negative control — anthropometric, no a-priori mtDNA-CN causal pathway.",
    )
    process_neg_ctrl(
        "NegCtrl_Height_Yengo_2022",
        "NegCtrl_height.gz",
        "Yengo L et al., Nature 2022 (GIANT standing height) [ref 45 in v10]",
        "EUR",
        "max per-variant N ~ 605,309",
        "Negative control — anthropometric, no a-priori mtDNA-CN causal pathway.",
    )

    # --- Emit TSV (tab-separated; pipe-delimited inside cell free text where commas may exist) ---
    columns = [
        "dataset", "role", "path", "citation", "ancestry", "headline_n",
        "nca_min", "nca_median", "nca_max",
        "nco_min", "nco_median", "nco_max",
        "n_min", "n_median", "n_max",
        "n_eff", "total_snps", "file_format", "notes",
    ]
    os.makedirs(os.path.dirname(OUTPUT_TSV), exist_ok=True)
    with open(OUTPUT_TSV, "w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in ROWS:
            line = "\t".join(str(row.get(c, "NA")) for c in columns)
            fh.write(line + "\n")
    print(f"WROTE {OUTPUT_TSV}", file=sys.stderr)
    print(f"ROWS: {len(ROWS)}", file=sys.stderr)


if __name__ == "__main__":
    main()
