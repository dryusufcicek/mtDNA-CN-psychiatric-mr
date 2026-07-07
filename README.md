# Cross-platform Mendelian randomization of mitochondrial DNA copy number across psychiatric disorders

Analysis code and derived results for a two-sample Mendelian randomization (MR) study evaluating
genetically proxied blood **mitochondrial DNA copy number (mtDNA-CN)** across psychiatric phenotypes,
using **four partially overlapping mtDNA-CN GWAS exposure panels** (qPCR, SNP-array, and WGS-derived,
with and without blood-cell-composition adjustment).

Yusuf Cicek et al.

## Summary of findings
- Higher genetically proxied mtDNA-CN was associated most robustly with **lower autism spectrum
  disorder (ASD) liability** (β = −0.182; design-effect-corrected *p* = 5.97 × 10⁻³; BH-FDR *q* = 0.030).
- The ASD–ADHD **neurodevelopmental factor (F3)** showed a directionally concordant nominal signal;
  ADHD alone did not survive correction.
- No robust effects for schizophrenia, bipolar disorder, major depression, general psychopathology, or
  other disorder-class factors; reverse-direction MR was null.
- Because the four exposure panels share UK Biobank participants, panel estimates are correlated
  (mean inter-panel *r* = 0.70) and were pooled with a **design-effect (effective-sample-size)
  correction** rather than as independent replications.
- Gene-set enrichment (MAGMA × MitoCarta3.0) and exploratory brain eQTL colocalization (GTEx v10) are
  hypothesis-generating; three-way (eQTL × mtDNA-CN × disorder) colocalization and SuSiE conditional
  fine-mapping did not establish a mediator.

## Repository structure
```
scripts/                 # full analysis pipeline (R / Python / bash)
results/                 # small derived summary results (TSV)
  coloc/                 #   coloc.abf, three-way, and SuSiE outputs
  forward_mr/  meta_mr/  reverse_mr/  mrlap/  ldsc_rho/  …
supplementary_tables/    # Supplementary_Tables.xlsx (S1–S18) + per-table TSVs
docs/ENVIRONMENT.md      # compute environment + R/Python package versions
```

## Pipeline overview (`scripts/`)
| Step | Script(s) | Purpose |
|---|---|---|
| 0 | `00_install_packages.R` | install R dependencies |
| 1 | `01_inspect_headers.R`, `01b_data_qc.R`, `01c_build_rsid_lookup.R` | QC; build GRCh37 chr:pos→rsID lookup for the WGS panels |
| 2 | `02_harmonize_exposures.R` | harmonize the four mtDNA-CN exposure panels |
| 3 | `03_forward_mr.R`, `add_mr_raps.R`, `lib_robust.R` | forward two-sample MR (IVW, MR-Egger, weighted median, MR-RAPS, MR-PRESSO) |
| 4 | `04_meta_mr.R`, `corr_aware_meta.R` | design-effect-corrected cross-panel synthesis + naïve REML sensitivity |
| 5 | `C1_ldsc_rho.R` | bivariate-LDSC sample-overlap ρ (GenomicSEM) |
| 6 | `06_mrlap.R` | MRlap sample-overlap-corrected sensitivity |
| 7 | `07_reverse_mr.R` | reverse-direction MR |
| 8 | `08a_extract_cis_eqtl.py`, `08b_liftover_eqtl.py`, `08b_coloc_fixed.R`, `08c_threeway_coloc.R`, `08d_coloc_susie_v2.R` | GTEx v10 brain eQTL colocalization: eQTL×disorder, three-way, and SuSiE fine-mapping |
| 9 | `09a_prep_magma_inputs.py`, `09_magma_mitocarta.sh` | MAGMA × MitoCarta3.0 gene-set / sub-pathway enrichment |
| 10 | `10_rare_variant_sensitivity.py` | rare / low-frequency-instrument (MAF) leave-out sensitivity of the forward MR (Supp. Table 17) |
| 11 | `11_instrument_vep_annotation.py` | Ensembl (GRCh38.p14) functional annotation of the instruments (Supp. Table 18) |
| – | `build_supp_table_1.py`, `rebuild_supp_t2.R` | supplementary table assembly |

`run_mr_pipeline.sh` orchestrates the core MR steps; `submit_*.sh` are SLURM wrappers for the HPC
environment. Paths are written for the TRUBA cluster (`/arf/...`) and must be adjusted for other systems.

## Data availability (raw GWAS summary statistics are **not** redistributed here)
The exposure and outcome summary statistics are governed by their original licenses/consortium terms
and must be obtained from the sources below. The pipeline expects them under `data/` (git-ignored).

**Exposures (mtDNA-CN):**
| Panel | Source | Accession / URL |
|---|---|---|
| Longchamps 2022 (qPCR meta) | Arking lab | `arkinglab.org/upload/mtDNA_CN_GWAS/mtDNA_CN.ALLm2.bgen.stats.gz` |
| Chong 2022 (SNP-array, AutoMitoC) | GWAS Catalog | GCST90026372 |
| Gupta 2023 (WGS, blood-cell-adjusted) | GWAS Catalog | GCST90268497 |
| Gupta 2023 (WGS, unadjusted) | GWAS Catalog | GCST90268498 |

**Outcomes:** PGC psychiatric GWAS (ADHD Demontis 2023; ASD Grove 2019; BD O'Connell 2025;
MDD Adams 2025; SCZ Trubetskoy 2022; anxiety, OCD, PTSD, anorexia, BPD, Tourette, postpartum depression,
antidepressant response) via the Psychiatric Genomics Consortium portal; cross-disorder factor GWAS
(Grotzinger 2025) via the associated figshare; height/BMI (Yengo 2018, GIANT). Brain cis-eQTL:
**GTEx v10** (GTEx Portal). Mitochondrial gene set: MitoCarta3.0.

Derived **summary results** that support the figures and tables are included under `results/` and
`supplementary_tables/`; they remain subject to the upstream data licenses.

## Environment
See `docs/ENVIRONMENT.md`. Core: R 4.4 (data.table, TwoSampleMR, MR-PRESSO, MRlap, MendelianRandomization,
metafor, GenomicSEM, coloc, susieR), Python 3.11 (pyarrow, pyliftover), MAGMA, PLINK v1.9.

## Citation
Cicek Y, Çelik MC. *A cross-platform Mendelian randomization study of mitochondrial DNA copy number
across psychiatric disorders.* (Manuscript; citation/DOI to be added on publication.)

## License
Code is released under the MIT License (see `LICENSE`). Derived results are provided for reproducibility
and remain subject to the licenses of the source GWAS datasets; raw summary statistics are not
redistributed in this repository.
