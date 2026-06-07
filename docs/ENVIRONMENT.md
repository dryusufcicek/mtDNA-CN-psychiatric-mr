# TRUBA Setup — iter_003

**Cluster:** TRUBA (arf-ui1)
**User:** ycicek
**Project root:** `/arf/home/ycicek/scz-research/iter_003_mtdna_psy_mr/`
**Setup date:** 2026-05-22

## Environment

### micromamba (~/bin/micromamba)
- Version: 2.6.2
- Pre-existing — Yusuf installed previously
- Used instead of full conda/miniforge (single binary, ~18 MB)

### R environment: `iter003_mr`
Created via:
```bash
~/bin/micromamba create -y -n iter003_mr -c conda-forge \
  r-base=4.4 r-data.table r-dplyr r-ggplot2 r-remotes r-devtools \
  r-meta r-metafor r-mendelianrandomization
```

### Activate
```bash
export PATH=$HOME/bin:$PATH
eval "$(micromamba shell hook --shell bash)"
micromamba activate iter003_mr
```

### GitHub R packages (installed manually after env creation)
```r
remotes::install_github("MRCIEU/TwoSampleMR", upgrade = "never")
remotes::install_github("rondolab/MR-PRESSO", upgrade = "never")
remotes::install_github("WSpiller/RadialMR", upgrade = "never")
remotes::install_github("MRCIEU/ieugwasr", upgrade = "never")
```

## Filesystem layout (TRUBA)

```
~/scz-research/iter_003_mtdna_psy_mr/
├── README.md
├── PROJECT_LOG.md
├── data/
│   ├── exposures/                     # 4 mtDNA-CN GWAS (~2.5 GB) — rsynced from local
│   │   ├── Longchamps2022_mtDNA_CN.ALLm2.bgen.stats.gz   (830 MB)
│   │   ├── Chong2022_mtDNA_CN_EUR_GCST90026372.tsv.gz   (205 MB)
│   │   ├── Gupta2023_mtDNA_CN_adjusted_GCST90268497.tsv.gz   (742 MB)
│   │   └── Gupta2023_mtDNA_CN_raw_GCST90268498.tsv.gz   (743 MB)
│   └── outcomes/                      # SYMLINKS to ~/EVOSCZ/data/raw/gwas/...
│       ├── SCZ_PGC3_EUR.tsv.gz → ~/EVOSCZ/.../pgc3_sumstats/PGC3_EUR_autosome.vcf.tsv.gz
│       ├── ADHD_Demontis2023.meta.gz → ~/EVOSCZ/.../adhd/...
│       ├── ASD_Grove2019.gz → ~/EVOSCZ/.../asd/...
│       ├── BD_Mullins2024_EUR.gz → ~/EVOSCZ/.../bd/...
│       ├── MDD_Adams2025_EUR.gz → ~/EVOSCZ/.../mdd/...
│       ├── CDG3_PFactor.tsv.gz → ~/EVOSCZ/.../cdg3/...
│       ├── CDG3_F3_Neurodev.tsv.gz → ~/EVOSCZ/.../cdg3/...
│       ├── CDG3_F4_Internalizing.tsv.gz → ~/EVOSCZ/.../cdg3/...
│       ├── NegCtrl_BMI.gz → ~/EVOSCZ/.../negative_control/...
│       └── NegCtrl_height.gz → ~/EVOSCZ/.../negative_control/...
├── scripts/                            # Pipeline R scripts (rsynced)
├── results/                            # Output (created by scripts)
├── logs/                               # Run logs
├── manuscript_draft/                   # Brief Report
└── docs/                               # Metadata + preregistration
```

## Running scripts on TRUBA

### Interactive (login node, light tasks):
```bash
cd ~/scz-research/iter_003_mtdna_psy_mr
export PATH=$HOME/bin:$PATH
eval "$(micromamba shell hook --shell bash)"
micromamba activate iter003_mr
Rscript scripts/01_inspect_headers.R
```

### SLURM batch (compute-heavy: MR pipeline):
```bash
sbatch scripts/run_mr_pipeline.sh
```

SLURM partition: `orfoz` (3-day limit)

## Disk

- Home quota: not yet established (Yusuf to verify with `lfs quota`)
- Current TRUBA usage: 151 GB
- This project adds: ~2.5 GB exposures + symlinks (0 GB outcomes — shared with EVOSCZ)

## Backup strategy

- Code → GitHub repo (public on submission)
- Results → rsync to local periodically + TRUBA scratch
- Data exposures → preserved in TRUBA project; original public URLs documented in DATA_ACQUISITION.md

## Local cleanup decision (pending)

Local copy at `/Users/yusuf/scz-research/iter_003_mtdna_psy_mr/data/exposures/` (~2.5 GB) can be **DELETED** after successful TRUBA rsync confirmed and md5sum verified. Yusuf to confirm before deletion (per `no-delete-without-asking` rule).
