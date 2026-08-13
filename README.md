# C9ORF72 expression in C9-ALS cerebellum — Figure 1b source

Data and code to regenerate **Figure 1b** and its statistic: total *C9ORF72*
normalized RNA abundance in human post-mortem cerebellum, control vs C9-ALS
(reanalysis of GEO **GSE67196**, Prudencio et al. 2015).

## Provenance and scope

Normalized per-donor counts were produced upstream from GSE67196 (STAR/GRCh38,
GENCODE v45 alignment, reverse-stranded; DESeq2 median-of-ratios normalization
computed across the full 26-sample GSE67196 cerebellum cohort). **This repository
takes those normalized counts as input and reproduces Figure 1b (the control and
C9-ALS subset, 16 donors) and its statistics; it does not regenerate counts from
raw FASTQ.** Figure 1b establishes lower bulk-tissue
normalized C9ORF72 RNA abundance in this cohort.

## Regenerate

```bash
python analysis/c9_stats.py         # -> results/stats.tsv   (two-sided Mann-Whitney U, exact, SciPy)
Rscript analysis/plot_c9_panel.R    # -> Fig1b_C9ORF72.{pdf,svg,png}  (reads the p-value from results/)
```

The figure reads the p-value from `results/stats.tsv`; the test is computed once,
in the Python step, and not recomputed in R.

## Data dictionary — `data/total_c9_normalized.tsv`

| column | definition |
|---|---|
| `sample` | GEO sample accession (GSM) |
| `group` | `control` or `c9ALS` |
| `norm_count` | DESeq2 median-of-ratios normalized count for human C9ORF72 (ENSG00000147894) |

`data/C9ORF72_control_vs_C9ALS.xlsx` holds the same 16 donors plus a group
summary whose mean/median/SD/SEM/fold are live formulas over the donor rows.

## Requirements
Statistics: Python 3.10.12 with NumPy 2.2.6 and SciPy 1.15.3 (pip install -r requirements.txt). Figure: R >= 4.2 with ggplot2.
## License

Code: MIT (`LICENSE`). Source data: GEO GSE67196.
