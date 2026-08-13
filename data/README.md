# data/

- `total_c9_normalized.tsv` — one row per donor (8 control, 8 C9-ALS). Columns:
  `sample` (GEO GSM accession), `group` (`control`/`c9ALS`), `norm_count`
  (DESeq2 median-of-ratios normalized count for human C9ORF72, ENSG00000147894).
  This is the Figure 1b input.
- `C9ORF72_control_vs_C9ALS.xlsx` — the same 16 donors plus a group summary whose
  mean/median/SD/SEM/fold are live formulas over the donor rows. The Mann-Whitney
  p-value is computed in analysis/c9_stats.py -> results/stats.tsv, not duplicated here.

Normalized counts were produced upstream (STAR/GRCh38, GENCODE v45; DESeq2 across
the full 26-sample cohort). Raw reads: GEO GSE67196 (via SRA).
