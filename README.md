# Reanalysis of Sehgal et al., *Nature Medicine* 2026

Code and data for a Matters Arising on:

> Sehgal, R. *et al.* Responsiveness of epigenetic aging biomarkers to longevity interventions in humans. *Nature Medicine* (2026). https://doi.org/10.1038/s41591-026-04562-9

Every number in the letter is recomputed from the article's own Supplementary Table 1, at the
article's own significance thresholds. Nothing here uses data the authors did not publish.

## Getting the data

The article's supplementary file is not redistributed here. Download **Supplementary Table 1**
(MOESM4) from the paper and save it in this folder as:

```
SupplementaryTable1_Sehgal2026.xlsx
```

## Reproducing everything

```bash
pip install pandas numpy scipy openpyxl statsmodels
python3 verify_all.py          # 41 named checks, expect 0 failures
python3 make_figure_data.py    # regenerates the figure inputs
Rscript make_figure1.R         # Figure 1
Rscript make_extended_data_fig1.R   # Extended Data Fig. 1
```

`verify_all.py` is the entry point. It re-derives the article's published counts, the strata, the
randomized differences and the exercise contrasts, and fails loudly if any of them moves. Every check
corresponds to a claim made in the letter; nothing is verified here that the letter does not assert.

## What each file does

| File | Purpose |
|---|---|
| `verify_all.py` | 41 named checks covering every numerical claim in the letter |
| `league_table.py` | All 51 arms ranked on that metric, with the article's own verdict per arm |
| `randomized_contrasts.py` | Treated-minus-control differences in the trials that permit them |
| `make_figure_data.py` | Regenerates the figure inputs; verifies the curated panel-a file against the source |
| `make_figure1.R` | Figure 1 (panels a, b, c) |
| `make_extended_data_fig1.R` | Extended Data Fig. 1 (exercise contrast, all 16 clocks) |
| `HOW_TO_VERIFY.md` | Checks the headline claims by hand in Excel, without running any code |
| `replication/` | Independent pre-registered replication (see below) |

`fig_data_a.csv` is a curated file: the display labels and category flags are editorial, but its
effect sizes, sample sizes and verdicts are checked against the source table by
`make_figure_data.py`, which reports any mismatch.

## Independent replication

`replication/` contains a pre-registered replication of the four arguments, written on a different
software stack (pandas / scipy / statsmodels rather than the hand-rolled t-distribution code used
elsewhere here).

- `PREREGISTRATION.md`: questions and predictions, written **before** the analysis was run
- `rep_analysis.py`: the analysis
- `FINDINGS.md`: results against the predictions

Sixteen of eighteen predictions were confirmed and **two failed**. Both failures concerned the
randomized-difference argument: with three trials the minimum detectable difference exceeds the
published effect sizes, and the exclusion of the published estimates from the confidence intervals
held only under fixed-effect intervals, not under Hartung-Knapp. The letter was changed
accordingly. The failures are documented rather than removed.

## Methods in brief

Arm-level metric: one-sample t-test (df = 15) across the 16 prominent-clock effect sizes per arm,
reproducing Fig. 3d of the article; threshold *P* < 0.00625 (Bonferroni at *M*ₑ = 8, from the
article's Supplementary Table 20). Standard errors are reconstructed as
SE = |*d*| / *t*(*P*, *n* − 1) from the reported effect size, *P* value and sample size; the
difference between two arms has variance SE₁² + SE₂².

## Licence

Code is released under the MIT Licence. The article's supplementary data remain the property of
their authors and publisher and are not included here.
