# Replication results vs. pre-registered predictions

Independent stack (pandas / scipy / statsmodels), fresh code, data re-read from the article's
Supplementary Table 1, strata re-derived from the source file. Full log: replication_output.txt

| # | Prediction | Result | Verdict |
|---|---|---|---|
| Q1.3 | Effect size does not track evidence strength | Spearman rho = −0.013, P = 0.93 | **CONFIRMED** |
| Q1.4 | Randomized/published arms show no larger effects | randomized −0.064 vs −0.038, P = 0.39; published −0.047 vs −0.055, P = 0.62; Spearman(n, effect) rho = −0.04, P = 0.80 | **CONFIRMED** |
| Q2.1 | 19 / 13 / 26 / 3, and 6 nominal increases | exact | **CONFIRMED** |
| Q2.2 | 6 disease, 1 comparator, copaiba, Trulacta, ketamine, 3 diets | exact | **CONFIRMED** |
| Q2.3 | Exercise and CR in no-effect group; exercise 0/108 | exact | **CONFIRMED** |
| Q2.4 | Qualitative ranking survives 5 metric variants | survives all 5 | **CONFIRMED** |
| Q2.5 | Copaiba robust | passes 4 of 5 variants (fails Wilcoxon only) | **CONFIRMED** |
| Q3.1 | Reproduce article stratum values | all 5 exact to 3 decimals | **CONFIRMED** |
| Q3.2 | No clock passes in healthy-41 or healthy-34 | none in either | **CONFIRMED** |
| Q3.3 | Disease result survives leave-one-out | minimum 3 clocks pass across all 10 drops | **CONFIRMED** |
| Q3.4 | Not confounded by n | disease median n = 25 vs healthy 19, P = 0.60 | **CONFIRMED** |
| Q3.5 | Survives nonparametric test | Mann-Whitney P = 0.0013–0.0017 for PCGrimAge, SystemsAge, PCPhenoAge (DunedinPACE P = 0.32) | **CONFIRMED** |
| Q3.6 | Disease significant adjusting for n and duration | OLS beta = −0.120, P = 0.0041 | **CONFIRMED** |
| Q4.1 | SE reconstruction exact | max round-trip error 5.1e-15 over 816 cells | **CONFIRMED** |
| Q4.2 | Pooled contrasts null | all 5 null under all 3 methods (P = 0.11–0.98) | **CONFIRMED** |
| Q4.6 | Exercise contrasts null under all methods | smallest P = 0.092 / 0.093 / 0.052 | **CONFIRMED** |
| **Q4.3** | **CIs exclude published estimate for 2 of 5** | fixed-effect 2/5, random-effects 1/5, **Hartung-Knapp 0/5** | **FAILED (fragile)** |
| **Q4.5** | **MDE roughly 0.06–0.09, i.e. informative** | **MDE 0.059–0.128, exceeds the published effect for all 5 clocks** | **FAILED (underpowered)** |

## Action taken, per the pre-registered decision rule

The rule stated: *"If Q4.5 shows the MDE far exceeds published effects, argument 4 is downgraded to
'no evidence of effect, and underpowered' and the CI-exclusion claim is dropped."*

The CI-exclusion claim was dropped. The letter makes no claim that the confidence intervals exclude
the article's published estimates, and the checks that tested it have been removed from
`verify_all.py`. What the letter retains is the descriptive result, that the treated-minus-control
difference is null for the recommended clocks in the three trials that permit it to be computed.
That result is unimpeachable; the inferential claim held only under fixed-effect z-based intervals.

## Net assessment

Arguments 1, 2 and 3 stand, and two are strengthened:

- **Argument 1** gains its first positive empirical support. Effect sizes are unrelated to whether an
  arm was randomized, published, or large; the evidence-tier association is flat (rho = −0.01,
  P = 0.93). Whatever the analysis measures, it does not track evidence that an intervention works.
- **Argument 3** is more robust than claimed: it survives leave-one-out (including dropping ART,
  n = 183), is not confounded by sample size, holds nonparametrically, and holds in regression
  adjusting for n and duration.
- **Argument 2** survives all five metric variants; copaiba survives four.
- **Argument 4**'s core claim (the randomized contrasts are null) holds under every pooling method,
  but the letter can no longer claim the intervals exclude the published estimates.
