# Pre-registration: independent replication of the four arguments
Written BEFORE any analysis in this folder was run. Predictions are stated so they can fail.

Independence measures: fresh code, pandas/scipy/statsmodels (the original used hand-rolled
t-distribution code and openpyxl loops); data re-read from the article's .xlsx; no import of
any earlier script; strata re-derived from the source file rather than from hard-coded lists.

---

## ARGUMENT 1 - The central premise cannot be met

**Q1.1** Does the inclusion criterion select on evidence of efficacy?
*Predicted:* No. Selection is on hypothesis plus data availability. Evidence: Methods text, and the
presence of arms with no efficacy literature at all (unpublished commercial supplements).

**Q1.2** Can the analysis fail? Is any outcome inconsistent with "clocks are responsive"?
*Predicted:* No. A pooled decrease is read as responsiveness; a null is attributable to the
intervention. The design has no falsifying outcome.

**Q1.3 (testable)** If the interventions genuinely targeted aging, effect size should track external
evidence strength. Assign each arm an evidence tier from its own metadata (randomized + published +
whole blood + n>=30 = strongest; unpublished single-arm = weakest) and test association with effect.
*Predicted:* No monotonic association; possibly inverse (weaker designs, larger effects).

**Q1.4 (testable)** Do effect sizes track study quality proxies (randomization, publication, n)?
*Predicted:* Randomized and published arms show effects no larger than non-randomized/unpublished.

---

## ARGUMENT 2 - The ranking, on their own metric

**Q2.1** Does the Fig 3d metric reproduce the published counts (19 / 13 / 26 nominal-dec / corrected-dec /
no-effect, and 3 corrected increases)?
*Predicted:* Yes exactly, except 6 rather than 5 nominal increases.

**Q2.2** Composition of the 13: how many treat active disease, how many are comparator arms, which
supplements appear?
*Predicted:* 6 disease, 1 comparator arm, copaiba oil, Trulacta human-milk, ketamine, 3 diets.

**Q2.3** Where do exercise and caloric restriction rank?
*Predicted:* Both in the no-effect group; exercise 0/108 biomarkers at nominal P.

**Q2.4 (robustness)** Does the ranking survive metric variation: median instead of mean; trimmed mean;
Wilcoxon signed-rank instead of t; Gen2+ clocks only; all 108 biomarkers instead of 16?
*Predicted:* The qualitative claim (disease arms and supplements above exercise/CR) survives all five.
Individual membership may shift by 1-3 arms.

**Q2.5 (falsification test)** Is copaiba's rank an artifact of the 16-clock choice?
*Predicted:* Copaiba remains a significant decrease under most variants; if it fails under >=3 of 5,
the letter's headline example is fragile and must be re-worded.

---

## ARGUMENT 3 - Decreases concentrate in disease, absent in healthy

**Q3.1** Reproduce the article's own stratum values (disease PCPhenoAge -0.320 P=0.0057 etc.).
*Predicted:* Exact reproduction.

**Q3.2** In 41 healthy arms, does any clock pass P<0.00833? In 34 (excluding reference arms)?
*Predicted:* None in either.

**Q3.3 (leave-one-out)** Is the disease result driven by ART (n=183)?
*Predicted:* At least 2 clocks still pass in the disease stratum after dropping any single arm.

**Q3.4 (confounding)** Is the disease/healthy difference explained by sample size or duration?
*Predicted:* No. Disease arms are not systematically larger; effect persists adjusting for n.

**Q3.5 (nonparametric)** Does the disease/healthy difference survive Mann-Whitney?
*Predicted:* Yes for the recommended clocks.

**Q3.6 (model-based)** In a regression of arm-level effect on disease status with n and duration as
covariates, is disease status significant?
*Predicted:* Yes, negative coefficient.

---

## ARGUMENT 4 - Randomized contrasts are null

**Q4.1** Is SE = |d| / t(P, n-1) exact? Round-trip every cell.
*Predicted:* Recovered P equals reported P to <1e-9 for all cells.

**Q4.2** Pooled 3-contrast estimates for the 5 clocks.
*Predicted:* All null (P>0.05); PCGrimAge ~+0.001, SystemsAge ~+0.008.

**Q4.3** Do the CIs exclude the article's published pooled estimates?
*Predicted:* For exactly 2 of 5 (PCGrimAge, SystemsAge).

**Q4.4 (method robustness)** Fixed-effect vs DerSimonian-Laird random-effects vs Hartung-Knapp;
t-based rather than z-based intervals.
*Predicted:* All null. Exclusion of the published estimate may weaken under Hartung-Knapp with 3 studies.

**Q4.5 (power)** What is the minimum detectable difference with 3 contrasts at 80% power?
*Predicted:* Roughly 0.06-0.09 s.d. -- i.e. comparable to the published effects, which is why the
comparison is informative rather than merely underpowered. IF the MDE is much larger than the
published effects, argument 4 is weaker than stated and must be softened.

**Q4.6 (falsification)** Do the exercise factorial contrasts remain null under all pooling methods?
*Predicted:* Yes.

---

## Decision rules, fixed in advance

- Any prediction that fails is reported in the summary, whether or not it helps the letter.
- If Q2.5 fails (copaiba fragile), the abstract's headline example changes.
- If Q4.5 shows the MDE far exceeds published effects, argument 4 is downgraded to "no evidence of
  effect, and underpowered" and the CI-exclusion claim is dropped.
- If Q3.3 shows the disease result depends on one arm, argument 3 is softened.
