# How to check every number yourself

Everything in the letter comes from **one file the authors published**: their Supplementary Table 1.
No private data, no access request, no reprocessing of methylation arrays.

---

## 1. Get their supplementary file

Article page: <https://www.nature.com/articles/s41591-026-04562-9>
Scroll to **Supplementary information** near the bottom. Supplementary Table 1 is the first .xlsx
("TranslAGE datasets and biomarker catalog").

Direct link. This is the file we used, byte-identical to `SupplementaryTable1_Sehgal2026.xlsx`:

    https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41591-026-04562-9/MediaObjects/41591_2026_4562_MOESM4_ESM.xlsx

Also useful:
- **Supplementary Table 20** (MOESM23): their multiple-testing thresholds. This is where P < 0.00625
  (Fig. 3d), P < 0.00833 (Fig. 4b/5b) and P < 0.0125 (Fig. 5c) come from.
- **Supplementary Table 7** (MOESM10): their published pooled per-clock effects.
- **Peer review file** (MOESM2): the three referee reports and both author rebuttals.

Numbering rule: Supplementary Table N = MOESM(N+3).

---

## 2. What is inside the file

Two sheets matter.

**"Effect Sizes"**: 51 data rows, one per intervention arm. Column A is the arm name. Then study
metadata (Access, Sample Type, Array, Study Design: Allocation, Condition, Length of intervention,
Total subjects), then one column per biomarker holding the standardized effect size. Effect sizes are
age-residualized and divided by the s.d. of that biomarker in the Health and Retirement Study, so
they are in s.d. units and comparable across clocks.

**"Pvalues"**: the same 51 arms, same biomarker columns, holding the paired t-test P value for each
arm-by-biomarker cell. Column index 17 is the sample size used in the test.

The 16 "prominent clocks" are: Horvath1, Horvath2, Hannum, PCHorvath1, PCHorvath2, PCHannum,
DNAmEMRAge, OMICmAge, PhenoAge, GrimAgeV1, GrimAgeV2, PCPhenoAge, PCGrimAge, SystemsAge,
DunedinPoAm38, DunedinPACE.

---

## 3. The one calculation everything rests on

Their Figure 3d metric, which produces the verdict list in our Figure 1a:

> For each arm, take that arm's 16 prominent-clock effect sizes. Run a one-sample two-sided
> t-test (df = 15) against zero. Significant after their correction means P < 0.00625.

**This is their analysis, not ours.** The proof is that it reproduces their published sentence exactly:
19 nominal decreases, 13 after correction, 26 with no significant effect, 3 corrected increases.
(It gives 6 nominal increases where they report 5; their stated counts sum to 50 of 51 arms.)

---

## 4. Check it in Excel, in five minutes

1. Open the "Effect Sizes" sheet.
2. In an empty column, for the row `ART Therapy`, compute the mean of its 16 clock cells,
   and the sample s.d. of the same 16 cells.
3. `t = mean / (sd / SQRT(16))`
4. `P = T.DIST.2T(ABS(t), 15)`
5. You should get mean = −0.350, P < 0.0001.

Repeat for `Botanical Supplement` (copaiba): mean = −0.050, P = 0.0050, which is below 0.00625.
Then for `Exercise`: mean = −0.024, P = 0.256. Those three rows are the letter's whole argument.

To check the exercise claim directly: on the "Pvalues" sheet, take the `Exercise` row and
`COUNTIF(<the 108 biomarker cells>, "<0.05")`. It returns 0.

---

## 5. Or run our script

    cd dnam_clock_reanalysis
    python3 verify_all.py

42 named checks, each printing PASS or FAIL with the computed value. It re-derives every number in
the letter from the .xlsx, including the reference values that test the t-distribution code itself.
`league_table.py` prints the full 51-arm ranking; `randomized_contrasts.py` prints the contrast
analysis with confidence intervals.

---

## 6. Two provenance points worth checking yourself

1. **"Botanical Supplement" is copaiba essential oil, sponsored by dōTERRA.** This is inside the
   spreadsheet, in the IRB column of the Effect Sizes sheet: search that row for "copaiba". It reads
   "For the botanical supplement (copaiba essential oil) study... sponsored by doTERRA International."
2. **"Lacta Supplement" is not lactoferrin.** Its Regimen field reads "Daily 75mg dose of Trulacta
   freeze dried human milk supplement". The word lactoferrin appears nowhere in the workbook.

## 7. The one thing you cannot check from the .xlsx

The verbatim quotations from the article's Results, Discussion and Methods. Those are in the article
text; search the HTML for each quoted phrase before submission.
