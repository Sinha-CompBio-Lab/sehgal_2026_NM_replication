"""Regenerate the figure input files from the article's Supplementary Table 1.

Outputs
  fig1b_v2_data.csv         treated / control / difference for the three randomized diet trials
  exercise_contrast_data.csv exercise vs no exercise (CENTRAL factorial), all 16 clocks

fig_data_a.csv is a curated file: the numbers are computed from the source table but the
display labels and category flags are editorial. This script verifies its numbers rather
than regenerating them, so any discrepancy is reported.
"""
import numpy as np, pandas as pd
from scipy import stats

SRC = 'SupplementaryTable1_Sehgal2026.xlsx'
CLOCKS = ['Horvath1','Horvath2','Hannum','PCHorvath1','PCHorvath2','PCHannum','DNAmEMRAge',
          'OMICmAge','PhenoAge','GrimAgeV1','GrimAgeV2','PCPhenoAge','PCGrimAge','SystemsAge',
          'DunedinPoAm38','DunedinPACE']

E = pd.read_excel(SRC, sheet_name='Effect Sizes'); E = E.set_index(E.columns[0])
P = pd.read_excel(SRC, sheet_name='Pvalues');      P = P.set_index(P.columns[0]).loc[E.index]
N = P[[c for c in P.columns if 'NSubjects' in str(c)][0]].astype(float)

def se(arm, clock):
    """Standard error implied by the reported effect size, P value and n."""
    d, p = float(E.loc[arm, clock]), float(P.loc[arm, clock])
    return abs(d) / stats.t.isf(p / 2, N[arm] - 1)

def pool(est, ses):
    w = 1 / np.asarray(ses) ** 2
    e = (w * np.asarray(est)).sum() / w.sum()
    return e, np.sqrt(1 / w.sum())

# ---- panel b: three randomized parallel-group diet trials -------------------
PAIRS = [('DIRECT-PLUS', 'Green Mediterranean Diet', 'Healthy Guidelines Diet'),
         ('CENTRAL',     'Low Carb Diet',            'Low Fat Diet'),
         ('Twin study',  'Vegan Diet',               'Omnivore Diet')]
rows = []
for trial, treated, control in PAIRS:
    for c in CLOCKS:
        dt, dc = float(E.loc[treated, c]), float(E.loc[control, c])
        st, sc = se(treated, c), se(control, c)
        diff, sd = dt - dc, np.hypot(st, sc)
        rows.append(dict(trial=trial, clock=c,
                         n_treated=int(N[treated]), n_control=int(N[control]),
                         treated=dt, treated_lo=dt - 1.96*st, treated_hi=dt + 1.96*st,
                         treated_p=float(P.loc[treated, c]),
                         control=dc, control_lo=dc - 1.96*sc, control_hi=dc + 1.96*sc,
                         control_p=float(P.loc[control, c]),
                         diff=diff, diff_lo=diff - 1.96*sd, diff_hi=diff + 1.96*sd,
                         diff_p=2 * stats.norm.sf(abs(diff / sd))))
pd.DataFrame(rows).to_csv('fig1b_v2_data.csv', index=False)
print('wrote fig1b_v2_data.csv')

# ---- panel c and Extended Data Fig. 1: CENTRAL factorial exercise contrast --
EX = [('Low Carb Diet + Exercise', 'Low Carb Diet + no Exercise'),
      ('Low Fat Diet + Exercise',  'Low Fat Diet + no Exercise')]
rows = []
for c in CLOCKS:
    te, ts = pool([float(E.loc[a, c]) for a, _ in EX], [se(a, c) for a, _ in EX])
    ce, cs = pool([float(E.loc[b, c]) for _, b in EX], [se(b, c) for _, b in EX])
    ds = [float(E.loc[a, c]) - float(E.loc[b, c]) for a, b in EX]
    ss = [np.hypot(se(a, c), se(b, c)) for a, b in EX]
    de, dsd = pool(ds, ss)
    rows.append(dict(clock=c,
                     treated=te, treated_lo=te - 1.96*ts, treated_hi=te + 1.96*ts,
                     control=ce, control_lo=ce - 1.96*cs, control_hi=ce + 1.96*cs,
                     diff=de, diff_lo=de - 1.96*dsd, diff_hi=de + 1.96*dsd,
                     diff_p=2 * stats.norm.sf(abs(de / dsd)),
                     sub1=ds[0], sub1_lo=ds[0] - 1.96*ss[0], sub1_hi=ds[0] + 1.96*ss[0],
                     sub2=ds[1], sub2_lo=ds[1] - 1.96*ss[1], sub2_hi=ds[1] + 1.96*ss[1]))
pd.DataFrame(rows).to_csv('exercise_contrast_data.csv', index=False)
print('wrote exercise_contrast_data.csv')

# ---- verify the curated panel-a file against the source --------------------
M = E[CLOCKS].astype(float)
t, p = stats.ttest_1samp(M.values, 0.0, axis=1)
arm = pd.DataFrame({'mean': M.mean(axis=1), 'p': p, 'n': N}, index=E.index)
arm['group'] = np.where(arm.p < 0.00625, np.where(arm['mean'] < 0, 'dec', 'inc'), 'ns')

a = pd.read_csv('fig_data_a.csv')
bad = 0
for _, r in a.iterrows():
    hit = arm[(np.isclose(arm['mean'], r['mean'], atol=5e-4)) & (arm.n == r['n'])]
    if hit.empty or r['group'] not in set(hit['group']):
        print(f"  MISMATCH: {r['label']} (mean={r['mean']}, n={r['n']}, group={r['group']})"); bad += 1
print(f'fig_data_a.csv: {len(a)} rows checked against the source table, {bad} mismatches')
print(f"counts at P < 0.00625 -> decreases {(arm.group=='dec').sum()}, "
      f"increases {(arm.group=='inc').sum()}, no effect {(arm.group=='ns').sum()}")
