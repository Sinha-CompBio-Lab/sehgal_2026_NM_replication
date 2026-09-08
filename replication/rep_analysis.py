"""Independent replication. Built from scratch: pandas + scipy + statsmodels.
Reads the article's Supplementary Table 1 directly. Imports nothing from prior work."""
import os
import pandas as pd, numpy as np
from scipy import stats
import statsmodels.api as sm
import warnings; warnings.filterwarnings('ignore')

SRC=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','SupplementaryTable1_Sehgal2026.xlsx')
E=pd.read_excel(SRC, sheet_name='Effect Sizes')
P=pd.read_excel(SRC, sheet_name='Pvalues')
E.columns=[str(c).strip() for c in E.columns]; P.columns=[str(c).strip() for c in P.columns]
armcol_E=E.columns[0]; armcol_P=P.columns[0]
E=E.set_index(armcol_E); P=P.set_index(armcol_P)
P=P.loc[E.index]                                   # align by NAME (row orders differ)
CL16=['Horvath1','Horvath2','Hannum','PCHorvath1','PCHorvath2','PCHannum','DNAmEMRAge','OMICmAge',
      'PhenoAge','GrimAgeV1','GrimAgeV2','PCPhenoAge','PCGrimAge','SystemsAge','DunedinPoAm38','DunedinPACE']
BIO=[c for c in P.columns if c in set(E.columns) and pd.api.types.is_numeric_dtype(P[c]) and P[c].notna().sum()>40]
BIO=[c for c in BIO if c not in ('NSubjects_ttest',)]
NCOL=[c for c in P.columns if 'NSubjects' in str(c)][0]
N=P[NCOL].astype(float)
_c=E['Condition']
COND={a: ('' if pd.isna(v) else str(v).strip()) for a,v in _c.items()}
DIS=[a for a in E.index if COND[a] not in ('','None','NA','nan','NaN')]
assert len(DIS)==10, f'disease stratum should be 10, got {len(DIS)}: {DIS}'
REF=['Healthy Guidelines Diet','Low Fat Diet','Low Fat Diet + Exercise','Low Fat Diet + no Exercise',
     'Omnivore Diet','HBOT - Mild Pressure','HBT - Mild Pressure']
H41=[a for a in E.index if a not in DIS]; H34=[a for a in H41 if a not in REF]
out=[]
def say(s=''): print(s); out.append(str(s))
def hdr(s): say('\n'+'='*76); say(s); say('='*76)

hdr('SANITY: data shape and alignment')
say(f'arms: {len(E)} | biomarker columns usable: {len(BIO)} | 16 clocks present: {all(c in E.columns for c in CL16)}')
say(f'row orders differ across sheets: {sum(a!=b for a,b in zip(E.index, pd.read_excel(SRC,sheet_name="Pvalues").iloc[:,0]))} of 51 mismatched (aligned by name)')
say(f'disease arms (non-empty Condition): {len(DIS)} | healthy: {len(H41)} | healthy minus reference: {len(H34)}')

# ---------------- ARGUMENT 2 ----------------
hdr('ARGUMENT 2 - ranking on the article Fig 3d metric')
M=E[CL16].astype(float)
t,p=stats.ttest_1samp(M.values, 0.0, axis=1)
arm=pd.DataFrame({'mean':M.mean(axis=1),'t':t,'p':p,'n':N}, index=E.index)
CUT=0.00625
arm['verdict']=np.where(arm.p<CUT, np.where(arm['mean']<0,'DEC','INC'),
                 np.where(arm.p<0.05, np.where(arm['mean']<0,'dec_nom','inc_nom'),'ns'))
say(f"Q2.1 nominal decreases: {((arm['mean']<0)&(arm.p<0.05)).sum()} (predicted 19)")
say(f"     corrected decreases: {(arm.verdict=='DEC').sum()} (predicted 13)")
say(f"     nominal increases: {((arm['mean']>0)&(arm.p<0.05)).sum()} (predicted 6; article says 5)")
say(f"     corrected increases: {(arm.verdict=='INC').sum()} (predicted 3)")
say(f"     no significant effect: {(arm.p>=0.05).sum()} (predicted 26)")
dec=arm[arm.verdict=='DEC'].sort_values('mean')
say('\nQ2.2 the corrected decreases:')
for a,r in dec.iterrows():
    tag='DISEASE' if a in DIS else ('COMPARATOR' if a in REF else '')
    say(f'   {a:34s} mean={r["mean"]:+.3f} P={r.p:.4g} n={r.n:.0f}  {tag}')
say(f'   -> disease: {sum(a in DIS for a in dec.index)} (pred 6) | comparator: {sum(a in REF for a in dec.index)} (pred 1)')
say('\nQ2.3 exercise and CR:')
for a in ['Exercise','CR / IMF']:
    say(f'   {a:12s} mean={arm.loc[a,"mean"]:+.3f} P={arm.loc[a,"p"]:.3f} verdict={arm.loc[a,"verdict"]}')
ex_hits=int((P.loc['Exercise',BIO].astype(float)<0.05).sum())
say(f'   Exercise biomarkers at nominal P<0.05: {ex_hits} of {len(BIO)} (expected by chance {len(BIO)*0.05:.1f})')

hdr('Q2.4 ROBUSTNESS of the ranking to metric choice')
GEN2=['PhenoAge','GrimAgeV1','GrimAgeV2','PCPhenoAge','PCGrimAge','SystemsAge','DunedinPoAm38','DunedinPACE']
def variant(name, mat, test='t'):
    if test=='t':
        tt,pp=stats.ttest_1samp(mat.values,0,axis=1); est=mat.mean(axis=1)
    elif test=='wilcoxon':
        pp=np.array([stats.wilcoxon(r)[1] if np.any(r!=0) else 1.0 for r in mat.values]); est=mat.median(axis=1)
    elif test=='median':
        tt,pp=stats.ttest_1samp(mat.values,0,axis=1); est=mat.median(axis=1)
    elif test=='trim':
        tt,pp=stats.ttest_1samp(stats.trimboth(mat.values,0.125,axis=1),0,axis=1); est=mat.mean(axis=1)
    d=pd.DataFrame({'est':est,'p':pp},index=mat.index)
    sig=d[(d.p<CUT)&(d.est<0)]
    return name, sorted(sig.index), d
VAR=[variant('16 clocks, t-test (their metric)', M, 't'),
     variant('16 clocks, Wilcoxon signed-rank', M, 'wilcoxon'),
     variant('16 clocks, trimmed mean', M, 'trim'),
     variant('Gen2+ clocks only (8)', E[GEN2].astype(float), 't'),
     variant('all 108 biomarkers', E[BIO].astype(float), 't')]
KEY=['Botanical Supplement','Lacta Supplement','Ketamine','Exercise','CR / IMF','Low Fat Diet + no Exercise']
say(f'{"variant":34s} {"#DEC":>5s}  ' + ' '.join(f'{k[:11]:>12s}' for k in KEY))
for name,sig,d in VAR:
    marks=' '.join(f'{("PASS" if k in sig else "-"):>12s}' for k in KEY)
    say(f'{name:34s} {len(sig):5d}  {marks}')
say('\nQ2.5 copaiba ("Botanical Supplement") passes in '
    f'{sum("Botanical Supplement" in s for _,s,_ in VAR)} of {len(VAR)} variants')
say(f'     exercise/CR pass in {sum(("Exercise" in s) or ("CR / IMF" in s) for _,s,_ in VAR)} of {len(VAR)} variants')

# ---------------- ARGUMENT 1 (testable parts) ----------------
hdr('ARGUMENT 1 - does effect size track study quality / evidence strength?')
alloc=E['Study Design: Allocation'].astype(str)
pub=E['Link to main publication (PMID)'].notna() & (E['Link to main publication (PMID)'].astype(str)!='None')
qual=pd.DataFrame({'mean':arm['mean'],'n':N,
                   'randomized':(alloc=='Randomized').astype(int),
                   'published':pub.astype(int),
                   'disease':[1 if a in DIS else 0 for a in E.index]}, index=E.index)
for g,lab in [('randomized','randomized'),('published','published')]:
    a1=qual[qual[g]==1]['mean']; a0=qual[qual[g]==0]['mean']
    u=stats.mannwhitneyu(a1,a0)
    say(f'Q1.4 {lab:12s} yes n={len(a1):2d} mean={a1.mean():+.4f} | no n={len(a0):2d} mean={a0.mean():+.4f} | MWU P={u.pvalue:.3f}')
r,pv=stats.spearmanr(qual['n'],qual['mean'])
say(f'     Spearman(n, effect) rho={r:+.3f} P={pv:.3f}  [negative rho would mean bigger studies -> bigger decreases]')
# evidence tier
tier=(qual.randomized+qual.published+(qual.n>=30).astype(int))
r2,p2=stats.spearmanr(tier, qual['mean'])
say(f'Q1.3 evidence tier (0-3) vs effect: Spearman rho={r2:+.3f} P={p2:.3f}')
say(f'     tier means: ' + ', '.join(f'{k}:{v:+.4f}(n={(tier==k).sum()})' for k,v in qual.groupby(tier)['mean'].mean().items()))

# ---------------- ARGUMENT 3 ----------------
hdr('ARGUMENT 3 - disease vs healthy')
def stratum(arms):
    m=E.loc[arms,CL16].astype(float)
    tt,pp=stats.ttest_1samp(m.values,0,axis=0)
    return pd.DataFrame({'mean':m.mean(axis=0),'p':pp},index=CL16)
s41=stratum(H41); s10=stratum(DIS); s34=stratum(H34)
say('Q3.1 reproduce article stratum values:')
for c,exp_m,exp_p in [('PCPhenoAge',-0.320,0.0057),('PCGrimAge',-0.223,0.0022),('SystemsAge',-0.213,0.0105)]:
    say(f'   disease {c:12s} ours {s10.loc[c,"mean"]:+.4f} P={s10.loc[c,"p"]:.4f} | article {exp_m:+.3f} P={exp_p}')
for c,exp_m,exp_p in [('DunedinPACE',-0.074,0.0143),('PCGrimAge',-0.043,0.0191)]:
    say(f'   healthy {c:12s} ours {s41.loc[c,"mean"]:+.4f} P={s41.loc[c,"p"]:.4f} | article {exp_m:+.3f} P={exp_p}')
say(f'\nQ3.2 clocks passing 0.00833 in healthy-41: {list(s41[s41.p<0.00833].index) or "NONE"}')
say(f'     clocks passing 0.00833 in healthy-34: {list(s34[s34.p<0.00833].index) or "NONE"}')
say(f'     clocks passing 0.0125 in disease-10:  {list(s10[s10.p<0.0125].index)}')
say(f'     mean over 16 clocks: H34={E.loc[H34,CL16].astype(float).mean(axis=1).mean():+.4f}  DIS={E.loc[DIS,CL16].astype(float).mean(axis=1).mean():+.4f}')
say('\nQ3.3 leave-one-out on the disease stratum (clocks passing 0.0125):')
mins=99
for drop in DIS:
    sub=[a for a in DIS if a!=drop]; s=stratum(sub); k=(s.p<0.0125).sum(); mins=min(mins,k)
    if k<=2 or drop=='ART Therapy': say(f'   drop {drop:26s} -> {k} clocks pass')
say(f'   minimum across all leave-one-out: {mins} clocks (predicted >=2)')
say('\nQ3.4/3.5/3.6 confounding and alternative tests:')
say(f'   n: disease median={N[DIS].median():.0f} vs healthy median={N[H41].median():.0f}; MWU P={stats.mannwhitneyu(N[DIS],N[H41]).pvalue:.3f}')
for c in ['PCGrimAge','SystemsAge','DunedinPACE','PCPhenoAge']:
    a=E.loc[DIS,c].astype(float); b=E.loc[H41,c].astype(float)
    say(f'   {c:12s} Mann-Whitney disease vs healthy P={stats.mannwhitneyu(a,b).pvalue:.4f}')
X=sm.add_constant(pd.DataFrame({'disease':qual['disease'],'logn':np.log(qual['n']),
                                'dur':pd.to_numeric(E['Length of intervention (Days)'],errors='coerce').fillna(180)}))
fit=sm.OLS(qual['mean'],X).fit()
say(f'   OLS: effect ~ disease + log(n) + duration | disease beta={fit.params["disease"]:+.4f} P={fit.pvalues["disease"]:.4f}')

# ---------------- ARGUMENT 4 ----------------
hdr('ARGUMENT 4 - randomized contrasts')
def se_from_p(d,p,n):
    if not np.isfinite(p) or p<=0 or p>=1 or n<2: return np.nan
    tt=stats.t.isf(p/2, n-1)
    return abs(d)/tt if tt>0 else np.nan
say('Q4.1 SE round-trip check across all arm x clock cells:')
err=[]
for a in E.index:
    for c in CL16:
        d=float(E.loc[a,c]); p=float(P.loc[a,c]); n=float(N[a])
        s=se_from_p(d,p,n)
        if np.isfinite(s) and s>0:
            back=2*stats.t.sf(abs(d)/s, n-1); err.append(abs(back-p))
say(f'   max |recovered P - reported P| = {max(err):.2e} over {len(err)} cells (predicted <1e-9)')
PAIRS=[('Green Mediterranean Diet','Healthy Guidelines Diet'),('Low Carb Diet','Low Fat Diet'),
       ('Vegan Diet','Omnivore Diet')]
EXP=[('Low Carb Diet + Exercise','Low Carb Diet + no Exercise'),('Low Fat Diet + Exercise','Low Fat Diet + no Exercise')]
PUBEST={'DunedinPACE':-0.0891,'PCGrimAge':-0.07843,'SystemsAge':-0.05859,'GrimAgeV2':-0.0814,'PCPhenoAge':-0.08876}
def contrasts(pairs,c):
    ds,vs=[],[]
    for tr,ct in pairs:
        s1=se_from_p(float(E.loc[tr,c]),float(P.loc[tr,c]),float(N[tr]))
        s2=se_from_p(float(E.loc[ct,c]),float(P.loc[ct,c]),float(N[ct]))
        ds.append(float(E.loc[tr,c])-float(E.loc[ct,c])); vs.append(s1**2+s2**2)
    return np.array(ds), np.array(vs)
def pool(ds,vs,method='fixed'):
    w=1/vs
    if method=='fixed':
        est=(w*ds).sum()/w.sum(); se=np.sqrt(1/w.sum())
        return est,se,2*stats.norm.sf(abs(est/se))
    Q=(w*(ds-(w*ds).sum()/w.sum())**2).sum(); k=len(ds)
    C=w.sum()-(w**2).sum()/w.sum(); tau2=max(0,(Q-(k-1))/C)
    w2=1/(vs+tau2); est=(w2*ds).sum()/w2.sum(); se=np.sqrt(1/w2.sum())
    if method=='hk':
        q=((w2*(ds-est)**2).sum())/(k-1); se=np.sqrt(q/w2.sum())
        return est,se,2*stats.t.sf(abs(est/se),k-1)
    return est,se,2*stats.norm.sf(abs(est/se))
say('\nQ4.2/4.3/4.4 pooled contrasts, three methods:')
say(f'{"clock":13s} {"method":10s} {"est":>8s} {"95% CI":>20s} {"P":>7s}  excl. published?')
excl_count={'fixed':0,'random':0,'hk':0}
for c in PUBEST:
    ds,vs=contrasts(PAIRS,c)
    for m in ['fixed','random','hk']:
        est,se,pv=pool(ds,vs,m)
        crit=stats.t.isf(0.025,len(ds)-1) if m=='hk' else 1.96
        lo,hi=est-crit*se, est+crit*se
        ex = not (lo<=PUBEST[c]<=hi); excl_count[m]+=ex
        say(f'{c:13s} {m:10s} {est:+8.4f} [{lo:+.3f}, {hi:+.3f}] {pv:7.3f}  {"YES" if ex else "no"}')
say(f'\n   exclusions: fixed {excl_count["fixed"]}/5, random {excl_count["random"]}/5, Hartung-Knapp {excl_count["hk"]}/5 (predicted 2 for fixed)')
say('\nQ4.5 POWER: minimum detectable difference (80% power, alpha .05, 2-sided):')
for c in PUBEST:
    ds,vs=contrasts(PAIRS,c); se=np.sqrt(1/ (1/vs).sum())
    mde=(1.96+0.84)*se
    say(f'   {c:13s} pooled SE={se:.4f}  MDE={mde:.4f}  published effect={abs(PUBEST[c]):.4f}  '
        f'{"MDE < published (informative)" if mde<abs(PUBEST[c]) else "MDE > published (underpowered)"}')
say('\nQ4.6 exercise factorial contrasts, all methods:')
worst={}
for m in ['fixed','random','hk']:
    ps=[]
    for c in CL16:
        ds,vs=contrasts(EXP,c); ps.append(pool(ds,vs,m)[2])
    worst[m]=min(ps)
    say(f'   {m:10s} smallest P across 16 clocks = {min(ps):.3f}  ({"all null" if min(ps)>0.05 else "SOME SIGNIFICANT"})')
open('replication_output.txt','w').write('\n'.join(out))
print('\n[saved replication_output.txt]')
