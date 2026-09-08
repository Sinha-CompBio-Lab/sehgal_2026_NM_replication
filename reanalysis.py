"""
Reanalysis of Sehgal et al., Nat Med 2026 (s41591-026-04562-9)
"Responsiveness of epigenetic aging biomarkers to longevity interventions in humans"

Input: the paper's own Supplementary Table 1 (MOESM4), sheets 'Effect Sizes' and 'Pvalues'.
Download:
  https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41591-026-04562-9/MediaObjects/41591_2026_4562_MOESM4_ESM.xlsx

Three analyses:
  A. Randomized between-arm contrasts (the estimand the paper does not report)
  B. Non-independence of the 51 "interventions" -> cluster-level re-test of the 16 clocks
  C. Multiple-testing / specificity audit at the study x biomarker level
"""
import openpyxl, math, statistics as st
from math import lgamma, log, exp, sqrt, erfc

XLSX = 'SupplementaryTable1_Sehgal2026.xlsx'
CLOCKS16 = ['Horvath1','Horvath2','Hannum','PCHorvath1','PCHorvath2','PCHannum',
            'DNAmEMRAge','OMICmAge','PhenoAge','GrimAgeV1','GrimAgeV2','PCPhenoAge',
            'PCGrimAge','SystemsAge','DunedinPoAm38','DunedinPACE']
HEADLINE = ['GrimAgeV2','PCPhenoAge','PCGrimAge','SystemsAge','DunedinPACE']

# ---------- t distribution helpers (no scipy dependency) ----------
def _betacf(a,b,x):
    MAXIT,EPS,FPMIN=300,3e-16,1e-300
    qab,qap,qam=a+b,a+1,a-1
    c,d=1.0,1-qab*x/qap
    d=1/(d if abs(d)>FPMIN else FPMIN); h=d
    for m in range(1,MAXIT+1):
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2)); d=1+aa*d
        d=1/(d if abs(d)>FPMIN else FPMIN); c=1+aa/c
        c=c if abs(c)>FPMIN else FPMIN; h*=d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2)); d=1+aa*d
        d=1/(d if abs(d)>FPMIN else FPMIN); c=1+aa/c
        c=c if abs(c)>FPMIN else FPMIN; de=d*c; h*=de
        if abs(de-1)<EPS: break
    return h
def _betai(a,b,x):
    if x<=0: return 0.0
    if x>=1: return 1.0
    bt=exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1-x))
    return bt*_betacf(a,b,x)/a if x<(a+1)/(a+b+2) else 1-bt*_betacf(b,a,1-x)/b
def t_pval(t,df): return _betai(df/2,0.5,df/(df+t*t))
def t_from_p(p,df):
    lo,hi=0.0,300.0
    for _ in range(200):
        mid=(lo+hi)/2
        lo,hi=(mid,hi) if t_pval(mid,df)>p else (lo,mid)
    return (lo+hi)/2
def z_pval(z): return erfc(abs(z)/sqrt(2))
def one_sample_t(v):
    n=len(v); m=st.mean(v); s=st.stdev(v); t=m/(s/sqrt(n))
    return m,t,t_pval(t,n-1),n
def bh_count(ps,q=0.05):
    s=sorted(ps); m=len(s); k=0
    for i,p in enumerate(s,1):
        if p<=i/m*q: k=i
    return k

# ---------- load ----------
wb=openpyxl.load_workbook(XLSX,data_only=True)
rE=[r for r in wb['Effect Sizes'].iter_rows(values_only=True) if r[0] is not None]
hE=list(rE[0]); E={r[0]:r for r in rE[1:]}
rP=[r for r in wb['Pvalues'].iter_rows(values_only=True) if any(c is not None for c in r)]
hP=list(rP[0]); P={r[0]:r for r in rP[1:]}
ARMS=[r[0] for r in rE[1:]]
N_COL=17          # NSubjects_ttest in 'Pvalues'
FOLDER_COL=2      # Folder.Name
BIOMARKERS=hP[19:]

def eff(arm,c): return E[arm][hE.index(c)]
def pv(arm,c):  return P[arm][hP.index(c)]
def nsub(arm):  return P[arm][N_COL]
def se(arm,c):
    d,p,n=eff(arm,c),pv(arm,c),nsub(arm)
    t=t_from_p(p,n-1)
    return abs(d)/t if t>1e-6 else None

# ================= A. randomized between-arm contrasts =================
print('='*78); print('A. RANDOMIZED BETWEEN-ARM CONTRASTS (arms present in the paper\'s own database)')
print('='*78)
PAIRS=[('Green Mediterranean Diet','Healthy Guidelines Diet','DIRECT-PLUS RCT, 18 mo'),
       ('Mediterranean Diet 1','Healthy Guidelines Diet','DIRECT-PLUS RCT, 18 mo'),
       ('Low Carb Diet','Low Fat Diet','CENTRAL RCT, 18 mo'),
       ('Vegan Diet','Omnivore Diet','Twin diet RCT, 8 wk'),
       ('HBOT - High Pressure','HBOT - Mild Pressure','HBOT, randomized triple-blind')]
for tr,ct,label in PAIRS:
    print(f'\n{tr}  minus  {ct}   [{label}]  n={nsub(tr):.0f} vs {nsub(ct):.0f}')
    for c in HEADLINE:
        d1,d2=eff(tr,c),eff(ct,c); s1,s2=se(tr,c),se(ct,c)
        if s1 and s2:
            diff=d1-d2; z=diff/sqrt(s1**2+s2**2)
            star='*' if pv(ct,c)<0.05 else ' '
            print(f'   {c:13s} treated {d1:+.3f}   comparator {d2:+.3f}{star}   contrast {diff:+.3f}  z={z:+.2f}  p={z_pval(z):.3f}')
    print('   (* = comparator arm itself passes the paper\'s p<0.05 "responsive" threshold)')

CONTROL_ARMS=['Healthy Guidelines Diet','Low Fat Diet','Low Fat Diet + Exercise',
              'Low Fat Diet + no Exercise','Omnivore Diet','HBOT - Mild Pressure','HBT - Mild Pressure']
m16=lambda a: st.mean([eff(a,c) for c in CLOCKS16])
print(f'\nMean effect over 16 clocks: comparator/sham arms counted among the 51 = '
      f'{st.mean([m16(a) for a in CONTROL_ARMS]):+.4f}  (n={len(CONTROL_ARMS)})')
print(f'                             all other arms                          = '
      f'{st.mean([m16(a) for a in ARMS if a not in CONTROL_ARMS]):+.4f}  (n={len(ARMS)-len(CONTROL_ARMS)})')

# ================= B. non-independence of the 51 arms =================
print('\n'+'='*78); print('B. THE 51 "INTERVENTIONS" ARE NOT 51 INDEPENDENT OBSERVATIONS')
print('='*78)
NEST=['E-MTAB-8956','E-MTAB-12527','E-MTAB-4931','GSE191297','TruD_HBOT','TruD_HBT',
      'TruD_Rapa','TruDiagnosticSenolytics1','TruDiagnosticSenolytics3','TruD_SRW']
def clus(a):
    f=P[a][FOLDER_COL]
    for p in NEST:
        if f.startswith(p): return p
    return f
G={}
for a in ARMS: G.setdefault(clus(a),[]).append(a)
print(f'51 arms -> {len(G)} independent study clusters. Clusters contributing >1 arm:')
for k,v in sorted(G.items(),key=lambda x:-len(x[1])):
    if len(v)>1:
        print(f'   {k:28s} {len(v)} arms: ' + ', '.join(f'{a} (n={nsub(a):.0f})' for a in v))
print('\nOne-sample t-test on effect sizes, as published (n=51) vs collapsed to clusters (n=%d):'%len(G))
print(f"{'clock':15s} {'mean51':>8s} {'p51':>8s} | {'mean_cl':>8s} {'t':>7s} {'p_cl':>8s}   Me=6 cutoff 0.0083 | Bonf/16 0.0031")
for c in CLOCKS16:
    m0,_,p0,_=one_sample_t([eff(a,c) for a in ARMS])
    m1,t1,p1,_=one_sample_t([st.mean([eff(a,c) for a in g]) for g in G.values()])
    tag=('PASS Me' if p1<0.05/6 else 'fail Me') + ' | ' + ('PASS Bonf16' if p1<0.05/16 else 'fail Bonf16')
    print(f'{c:15s} {m0:+8.4f} {p0:8.4f} | {m1:+8.4f} {t1:7.3f} {p1:8.4f}   {tag}')

# ================= C. multiple testing / specificity =================
print('\n'+'='*78); print('C. MULTIPLE TESTING AND SPECIFICITY AT THE ARM x BIOMARKER LEVEL')
print('='*78)
cells=[pv(a,c) for a in ARMS for c in CLOCKS16 if isinstance(pv(a,c),(int,float))]
print(f'51 arms x 16 clocks = {len(cells)} tests: nominal p<0.05 = {sum(p<0.05 for p in cells)} '
      f'(null expectation {len(cells)*0.05:.0f}); BH-FDR 5% = {bh_count(cells)}')
allc=[pv(a,c) for a in ARMS for c in BIOMARKERS if isinstance(pv(a,c),(int,float))]
print(f'51 arms x {len(BIOMARKERS)} biomarkers = {len(allc)} tests: nominal p<0.05 = {sum(p<0.05 for p in allc)} '
      f'(null expectation {len(allc)*0.05:.0f}); BH-FDR 5% = {bh_count(allc)}')
print('\n=> There IS more signal than chance. The question is whether it is intervention signal.')
print('   Per-study specificity audit (all biomarkers, within-study BH-FDR at 5%):\n')
print(f"{'arm':32s} {'n':>4s} {'nom p<.05':>10s} {'expected':>9s} {'BH-FDR':>7s}")
for a in ['Lavendar Oil Supplement','Botanical Supplement','Lacta Supplement','Metformin',
          'Rapamycin dose 1','Rapamycin dose 3','DQ Senolytics 1','Follistatin Gene Therapy',
          'Green Mediterranean Diet','CR / IMF']:
    ps=[pv(a,c) for c in BIOMARKERS if isinstance(pv(a,c),(int,float))]
    print(f'{a:32s} {nsub(a):4.0f} {sum(p<0.05 for p in ps):10d} {len(ps)*0.05:9.1f} {bh_count(ps):7d}')
print('\nTissue of origin (clocks are whole-blood trained):')
for a in ARMS:
    ttype=E[a][hE.index('Sample Type')]
    if ttype!='Whole Blood': print(f'   {a:32s} {ttype}')

# ================= D. DIRECTION-AWARE AND STRATIFIED =================
print('\n'+'='*78); print('D. DIRECTION OF THE SIGNIFICANT RESULTS, AND WHERE IT COMES FROM')
print('='*78)
def bh_keys(items,q=0.05):
    s=sorted(items); m=len(s); k=0
    for i,(p,_) in enumerate(s,1):
        if p<=i/m*q: k=i
    return [x[1] for x in s[:k]]
def lchoose(n,k): return lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1)
def binom2(k,n):
    if n==0: return 1.0
    lp=lambda i: lchoose(n,i)-n*log(2); ref=lp(k)+1e-12
    return min(1.0,sum(exp(lp(i)) for i in range(n+1) if lp(i)<=ref))
def report(tag,arms,cols):
    it=[(pv(a,c),(a,c)) for a in arms for c in cols if isinstance(pv(a,c),(int,float))]
    nom=[k for p,k in it if p<0.05]; fdr=bh_keys(it)
    for nm,ks in [('nominal',nom),('BH-FDR',fdr)]:
        d=sum(1 for a,c in ks if eff(a,c)<0); u=len(ks)-d
        pct=f'{d/len(ks)*100:3.0f}% down' if ks else '   n/a   '
        print(f'  {tag:34s} {nm:8s} {len(ks):4d}/{len(it):5d}  down {d:4d}  up {u:4d}  {pct}  binom p={binom2(d,len(ks)):.1e}')

report('all 51 arms',ARMS,CLOCKS16)
report('all 51 arms, 108 biomarkers',ARMS,BIOMARKERS)
print()
DISEASE=['ART Therapy','Anti-TNF therapy 1a','Anti-TNF therapy 1b','Anti-TNF therapy 2',
         'Kidney Transplant','Kidney Dialysis','Metformin','Gastric Bypass']
report('disease-treatment arms (8)',DISEASE,CLOCKS16)
report('control / sham arms (7)',CONTROL_ARMS,CLOCKS16)
LONG=[a for a in ARMS if a not in DISEASE and a not in CONTROL_ARMS]
report('longevity-intervention arms (36)',LONG,CLOCKS16)
print()
for nm,sel in [('disease treatment',DISEASE),('control / sham',CONTROL_ARMS),('longevity interventions',LONG)]:
    v=[m16(a) for a in sel]
    print(f'  mean effect over 16 clocks: {nm:24s} n={len(v):2d}  {st.mean(v):+.4f}   median {st.median(v):+.4f}')

print('\nFigure 4b re-run, dropping strata (nominal p; paper Me cutoff 0.0083):')
DIS_CORE=DISEASE[:6]; CTRL_STRICT=['Healthy Guidelines Diet','Omnivore Diet','HBOT - Mild Pressure','HBT - Mild Pressure']
scen=[('all 51',[]),('-disease(6)',DIS_CORE),('-disease(8)',DISEASE),
      ('-ctrl(4)',CTRL_STRICT),('-dis6-ctrl4',DIS_CORE+CTRL_STRICT),('-dis8-ctrl7',DISEASE+CONTROL_ARMS)]
print(f"{'clock':15s}"+''.join(f'{n:>13s}' for n,_ in scen))
for c in CLOCKS16:
    row=f'{c:15s}'
    for _,ex in scen:
        _,_,p,_=one_sample_t([eff(a,c) for a in ARMS if a not in ex])
        row+=f"{('%.4f'%p)+('*' if p<0.05 else ' '):>13s}"
    print(row)
