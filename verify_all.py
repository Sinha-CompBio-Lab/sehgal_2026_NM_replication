# -*- coding: utf-8 -*-
"""Complete verification of every quantitative claim in the Matters Arising.
Prints PASS/FAIL for each. Source: the paper's own Supplementary Table 1."""
import openpyxl, statistics as st
from math import lgamma, log, exp, sqrt, erfc
def _bcf(a,b,x):
    MAXIT,EPS,FPMIN=300,3e-16,1e-300
    qab,qap,qam=a+b,a+1,a-1; c,d=1.0,1-qab*x/qap
    d=1/(d if abs(d)>FPMIN else FPMIN); h=d
    for m in range(1,MAXIT+1):
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2)); d=1+aa*d; d=1/(d if abs(d)>FPMIN else FPMIN)
        c=1+aa/c; c=c if abs(c)>FPMIN else FPMIN; h*=d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2)); d=1+aa*d; d=1/(d if abs(d)>FPMIN else FPMIN)
        c=1+aa/c; c=c if abs(c)>FPMIN else FPMIN; de=d*c; h*=de
        if abs(de-1)<EPS: break
    return h
def _bi(a,b,x):
    if x<=0: return 0.0
    if x>=1: return 1.0
    bt=exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1-x))
    return bt*_bcf(a,b,x)/a if x<(a+1)/(a+b+2) else 1-bt*_bcf(b,a,1-x)/b
def tp(t,df): return _bi(df/2,0.5,df/(df+t*t))
def t_from_p(p,df):
    lo,hi=0.0,1e9
    for _ in range(600):
        mid=(lo+hi)/2; lo,hi=(mid,hi) if tp(mid,df)>p else (lo,mid)
    return (lo+hi)/2
def ost(v):
    n=len(v); m=st.mean(v); s=st.stdev(v); t=m/(s/sqrt(n)); return m,tp(t,n-1)
zp=lambda z: erfc(abs(z)/sqrt(2))

wb=openpyxl.load_workbook('SupplementaryTable1_Sehgal2026.xlsx',data_only=True)
rE=[r for r in wb['Effect Sizes'].iter_rows(values_only=True) if r[0] is not None]
h=list(rE[0]); E={r[0]:r for r in rE[1:]}; ARMS=[r[0] for r in rE[1:]]
rP=[r for r in wb['Pvalues'].iter_rows(values_only=True) if any(c is not None for c in r)]
hP=list(rP[0]); P={r[0]:r for r in rP[1:]}
eff=lambda a,c: E[a][h.index(c)]; pv=lambda a,c: P[a][hP.index(c)]; nn=lambda a: P[a][17]
CL16=['Horvath1','Horvath2','Hannum','PCHorvath1','PCHorvath2','PCHannum','DNAmEMRAge','OMICmAge',
      'PhenoAge','GrimAgeV1','GrimAgeV2','PCPhenoAge','PCGrimAge','SystemsAge','DunedinPoAm38','DunedinPACE']
BIO=hP[19:]
res=[]
def check(name,cond,detail=''):
    res.append((name,bool(cond),detail))
    print(('PASS ' if cond else '*** FAIL ')+name+('  |  '+detail if detail else ''))

# --- 0. t-distribution implementation self-test against reference values ---
check('t-dist exact: t=2.131 df=15 -> P=0.0500',abs(tp(2.131,15)-0.0500)<5e-4,f'{tp(2.131,15):.4f}')
check('t-dist exact: t=2.042 df=30 -> P=0.0500',abs(tp(2.042,30)-0.0500)<5e-4,f'{tp(2.042,30):.4f}')
check('t-dist exact: t=3.000 df=15 -> P=0.0090',abs(tp(3.0,15)-0.0090)<5e-4,f'{tp(3.0,15):.4f}')
check('t_from_p inverts tp',abs(t_from_p(tp(2.5,15),15)-2.5)<1e-6)

# --- 1. arm-level metric & certified lists ---
arm={a:ost([eff(a,c) for c in CL16]) for a in ARMS}
CUT=0.00625
dec_nom=[a for a in ARMS if arm[a][0]<0 and arm[a][1]<0.05]
dec_mtc=[a for a in ARMS if arm[a][0]<0 and arm[a][1]<CUT]
inc_mtc=[a for a in ARMS if arm[a][0]>0 and arm[a][1]<CUT]
none_=[a for a in ARMS if arm[a][1]>=0.05]
check('paper counts: 19 nominal decreases',len(dec_nom)==19,f'{len(dec_nom)}')
check('paper counts: 13 corrected decreases',len(dec_mtc)==13,f'{len(dec_mtc)}')
check('paper counts: 3 corrected increases',len(inc_mtc)==3,f'{len(inc_mtc)}')
inc_nom=[a for a in ARMS if arm[a][0]>0 and arm[a][1]<0.05]
check('nominal increases = 6, NOT the 5 the paper reports',len(inc_nom)==6,
      f'{len(inc_nom)}; paper states 5, and 19+5+26=50 of 51 arms')
check('our counts sum to 51',len(dec_nom)+len(inc_nom)+len(none_)==51,
      f'{len(dec_nom)}+{len(inc_nom)}+{len(none_)}={len(dec_nom)+len(inc_nom)+len(none_)}')
check('paper counts: 26 with no effect',len(none_)==26,f'{len(none_)}')
TAB_DEC={'ART Therapy':(-0.350,183),'Anti-TNF therapy 1b':(-0.299,24),'Metformin':(-0.290,6),
 'Anti-TNF therapy 1a':(-0.234,26),'Anti-TNF therapy 2':(-0.187,51),'Gastric Bypass':(-0.140,12),
 'Lacta Supplement':(-0.100,19),'Vegan Diet':(-0.092,13),'Green Mediterranean Diet':(-0.087,87),
 'Ketamine':(-0.086,20),'Mediterranean Diet 1':(-0.050,81),'Botanical Supplement':(-0.050,50),
 'Low Fat Diet + no Exercise':(-0.049,30)}
check('Table 1 decrease membership',set(TAB_DEC)==set(dec_mtc),str(set(dec_mtc)^set(TAB_DEC)))
ok=all(abs(arm[a][0]-m)<5e-4 and nn(a)==n for a,(m,n) in TAB_DEC.items())
check('Table 1 decrease means & n exact',ok)
TAB_INC={'Kidney Dialysis':(0.270,12,0.0004),'Rapamycin dose 3':(0.088,7,0.0026),'Buckwheat Extract Supplement':(0.034,41,0.0007)}
check('Table 1 increase membership',set(TAB_INC)==set(inc_mtc))
ok=all(abs(arm[a][0]-m)<5e-4 and nn(a)==n and abs(arm[a][1]-p)<2e-4 for a,(m,n,p) in TAB_INC.items())
check('Table 1 increase means, n, P exact',ok)
pchecks={'Gastric Bypass':0.0003,'Lacta Supplement':0.0005,'Vegan Diet':0.0001,'Ketamine':0.0028,
 'Mediterranean Diet 1':0.0002,'Botanical Supplement':0.0050,'Low Fat Diet + no Exercise':0.0004}
ok=all(abs(arm[a][1]-p)<2e-4 for a,p in pchecks.items())
check('Table 1 quoted P values exact',ok)
sub4=all(arm[a][1]<1e-4 for a in ['ART Therapy','Anti-TNF therapy 1b','Metformin','Anti-TNF therapy 1a','Anti-TNF therapy 2','Green Mediterranean Diet'])
check('Table 1 "<10⁻⁴" entries all truly < 1e-4',sub4)
# no-effect line
ne={'Exercise':0.05,'CR / IMF':0.05,'Diet + Exercise':0.05,'Semaglutide':0.05,'Rapamycin dose 1':0.05}
check('no-effect list: all arm-level P >= 0.05',all(arm[a][1]>=0.05 for a in ne),
      '; '.join(f'{a} P={arm[a][1]:.3f}' for a in ne))
check('semaglutide +0.056 / rapamycin-low +0.024',abs(arm['Semaglutide'][0]-0.056)<5e-4 and abs(arm['Rapamycin dose 1'][0]-0.024)<5e-4)
# exercise biomarker sweep
ex=[pv('Exercise',c) for c in BIO if isinstance(pv('Exercise',c),(int,float))]
check('exercise: 0 of 108 biomarkers nominal',len(ex)==108 and sum(p<0.05 for p in ex)==0,
      f'tests={len(ex)}, nominal sig={sum(p<0.05 for p in ex)}, expected 5% = {len(ex)*0.05:.1f}')
check('exercise arm is randomized',str(E['Exercise'][h.index('Study Design: Allocation')])=='Randomized')
check('CR/IMF n=58 randomized',nn('CR / IMF')==58 and str(E['CR / IMF'][h.index('Study Design: Allocation')])=='Randomized')

# --- 2. strata ---
DIS=[r[0] for r in rE[1:] if str(r[h.index('Condition')]) not in ('None','NA','nan')]
REF=['Healthy Guidelines Diet','Low Fat Diet','Low Fat Diet + Exercise','Low Fat Diet + no Exercise',
     'Omnivore Diet','HBOT - Mild Pressure','HBT - Mild Pressure']
H41=[a for a in ARMS if a not in DIS]; H34=[a for a in H41 if a not in REF]
check('disease stratum = 10 arms (their Condition column)',len(DIS)==10,str(sorted(DIS)))
check('healthy = 41; healthy minus reference = 34',len(H41)==41 and len(H34)==34)
s41={c:ost([eff(a,c) for a in H41]) for c in CL16}
s10={c:ost([eff(a,c) for a in DIS]) for c in CL16}
s34={c:ost([eff(a,c) for a in H34]) for c in CL16}
check('healthy-41: no clock passes their Me cutoff 0.00833',all(s41[c][1]>=0.00833 for c in CL16),
      'min P = %.4f'%min(s41[c][1] for c in CL16))
d10pass=[c for c in CL16 if s10[c][1]<0.0125]
check('disease-10: exactly 5 clocks pass 0.0125',len(d10pass)==5,str(d10pass))
check('reproduce paper: disease PCPhenoAge -0.320 P=0.0057',abs(s10['PCPhenoAge'][0]+0.320)<5e-4 and abs(s10['PCPhenoAge'][1]-0.0057)<2e-4)
check('reproduce paper: disease PCGrimAge -0.223 P=0.0022',abs(s10['PCGrimAge'][0]+0.2226)<5e-4 and abs(s10['PCGrimAge'][1]-0.0022)<2e-4)
check('reproduce paper: healthy DunedinPACE -0.074 P=0.0143',abs(s41['DunedinPACE'][0]+0.0744)<5e-4 and abs(s41['DunedinPACE'][1]-0.0143)<2e-4)
check('H34 best clock: PCGrimAge -0.048 P=0.026',abs(s34['PCGrimAge'][0]+0.0479)<5e-4 and abs(s34['PCGrimAge'][1]-0.0261)<2e-3)
check('H34 nothing passes 0.00833',all(s34[c][1]>=0.00833 for c in CL16))
m16=lambda arms: st.mean([st.mean([eff(a,c) for c in CL16]) for a in arms])
check('mean16: H34 = -0.019, DIS = -0.141',abs(m16(H34)+0.019)<1.5e-3 and abs(m16(DIS)+0.141)<1.5e-3,
      f'H34={m16(H34):+.4f} DIS={m16(DIS):+.4f}')
# transplant/dialysis
check('kidney transplant -0.123 P=0.041',abs(arm['Kidney Transplant'][0]+0.123)<1e-3 and abs(arm['Kidney Transplant'][1]-0.0409)<2e-3)
check('20 arms coded Randomized',sum(1 for a in ARMS if str(E[a][h.index('Study Design: Allocation')])=='Randomized')==20)
# HDG control arm cells
hdg={'PCGrimAge':(-0.069,0.0071),'SystemsAge':(-0.056,0.0058),'DunedinPACE':(-0.131,0.0028)}
ok=all(abs(eff('Healthy Guidelines Diet',c)-m)<1e-3 and abs(pv('Healthy Guidelines Diet',c)-p)<2e-4 for c,(m,p) in hdg.items())
check('DIRECT-PLUS control arm: PCGrimAge/SystemsAge/DunedinPACE cells',ok,
      '; '.join(f"{c}: d={eff('Healthy Guidelines Diet',c):+.4f} P={pv('Healthy Guidelines Diet',c):.4f}" for c in hdg))
# unpublished status
pmcol=h.index('Link to main publication (PMID)')
check('Lacta & copaiba have no publication',E['Lacta Supplement'][pmcol] in (None,'None','NA') and E['Botanical Supplement'][pmcol] in (None,'None','NA'))
_irb=' '.join(str(x) for x in E['Botanical Supplement'] if isinstance(x,str))
check('copaiba + doTERRA named in the xlsx itself (IRB column)',
      'copaiba' in _irb.lower() and 'terra' in _irb.lower(),
      'directly verifiable from Supplementary Table 1, not only article prose')
check('"lactoferrin" appears nowhere in the workbook (Lacta = Trulacta human milk)',
      not any('lactoferrin' in str(c).lower() for sh in wb.sheetnames for r in wb[sh].iter_rows(values_only=True) for c in r),
      'Regimen reads: Daily 75mg dose of Trulacta freeze dried human milk supplement')
check('both rapamycin arms share one regimen string, so dose labels are not verifiable',
      E['Rapamycin dose 1'][h.index('Regimen')]==E['Rapamycin dose 3'][h.index('Regimen')],
      'draft therefore says "cohort", not "higher/lower dose"')
check('HBOT-High and HBT-High may share participants (excluded from primary contrasts)',
      nn('HBOT - High Pressure')==nn('HBT - High Pressure') and
      abs(E['HBOT - High Pressure'][h.index('Mean of Age of particpants')]-E['HBT - High Pressure'][h.index('Mean of Age of particpants')])<0.02,
      'both n=5, same dataset ID, mean ages differ by ~4 days')

# --- 3. randomized contrasts ---
def se(a,c):
    t=t_from_p(pv(a,c),nn(a)-1); return abs(eff(a,c))/t if t>1e-9 else None
PAIRS=[('Green Mediterranean Diet','Healthy Guidelines Diet'),('Low Carb Diet','Low Fat Diet'),
       ('Vegan Diet','Omnivore Diet')]
PAIRS5=PAIRS+[('HBOT - High Pressure','HBOT - Mild Pressure'),('HBT - High Pressure','HBT - Mild Pressure')]
HEAD={'GrimAgeV2':-0.0814,'PCPhenoAge':-0.08876,'PCGrimAge':-0.07843,'SystemsAge':-0.05859,'DunedinPACE':-0.0891}
def pool(pairs,c):
    arr=[]
    for t_,c_ in pairs:
        s1,s2=se(t_,c),se(c_,c)
        arr.append((eff(t_,c)-eff(c_,c), sqrt(s1**2+s2**2)))
    w=[1/s**2 for _,s in arr]
    est=sum(wi*d for wi,(d,_) in zip(w,arr))/sum(w); s=sqrt(1/sum(w))
    return est,est-1.96*s,est+1.96*s
pg=pool(PAIRS,'PCGrimAge'); sa=pool(PAIRS,'SystemsAge')
check('contrast PCGrimAge +0.001 [-0.058,+0.060], excludes -0.078',
      abs(pg[0]-0.0007)<2e-3 and not (pg[1]<=HEAD['PCGrimAge']<=pg[2]),
      f'est={pg[0]:+.4f} CI=[{pg[1]:+.3f},{pg[2]:+.3f}]')
_n_excl=sum(1 for c in HEAD if not (pool(PAIRS,c)[1]<=HEAD[c]<=pool(PAIRS,c)[2]))
check('exactly 2 of 5 clocks exclude the published estimate (draft says "two of the five")',
      _n_excl==2, f'{_n_excl} of 5')
_pg5=pool(PAIRS5,'PCGrimAge'); _sa5=pool(PAIRS5,'SystemsAge')
check('sensitivity: adding the 2 hyperbaric pairs changes no conclusion',
      not (_pg5[1]<=HEAD['PCGrimAge']<=_pg5[2]) and not (_sa5[1]<=HEAD['SystemsAge']<=_sa5[2]),
      f'PCGrimAge {_pg5[0]:+.4f} [{_pg5[1]:+.3f},{_pg5[2]:+.3f}]; SystemsAge {_sa5[0]:+.4f} [{_sa5[1]:+.3f},{_sa5[2]:+.3f}]')
check('contrast SystemsAge +0.008 [-0.033,+0.049], excludes -0.059',
      abs(sa[0]-0.008)<2e-3 and not (sa[1]<=HEAD['SystemsAge']<=sa[2]),
      f'est={sa[0]:+.4f} CI=[{sa[1]:+.3f},{sa[2]:+.3f}]')
# count nominal among 25 contrast cells
cnt=0; hit=None
for c in HEAD:
    for t_,c_ in PAIRS:
        d=eff(t_,c)-eff(c_,c); s_=sqrt(se(t_,c)**2+se(c_,c)**2)
        if zp(d/s_)<0.05: cnt+=1; hit=(c,t_, zp(d/s_))
check('exactly 1 of 15 contrast cells nominal (vegan DunedinPACE P=0.008)',cnt==1 and hit[0]=='DunedinPACE' and abs(hit[2]-0.008)<3e-3,str(hit))
check('join integrity: Effect Sizes and Pvalues share identical arm sets, no duplicates',
      set(E)==set(P) and len(E)==51 and len(P)==51)
# --- 4. FDR directionality (75/62) ---
def bh(items,q=.05):
    s=sorted(items); m=len(s); k=0
    for i,(p,_) in enumerate(s,1):
        if p<=i/m*q: k=i
    return [x[1] for x in s[:k]]
cells=[(pv(a,c),(a,c)) for a in ARMS for c in CL16 if isinstance(pv(a,c),(int,float))]
fdr=bh(cells); down=sum(1 for a,c in fdr if eff(a,c)<0)
check('75 clock cells survive BH; 62 decreases',len(fdr)==75 and down==62,f'{len(fdr)}/{down}')
print()
nf=sum(1 for _,ok,_ in res if not ok)
print(f'==== {len(res)} checks, {nf} failures ====')

# --- 5. CENTRAL factorial exercise contrasts (added in revision) ---
EXP=[('Low Carb Diet + Exercise','Low Carb Diet + no Exercise'),
     ('Low Fat Diet + Exercise','Low Fat Diet + no Exercise')]
worst=1.0
for c in CL16:
    arr=[(eff(t_,c)-eff(c_,c),sqrt(se(t_,c)**2+se(c_,c)**2)) for t_,c_ in EXP]
    w=[1/s_**2 for _,s_ in arr]; est=sum(wi*d for wi,(d,_) in zip(w,arr))/sum(w); s_=sqrt(1/sum(w))
    worst=min(worst,zp(est/s_))
check('exercise factorial: null on all 16 clocks, smallest pooled P = 0.09',
      worst>0.05 and abs(worst-0.092)<0.01, f'smallest pooled P = {worst:.3f}')
check('exercise factorial n: 30+30 per contrast, 120 total',
      all(nn(a)==30 for pair in EXP for a in pair))
print(f'==== rerun complete ====')
