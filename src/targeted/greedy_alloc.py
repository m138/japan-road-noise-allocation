# version 5: sequential selection (greedy) by the marginal number of avoided cases of each section (excess high sleep disturbance in adults, counterfactual 40 dB)
#  - pavement (source of the selected sections -3 dB) and heavy vehicles (nighttime heavy vehicles on the selected sections -50%, no shift) are optimized separately
#  - each round selects 1% of the length of ordinary roads; the marginal benefit is recalculated after the effect of the selected sections is applied (overlap of effects is considered)
#  - equity weight: maximize avoided cases with weight 1+gamma(q-1)/4 for adults in within-municipality quintile q of renter share (gamma=0 gives the largest benefit)
#  - the ranking by the marginal benefit of round 0 (effect of each section alone) is also kept as the section-level stand-alone efficiency
# usage: greedy_alloc.py <M dir> <out dir for mcrec>
import numpy as np, pandas as pd, scipy.sparse as sp, glob, sys, os, time, json
B=os.environ.get('F_ROOT','.')+'/'; MAT=sys.argv[1]; OUT=sys.argv[2]; os.makedirs(OUT,exist_ok=True)
SC=pd.read_parquet(B+'data/alloc/section_scenarios.parquet'); nsec=len(SC); LEN=SC.len_m.values.astype(float); TOT=LEN.sum()
FL=np.load(MAT+'/section_fl.npy').astype(np.float64); RED=1-10**(-0.3)
fW=lambda L:(19.4312-0.9336*L+0.0126*L**2)/100
def LLc(x): return np.clip(np.where(x>0,10*np.log10(np.maximum(x,1e-30)),40),40,None)
def exc(E): return fW(LLc(E))-fW(40)
P=[]
for f in sorted(glob.glob(MAT+'/M_*.npz')):
    pc=f[-6:-4]; z=np.load(f)
    M=sp.csr_matrix((z['data'],z['indices'],z['indptr']),shape=(len(z['indptr'])-1,nsec))
    P.append(dict(pc=pc,M=M,row=np.repeat(np.arange(M.shape[0],dtype=np.int32),np.diff(M.indptr)),E=z['E'],pw=z['pw'].astype(np.float64),
                  kc=pd.read_parquet(f'{MAT}/kc_{pc}.parquet').KEY_CODE.values))
print('prefectures',len(P),'nnz',sum(p['M'].nnz for p in P),flush=True)
# adults and within-city rent quintile (same as alloc_summary.py)
kc=np.concatenate([p['kc'] for p in P]); cat=pd.Categorical(kc); cats=pd.Series(cat.categories); code=cat.codes
def mp(s):
    v=cats.map(s); m=v.isna(); v[m]=cats[m].str[:9].map(s); return v.values
AG=pd.read_csv(B+'data/age_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0).share20plus
adi=pd.read_csv(B+'data/adi_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0).p_rent
s20=mp(AG); s20=np.where(np.isnan(s20),np.nanmedian(s20),s20)
pw=np.concatenate([p['pw'] for p in P])
SA=pd.DataFrame({'KEY_CODE':cats.values,'w':np.bincount(code,pw,len(cats)),'p_rent':mp(adi)}); SA['city']=SA.KEY_CODE.str[:5]
X=SA.dropna(subset=['p_rent']).sort_values(['city','p_rent']); cw=X.groupby('city').w.cumsum()/X.groupby('city').w.transform('sum')
SA['q']=0; SA.loc[X.index,'q']=np.minimum((cw.values-1e-12)*5,4.999).astype(int)+1
q=SA.q.values[code]; w20=pw*s20[code]; o=0
for p in P:
    n=len(p['pw']); p['w20']=w20[o:o+n]; p['q']=q[o:o+n]; o+=n
def coef(measure): return np.full(nsec,RED) if measure=='pav' else 0.5*FL
def evaluate(sel,measure):
    x=coef(measure)*sel; t=0; a1=a5=0
    for p in P:
        h=p['w20']*exc(p['E']-p['M']@x); t+=h.sum(); a1+=h[p['q']==1].sum(); a5+=h[p['q']==5].sum()
    return t,a1,a5
W1=sum(p['w20'][p['q']==1].sum() for p in P); W5=sum(p['w20'][p['q']==5].sum() for p in P)
def greedy(measure,gamma,rounds):
    c=coef(measure); sel=np.zeros(nsec,bool); order=np.full(nsec,-1); Ecur=[p['E'].copy() for p in P]; m0=None
    for r in range(rounds):
        mg=np.zeros(nsec)
        for p,Ec in zip(P,Ecur):
            om=p['w20']*(1+gamma*np.maximum(p['q']-1,0)/4)
            M=p['M']; ex=exc(Ec)
            v=om[p['row']]*(ex[p['row']]-exc(Ec[p['row']]-c[M.indices]*M.data))
            mg+=np.bincount(M.indices,v,nsec)
        if r==0: m0=mg.copy()
        eff=np.where(sel,-1,mg/np.maximum(LEN,1)); idx=np.argsort(-eff); cl=np.cumsum(LEN[idx]); k=int(np.searchsorted(cl,0.01*TOT))+1
        new=idx[:k]; new=new[eff[new]>0]; sel[new]=True; order[new]=r
        x=np.zeros(nsec); x[new]=c[new]
        for i,p in enumerate(P): Ecur[i]=Ecur[i]-p['M']@x
    return order,m0
res={}; t00=time.time()
RUNS=[('pav',0,30),('pav',0.5,20),('pav',1,20),('pav',2,20),('pav',4,20),('pav',8,20),('pav',16,20),('hv',0,30),('hv',1,20),('hv',2,20),('hv',4,20),('hv',8,20),('hv',16,20)]
for meas,g,R in RUNS:
    order,m0=greedy(meas,g,R); res[(meas,g)]=order
    if g==0: res[(meas,'solo')]=m0
    print(meas,g,'done',round(time.time()-t00),'s',flush=True)
# selections
SEL={}
def topk(score,b):
    idx=np.argsort(-score); cl=np.cumsum(LEN[idx]); k=int(np.searchsorted(cl,b/100*TOT))+1; s=np.zeros(nsec,bool); s[idx[:k]]=True; return s
for meas in ('pav','hv'):
    for (m,g),order in res.items():
        if m!=meas or g=='solo': continue
        for b in (10,20,30):
            if order.max()+1>=b: SEL[(meas,f'opt_g{g}',b)]=(order>=0)&(order<b)
    for b in (10,20,30): SEL[(meas,'solo',b)]=topk(res[(meas,'solo')]/np.maximum(LEN,1),b)
    for old in ('efficiency','mix0.25','mix0.5','rent','exposure','random00'):
        for b in (10,20,30): SEL[(meas,'old_'+old,b)]=SC[f'{old}_b{b}'].values.astype(bool)
    SEL[(meas,'all',100)]=SC['all'].values.astype(bool)
rows=[]; tb,_,_=evaluate(np.zeros(nsec),'pav')
_,b1,b5=evaluate(np.zeros(nsec),'pav'); rows.append(dict(measure='none',method='base',gamma=np.nan,budget=0,HSD_k=tb/1e3,reduction_pct=0,Q1=b1/W1*1e3,Q5=b5/W5*1e3))
for (meas,meth,b),s in SEL.items():
    t,a1,a5=evaluate(s.astype(float),meas)
    rows.append(dict(measure=meas,method=meth,budget=b,len_pct=LEN[s].sum()/TOT*100,HSD_k=t/1e3,reduction_pct=(1-t/tb)*100,Q1=a1/W1*1e3,Q5=a5/W5*1e3))
T=pd.DataFrame(rows); T['diff']=T.Q5-T.Q1; T['ratio']=T.Q5/T.Q1
pd.set_option('display.width',220); print(T.round(3).to_string(index=False))
T.round(4).to_csv(B+'data/v5_alloc_opt.csv',index=False)
S=pd.DataFrame({'id':SC.id.values,'len_m':LEN})
for (meas,g),order in res.items(): S[f'{meas}_{"solo" if g=="solo" else "order_g"+str(g)}']=order
S.to_parquet(B+'data/alloc/v5_opt_orders.parquet')
# mcrec for Monte Carlo
PK=[('pav','opt_g0',10),('pav','opt_g0.5',10),('pav','opt_g1',10),('pav','opt_g2',10),('pav','opt_g4',10),('pav','opt_g16',10),('pav','solo',10),('pav','old_efficiency',10),('pav','old_mix0.25',10),('pav','old_rent',10),('pav','old_random00',10),('pav','opt_g0',20),('pav','all',100)]
HK=[('hv','opt_g0',10),('hv','opt_g1',10),('hv','opt_g4',10),('hv','opt_g16',10),('hv','solo',10),('hv','old_efficiency',10),('hv','old_rent',10),('hv','old_random00',10),('hv','all',100)]
for p in P:
    RC=pd.DataFrame({'KEY_CODE':p['kc'],'w':p['pw'].astype(np.float32),'E':p['E'].astype(np.float32),'Egen':np.asarray(p['M'].sum(axis=1)).ravel().astype(np.float32)})
    for meas,meth,b in PK: RC[f'Ep_{meth}_b{b}']=(p['M']@SEL[(meas,meth,b)].astype(float)).astype(np.float32)
    for meas,meth,b in HK: RC[f'Eh_{meth}_b{b}']=(p['M']@(SEL[(meas,meth,b)]*FL)).astype(np.float32)
    RC.to_parquet(f"{OUT}/mcrec_{p['pc']}.parquet")
print('mcrec written',round(time.time()-t00),'s')
