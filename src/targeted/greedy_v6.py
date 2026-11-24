# version 6: additional analyses of the sequential selection (greedy marginal-benefit optimization)
#  A. Budget unit: select by 10% of lane length (lanes x length) or paved area (carriageway width x length) instead of road length
#  B. Dependence of the selected sections on the model: reselect with 20 models that change the exposure error, the exposure-response slope and the pavement effect, and with the Japanese curve;
#     loss of benefit (regret) when the sections selected with the central model are evaluated with each model, and overlap of the selected sections
#  C. Upper bound of the additive approximation: compare with the fractional knapsack (LP relaxation) whose value is the sum of stand-alone effects
# usage: greedy_v6.py <M dir> <noise data dir>
import os
import numpy as np, pandas as pd, scipy.sparse as sp, glob, sys, time
B=os.environ.get('F_ROOT','.')+'/'; MAT=sys.argv[1]; S=sys.argv[2]
SC=pd.read_parquet(B+'data/alloc/section_scenarios.parquet',columns=['id','len_m']); nsec=len(SC); LEN=SC.len_m.values.astype(float)
FL=np.load(MAT+'/section_fl.npy').astype(np.float64)
G=pd.read_parquet(f'{S}/geo/traffic_census_2021_converted.parquet',columns=['census','車線数','幅員構成／車道部幅員（ｍ）']).drop_duplicates('census')
G['census']=G.census.astype(str); G=G.set_index('census')
lanes=pd.to_numeric(G['車線数'],errors='coerce').reindex(SC.id.astype(str).values).values; carw=pd.to_numeric(G['幅員構成／車道部幅員（ｍ）'],errors='coerce').reindex(SC.id.astype(str).values).values
print('lanes missing',np.isnan(lanes).mean().round(4),'carw missing',np.isnan(carw).mean().round(4))
lanes=np.where(np.isnan(lanes)|(lanes<=0),2,lanes); carw=np.where(np.isnan(carw)|(carw<=0),lanes*3.25,carw)
COST={'len':LEN,'lane':LEN*lanes,'area':LEN*carw}
fW=lambda L:(19.4312-0.9336*L+0.0126*L**2)/100; fJ=lambda L:1/(1+np.exp(5.6431-0.0602*L))
def LLc(x): return np.clip(np.where(x>0,10*np.log10(np.maximum(x,1e-30)),40),40,None)
bc=np.log(2.13); bse=(np.log(2.48)-np.log(1.82))/(2*1.96)
lg=lambda L,b: 1/(1+np.exp(-(np.log(fW(40)/(1-fW(40)))+b*(L-40)/10)))
def make_exc(kind='who',b=bc):
    if kind=='J': return lambda E: fJ(LLc(E))-fJ(40)
    def f(E):
        L=LLc(E); base=fW(L)-fW(40); den=lg(L,bc)-fW(40); num=lg(L,b)-fW(40)
        return base*np.where(den>1e-9,num/np.maximum(den,1e-9),1.0)
    return f
P=[]
for f in sorted(glob.glob(MAT+'/M_*.npz')):
    pc=f[-6:-4]; z=np.load(f); M=sp.csr_matrix((z['data'],z['indices'],z['indptr']),shape=(len(z['indptr'])-1,nsec))
    P.append(dict(pc=pc,M=M,row=np.repeat(np.arange(M.shape[0],dtype=np.int32),np.diff(M.indptr)),E=z['E'],Eg=np.asarray(M.sum(axis=1)).ravel().astype(np.float64),
                  pw=z['pw'].astype(np.float64),kc=pd.read_parquet(f'{MAT}/kc_{pc}.parquet').KEY_CODE.values))
kc=np.concatenate([p['kc'] for p in P]); cat=pd.Categorical(kc); cats=pd.Series(cat.categories); code=cat.codes
def mp(s):
    v=cats.map(s); m=v.isna(); v[m]=cats[m].str[:9].map(s); return v.values
AG=pd.read_csv(B+'data/age_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0).share20plus
adi=pd.read_csv(B+'data/adi_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0).p_rent
s20=mp(AG); s20=np.where(np.isnan(s20),np.nanmedian(s20),s20); pw=np.concatenate([p['pw'] for p in P])
SA=pd.DataFrame({'KEY_CODE':cats.values,'w':np.bincount(code,pw,len(cats)),'p_rent':mp(adi)}); SA['city']=SA.KEY_CODE.str[:5]
X=SA.dropna(subset=['p_rent']).sort_values(['city','p_rent']); cw=X.groupby('city').w.cumsum()/X.groupby('city').w.transform('sum')
SA['q']=0; SA.loc[X.index,'q']=np.minimum((cw.values-1e-12)*5,4.999).astype(int)+1
q=SA.q.values[code]; w20=pw*s20[code]; city_codes,city=np.unique(cats.str[:5].values[code],return_inverse=True); pref=cats.str[:2].values[code].astype(int); o=0
for p in P:
    n=len(p['pw']); p['w20']=w20[o:o+n]; p['q']=q[o:o+n]; p['sl']=slice(o,o+n); o+=n
W1=w20[q==1].sum(); W5=w20[q==5].sum()
def ctx_central(): return dict(fg=[np.ones(len(p['pw'])) for p in P],E=[p['E'] for p in P],exc=make_exc(),red=1-10**(-0.3))
def ctx_draw(rng):
    gp=rng.normal(0,0.97,48); gc=rng.normal(0,1.37,len(city_codes)); mpf=rng.normal(0,2.87,48); b=rng.normal(bc,bse); pav=rng.triangular(1.5,3.0,4.5)
    fg=10**((gp[pref]+gc[city])/10); fm=10**(mpf[pref]/10)
    return dict(fg=[fg[p['sl']] for p in P],E=[p['Eg']*fg[p['sl']]+np.maximum(p['E']-p['Eg'],0)*fm[p['sl']] for p in P],exc=make_exc('who',b),red=1-10**(-pav/10),b=b,pav=pav)
def coef(meas,ctx): return np.full(nsec,ctx['red']) if meas=='pav' else 0.5*FL
def evaluate(sel,meas,ctx):
    x=coef(meas,ctx)*sel; t=t0=a1=a5=0
    for p,fg,E in zip(P,ctx['fg'],ctx['E']):
        h0=p['w20']*ctx['exc'](E); h=p['w20']*ctx['exc'](E-fg*(p['M']@x)); t+=h.sum(); t0+=h0.sum(); a1+=h[p['q']==1].sum(); a5+=h[p['q']==5].sum()
    return dict(red=(1-t/t0)*100,ratio=(a5/W5)/(a1/W1),avoided=t0-t)
def greedy(meas,gamma,budget,cost,ctx,step=0.01):
    c=coef(meas,ctx); C=COST[cost]; TOT=C.sum(); sel=np.zeros(nsec,bool); Ecur=[E.copy() for E in ctx['E']]; m0=None
    for r in range(int(round(budget/step))):
        mg=np.zeros(nsec)
        for p,Ec,fg in zip(P,Ecur,ctx['fg']):
            om=p['w20']*(1+gamma*np.maximum(p['q']-1,0)/4); M=p['M']; ex=ctx['exc'](Ec)
            v=om[p['row']]*(ex[p['row']]-ctx['exc'](Ec[p['row']]-c[M.indices]*M.data*fg[p['row']]))
            mg+=np.bincount(M.indices,v,nsec)
        if r==0: m0=mg.copy()
        eff=np.where(sel,-1,mg/np.maximum(C,1)); idx=np.argsort(-eff); k=int(np.searchsorted(np.cumsum(C[idx]),step*TOT))+1
        new=idx[:k]; new=new[eff[new]>0]; sel[new]=True; x=np.zeros(nsec); x[new]=c[new]
        for i,(p,fg) in enumerate(zip(P,ctx['fg'])): Ecur[i]=Ecur[i]-fg*(p['M']@x)
    return sel,m0
t00=time.time(); rows=[]; CEN=ctx_central()
# A. budget unit
for cost in ('len','lane','area'):
    for g in (0,4,16):
        s,_=greedy('pav',g,0.10,cost,CEN); e=evaluate(s.astype(float),'pav',CEN)
        rows.append(dict(part='A_cost',cost=cost,gamma=g,len_pct=LEN[s].sum()/LEN.sum()*100,lane_pct=COST['lane'][s].sum()/COST['lane'].sum()*100,area_pct=COST['area'][s].sum()/COST['area'].sum()*100,**e))
    print('A',cost,round(time.time()-t00),flush=True)
# all sections reference
ALL=np.ones(nsec); e=evaluate(ALL,'pav',CEN); rows.append(dict(part='A_cost',cost='all',gamma=np.nan,**e))
# C. LP bound with standalone values (gamma 0, length budget)
s0,m0=greedy('pav',0,0.10,'len',CEN); e0=evaluate(s0.astype(float),'pav',CEN)
r=m0/np.maximum(LEN,1); idx=np.argsort(-r); cl=np.cumsum(LEN[idx]); k=int(np.searchsorted(cl,0.10*LEN.sum())); frac=(0.10*LEN.sum()-(cl[k-1] if k>0 else 0))/LEN[idx[k]]
UB=m0[idx[:k]].sum()+frac*m0[idx[k]]
rows.append(dict(part='C_bound',cost='len',gamma=0,avoided=e0['avoided'],sum_solo_selected=m0[s0].sum(),LP_upper_additive=UB,gap_pct=(1-e0['avoided']/UB)*100))
print('C gap %',round((1-e0['avoided']/UB)*100,3),'additivity',round(e0['avoided']/m0[s0].sum(),4),flush=True)
# B. selection robustness
for meas in ('pav','hv'):
    sc,_=greedy(meas,0,0.10,'len',CEN) if meas=='hv' else (s0,None)
    rng=np.random.default_rng(606)
    for d in range(21):
        ctx=ctx_draw(rng) if d<20 else dict(fg=CEN['fg'],E=CEN['E'],exc=make_exc('J'),red=CEN['red'])
        sd,_=greedy(meas,0,0.10,'len',ctx); ec=evaluate(sc.astype(float),meas,ctx); ed=evaluate(sd.astype(float),meas,ctx)
        rows.append(dict(part='B_robust',measure=meas,draw=d if d<20 else 'Jcurve',red_central_sel=ec['red'],red_reopt=ed['red'],regret_pct=(1-ec['avoided']/ed['avoided'])*100,
                         overlap_len_pct=LEN[sc&sd].sum()/LEN[sc].sum()*100,ratio_central_sel=ec['ratio'],ratio_reopt=ed['ratio']))
        print('B',meas,d,round(rows[-1]['regret_pct'],3),round(rows[-1]['overlap_len_pct'],1),round(time.time()-t00),flush=True)
T=pd.DataFrame(rows); pd.set_option('display.width',250); print(T.round(3).to_string(index=False)); T.round(4).to_csv(B+'data/v6_greedy_checks.csv',index=False)
