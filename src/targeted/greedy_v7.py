# version 7: additional checks of the sequential maximization of the marginal benefit (same preparation as greedy_v6.py)
#  1. Pavement effect that depends on road characteristics (speed, share of heavy vehicles): the pavement lowers only tire-road (rolling) noise by X dB and leaves propulsion noise unchanged.
#     The share of rolling noise by vehicle category and speed comes from CNOSSOS-EU (Directive (EU) 2015/996, coefficients as amended in 2021, A-weighted octave bands 63 Hz-8 kHz),
#     equations for passenger cars (category 1) and heavy vehicles (category 3). Speed is the posted speed limit used in the noise model (ordinary roads 10-60 km/h, 40 if missing, lower limit 20 for CNOSSOS).
#     Reduction coefficient of a section c = (1-fl)*rs1*(1-10^(-X/10)) + fl*rs3*(1-10^(-X/10)) (fl: share of heavy vehicles in the nighttime energy).
#     V1: X=3 dB (as it is). V2: the section-level reductions (dB) of V1 multiplied by a constant, keeping their relative shape, so that treating all sections has the same effect as the uniform 3 dB (to see only the difference in allocation).
#     (rolling noise alone does not reach the overall effect of a uniform 3 dB on ordinary roads with many heavy vehicles at night, so a constant factor was used instead of adjusting X)
#  2. Definition of equity: absolute gap of renter share (per 1,000 adults, Q5-Q1), weighting by the areal deprivation index (ADI, within-municipality quintiles), and an allocation that aims only to narrow the absolute gap
#  3. Continuity of treatment: consecutive sections of the same route are joined into treatment units of 2 km or 5 km or more and selected by unit (shorter units are excluded). The degree of scatter is also summarized
# usage: greedy_v7.py <M dir> <noise data dir>
import os
import sys, numpy as np, pandas as pd, time
src=open(__file__.replace('greedy_v7.py','greedy_v6.py'),encoding='utf-8').read().split('t00=time.time()')[0]
exec(src)
import geopandas as gpd
B=os.environ.get('F_ROOT','.')+'/'; CEN=ctx_central(); t00=time.time(); rows=[]
# --- equity indices
ADI=pd.read_csv(B+'data/adi_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0).ADI_raw
SA['adi']=mp(ADI); X=SA.dropna(subset=['adi']).sort_values(['city','adi']); cw=X.groupby('city').w.cumsum()/X.groupby('city').w.transform('sum')
SA['qa']=0; SA.loc[X.index,'qa']=np.minimum((cw.values-1e-12)*5,4.999).astype(int)+1; qa=SA.qa.values[code]
for p in P: p['qa']=qa[p['sl']]
WA1=w20[qa==1].sum(); WA5=w20[qa==5].sum()
def ev(sel,c,ctx=CEN):
    x=c*sel; t=t0=0; a={}
    for p,fg,E in zip(P,ctx['fg'],ctx['E']):
        h0=p['w20']*ctx['exc'](E); h=p['w20']*ctx['exc'](E-fg*(p['M']@x)); t+=h.sum(); t0+=h0.sum()
        for nm,qq in (('r',p['q']),('a',p['qa'])):
            for k in (1,5):
                a[f'{nm}{k}']=a.get(f'{nm}{k}',0)+h[qq==k].sum(); a[f'{nm}{k}b']=a.get(f'{nm}{k}b',0)+h0[qq==k].sum()
    r1,r5,a1,a5=a['r1']/W1*1e3,a['r5']/W5*1e3,a['a1']/WA1*1e3,a['a5']/WA5*1e3
    r1b,r5b,a1b,a5b=a['r1b']/W1*1e3,a['r5b']/W5*1e3,a['a1b']/WA1*1e3,a['a5b']/WA5*1e3
    return dict(red=(1-t/t0)*100,avoided=t0-t,rent_ratio=r5/r1,rent_gap=r5-r1,rent_gap0=r5b-r1b,rent_ratio0=r5b/r1b,adi_ratio=a5/a1,adi_gap=a5-a1,adi_gap0=a5b-a1b,adi_ratio0=a5b/a1b,len_pct=LEN[sel>0].sum()/LEN.sum()*100)
def greedy_g(c,omf,budget,cost,ctx=CEN,units=None,step=0.01):
    C=COST[cost]; TOT=C.sum(); sel=np.zeros(nsec,bool); Ecur=[E.copy() for E in ctx['E']]
    if units is not None: ok=units>=0; nu=units.max()+1; UC=np.bincount(units[ok],C[ok],nu); usel=np.zeros(nu,bool)
    for r in range(int(round(budget/step))):
        mg=np.zeros(nsec)
        for p,Ec,fg in zip(P,Ecur,ctx['fg']):
            om=omf(p); M=p['M']; ex=ctx['exc'](Ec)
            mg+=np.bincount(M.indices,om[p['row']]*(ex[p['row']]-ctx['exc'](Ec[p['row']]-c[M.indices]*M.data*fg[p['row']])),nsec)
        if units is None:
            eff=np.where(sel,-np.inf,mg/np.maximum(C,1)); idx=np.argsort(-eff); k=int(np.searchsorted(np.cumsum(C[idx]),step*TOT))+1; new=idx[:k]; new=new[eff[new]>0]
        else:
            um=np.bincount(units[ok],mg[ok],nu); eff=np.where(usel,-np.inf,um/np.maximum(UC,1)); idx=np.argsort(-eff); k=int(np.searchsorted(np.cumsum(UC[idx]),step*TOT))+1
            nn=idx[:k]; nn=nn[eff[nn]>0]; usel[nn]=True; new=np.where(ok&np.isin(units,nn))[0]; new=new[~sel[new]]
        if len(new)==0: break
        sel[new]=True; x=np.zeros(nsec); x[new]=c[new]
        for i,(p,fg) in enumerate(zip(P,ctx['fg'])): Ecur[i]=Ecur[i]-fg*(p['M']@x)
    return sel
omR=lambda g: (lambda p: p['w20']*(1+g*np.maximum(p['q']-1,0)/4))
omA=lambda g: (lambda p: p['w20']*(1+g*np.maximum(p['qa']-1,0)/4))
omGap=lambda p: p['w20']*((p['q']==5)/W5-(p['q']==1)/W1)
def rec(part,name,sel,c,**kw):
    e=ev(sel.astype(float),c); rows.append(dict(part=part,name=name,**kw,**e)); print(part,name,{k:round(v,3) for k,v in e.items() if k in('red','rent_ratio','rent_gap','adi_ratio','adi_gap','len_pct')},round(time.time()-t00),flush=True); return e
cU=np.full(nsec,1-10**(-0.3)); ALL=np.ones(nsec,bool)
s0=greedy_g(cU,omR(0),0.10,'len'); s16=greedy_g(cU,omR(16),0.10,'len'); s4=greedy_g(cU,omR(4),0.10,'len')
# --- 1. pavement effect by road characteristics
Aw=np.array([-26.2,-16.1,-8.6,-3.2,0,1.2,1.0,-1.1])
CN={'1':dict(ar=[83.1,89.2,87.7,93.1,100.1,96.7,86.8,76.2],br=[30.0,41.5,38.9,25.7,32.5,37.2,39.0,40.0],ap=[97.9,92.5,90.7,87.2,84.7,88.0,84.4,77.1],bp=[-1.3,7.2,7.7,8.0,8.0,8.0,8.0,8.0]),
    '3':dict(ar=[91.7,96.2,98.2,104.9,105.1,98.5,91.1,85.6],br=[30.0,33.5,31.3,25.4,31.8,37.1,38.6,40.6],ap=[108.8,104.2,103.5,102.9,102.6,98.5,93.8,87.5],bp=[0.0,3.0,4.6,5.0,5.0,5.0,5.0,5.0])}
def rshare(cat,v):
    k={x:np.array(y) for x,y in CN[cat].items()}; v=np.maximum(np.asarray(v,float),20)[:,None]
    R=(10**((k['ar']+k['br']*np.log10(v/70)+Aw)/10)).sum(1); Pp=(10**((k['ap']+k['bp']*(v-70)/70+Aw)/10)).sum(1); return R/(R+Pp)
Gs=pd.read_parquet(f'{S}/geo/traffic_census_2021_converted.parquet',columns=['census','指定最高速度（ｋｍ／ｈ）']).drop_duplicates('census'); Gs['census']=Gs.census.astype(str)
vl=pd.to_numeric(Gs.set_index('census')['指定最高速度（ｋｍ／ｈ）'],errors='coerce').reindex(SC.id.astype(str).values).values
V=np.clip(np.where(np.isnan(vl)|(vl<=0),40,vl),10,60); RS1=rshare('1',V); RS3=rshare('3',V)
print('rolling share light',np.round(np.percentile(RS1,[5,50,95]),3),'heavy',np.round(np.percentile(RS3,[5,50,95]),3),'fl median',np.round(np.median(FL),3),flush=True)
def cvar(Xdb): return ((1-FL)*RS1+FL*RS3)*(1-10**(-Xdb/10))
dB3=-10*np.log10(1-cvar(3.0))   # section reduction (dB) when only rolling noise is cut by 3 dB
def cscale(k): return 1-10**(-k*dB3/10)
eU=ev(ALL.astype(float),cU)['red']; lo,hi=1.0,6.0
for _ in range(14):
    mid=(lo+hi)/2; (lo,hi)=(mid,hi) if ev(ALL.astype(float),cscale(mid))['red']<eU else (lo,mid)
ks=(lo+hi)/2; print('scale k',round(ks,3),'max possible all-section effect (rolling removed entirely)',round(ev(ALL.astype(float),cvar(60.0))['red'],2),flush=True)
for scen,c in (('uniform3',cU),('rolling3',cvar(3.0)),(f'rolling_pattern_scaled{ks:.2f}',cscale(ks))):
    eff_db=-10*np.log10(1-c); print(scen,'section reduction dB pct',np.round(np.percentile(eff_db,[5,25,50,75,95]),2),flush=True)
    ea=rec('1_pav',scen+'_all',ALL,c,scen=scen)
    for g in (0,16):
        sg=greedy_g(c,omR(g),0.10,'len'); e=rec('1_pav',f'{scen}_reopt_g{g}_len',sg,c,scen=scen,gamma=g,share_all=None)
        rows[-1]['share_all']=e['red']/ea['red']*100
        if g==0:
            sref=sg; eu=rec('1_pav',f'{scen}_uniformsel_g0_len',s0,c,scen=scen,gamma=0); rows[-1]['regret_pct']=(1-eu['avoided']/e['avoided'])*100; rows[-1]['overlap_len_pct']=LEN[s0&sg].sum()/LEN[s0].sum()*100; rows[-1]['share_all']=eu['red']/ea['red']*100
    sl=greedy_g(c,omR(0),0.10,'lane'); e=rec('1_pav',f'{scen}_reopt_g0_lane',sl,c,scen=scen,gamma=0); rows[-1]['share_all']=e['red']/ea['red']*100
# --- 2. equity definitions (uniform 3 dB)
rec('2_eq','eff_g0',s0,cU); rec('2_eq','rent_g4',s4,cU); rec('2_eq','rent_g16',s16,cU)
SCa=pd.read_parquet(B+'data/alloc/section_scenarios.parquet',columns=['rent_b10','random00_b10'])
rec('2_eq','rent_only_v4',SCa.rent_b10.values,cU); rec('2_eq','random',SCa.random00_b10.values,cU); rec('2_eq','all',ALL,cU)
for g in (4,16): rec('2_eq',f'adi_g{g}',greedy_g(cU,omA(g),0.10,'len'),cU)
rec('2_eq','absgap_only',greedy_g(cU,omGap,0.10,'len'),cU)
# --- 3. contiguity / construction units
GG=gpd.read_parquet(f'{S}/geo/traffic_census_2021_converted.parquet',columns=['census','geometry']); GG['census']=GG.census.astype(str)
GG=GG[GG.census.isin(set(SC.id.astype(str)))].dissolve(by='census')
geo=GG.geometry.reindex(SC.id.astype(str).values).values
ids=SC.id.astype(str).values; route=np.array([s[:7] for s in ids]); order=np.lexsort((ids,route))
chain=np.full(nsec,-1); cid=-1; prev=None
for i in order:
    adj=prev is not None and route[i]==route[prev] and geo[i] is not None and geo[prev] is not None and geo[i].distance(geo[prev])<0.0005
    if not adj: cid+=1
    chain[i]=cid; prev=i
def blocks(Lmin):
    u=np.full(nsec,-1); bid=-1
    for c_,grp in pd.Series(order).groupby(chain[order]):
        idx=grp.values; L=LEN[idx]
        if L.sum()<Lmin: continue
        cur=[]; acc=0; out=[]
        for j,l in zip(idx,L):
            cur.append(j); acc+=l
            if acc>=Lmin: out.append(cur); cur=[]; acc=0
        if cur: out[-1].extend(cur)
        for b in out: bid+=1; u[b]=bid
    return u
KC5=pd.read_parquet(B+'data/alloc/section_scenarios.parquet',columns=['KEY_CODE']).KEY_CODE.astype(str).str[:5].values
def runs_stats(sel,name):
    s=pd.DataFrame({'chain':chain[sel],'len':LEN[sel],'i':np.where(sel)[0]}).sort_values('i')
    # contiguous selected runs within a chain (consecutive in order)
    pos=np.empty(nsec,int); pos[order]=np.arange(nsec); s['pos']=pos[s.i.values]; s=s.sort_values('pos'); s['run']=((s.chain.diff()!=0)|(s.pos.diff()!=1)).cumsum()
    rl=s.groupby('run').len.sum()
    muni=pd.Series(KC5[sel]).nunique()
    d=dict(part='3_stats',name=name,n_sections=int(sel.sum()),n_runs=int(len(rl)),median_run_km=rl.median()/1e3,share_len_runs_ge2km=rl[rl>=2000].sum()/rl.sum()*100,share_len_runs_ge5km=rl[rl>=5000].sum()/rl.sum()*100,n_muni=muni)
    rows.append(d); print(d,flush=True)
print('chains',cid+1,flush=True)
runs_stats(s0,'free_g0_len')
for Lmin in (2000,5000):
    u=blocks(Lmin); print('Lmin',Lmin,'eligible len %',round(LEN[u>=0].sum()/LEN.sum()*100,1),'units',u.max()+1,flush=True)
    for g,cost in ((0,'len'),(16,'len'),(0,'area')):
        sb=greedy_g(cU,omR(g),0.10,cost,units=u); rec('3_unit',f'units{Lmin//1000}km_g{g}_{cost}',sb,cU,Lmin=Lmin,gamma=g,cost=cost,eligible_len_pct=LEN[u>=0].sum()/LEN.sum()*100)
        if g==0 and cost=='len': runs_stats(sb,f'units{Lmin//1000}km_g0_len')
sa=greedy_g(cU,omR(0),0.10,'area'); rec('3_unit','free_g0_area',sa,cU)
T=pd.DataFrame(rows); T.round(5).to_csv(B+'data/v7_checks.csv',index=False); print('done',round(time.time()-t00),flush=True)
