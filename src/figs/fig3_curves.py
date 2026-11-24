# For Fig. 3: reduction in sleep disturbance against the treated share of road length and lane length
#  Efficient allocation (gamma=0) and equity-weighted allocation (gamma=16): sequential selection in steps of 1% of length up to 100%; reduction and disparity after each round
#  Renter-share-only (sections in descending order of the renter share of the small area of the representative point, as in version 4) and random (seeds 1000-1002): evaluated in steps of 5%
#  Pavement effect of the main analysis (uniform 3 dB for all vehicle types)
# usage: fig3_curves.py <M dir> <noise data dir>
import os
import sys, numpy as np, pandas as pd, time
exec(open(os.environ.get('F_ROOT','.')+'/src/targeted/greedy_v6.py').read().split('t00=time.time()')[0])
CEN=ctx_central(); c=coef('pav',CEN); t00=time.time()
def summ(Es):
    t=a1=a5=0
    for p,E in zip(P,Es):
        h=p['w20']*CEN['exc'](E); t+=h.sum(); a1+=h[p['q']==1].sum(); a5+=h[p['q']==5].sum()
    return t,a1,a5
T0,A10,A50=summ(CEN['E'])
def row(method,cost,gamma,sel,t,a1,a5):
    return dict(method=method,cost=cost,gamma=gamma,len_pct=LEN[sel].sum()/LEN.sum()*100,lane_pct=COST['lane'][sel].sum()/COST['lane'].sum()*100,
                red=(1-t/T0)*100,ratio=(a5/W5)/(a1/W1),gap=(a5/W5-a1/W1)*1e3)
rows=[row('base','-',np.nan,np.zeros(nsec,bool),T0,A10,A50)]
def evsel(sel):
    x=c*sel; return summ([E-(p['M']@x) for p,E in zip(P,CEN['E'])])
for cost in ('len','lane'):
    C=COST[cost]; TOT=C.sum()
    for g in (0,16):
        sel=np.zeros(nsec,bool); Ecur=[E.copy() for E in CEN['E']]
        for r in range(100):
            mg=np.zeros(nsec)
            for p,Ec in zip(P,Ecur):
                om=p['w20']*(1+g*np.maximum(p['q']-1,0)/4); M=p['M']; ex=CEN['exc'](Ec)
                mg+=np.bincount(M.indices,om[p['row']]*(ex[p['row']]-CEN['exc'](Ec[p['row']]-c[M.indices]*M.data)),nsec)
            eff=np.where(sel,-np.inf,mg/np.maximum(C,1)); idx=np.argsort(-eff,kind='stable'); k=int(np.searchsorted(np.cumsum(C[idx]),0.01*TOT))+1
            new=idx[:k]; new=new[eff[new]>0]
            if len(new)==0: break
            sel[new]=True; x=np.zeros(nsec); x[new]=c[new]
            for i,p in enumerate(P): Ecur[i]=Ecur[i]-(p['M']@x)
            rows.append(row('greedy',cost,g,sel,*summ(Ecur)))
        print(cost,g,r,round(rows[-1]['red'],2),round(time.time()-t00),flush=True)
        pd.DataFrame(rows).to_csv(B+'data/fig3_curves.csv',index=False)
    # renter-share-only and random (fixed order, evaluated in steps of 5% of the budget)
    SS=pd.read_parquet(B+'data/alloc/section_scenarios.parquet',columns=['id','KEY_CODE'])
    adi=pd.read_csv(B+'data/adi_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0).p_rent
    pr=SS.KEY_CODE.map(adi).fillna(SS.KEY_CODE.str[:9].map(adi)).values
    orders={'rent_only':np.argsort(-np.nan_to_num(pd.Series(pr).rank(pct=True).values,nan=-1e9),kind='stable')}
    for i in range(3): orders[f'random{i}']=np.argsort(-np.random.default_rng(1000+i).random(nsec),kind='stable')
    for name,o in orders.items():
        cl=np.cumsum(C[o])
        for b in range(5,101,5):
            sel=np.zeros(nsec,bool); sel[o[cl<=b/100*TOT]]=True
            rows.append(row(name,cost,np.nan,sel,*evsel(sel)))
        print(cost,name,round(time.time()-t00),flush=True)
        pd.DataFrame(rows).to_csv(B+'data/fig3_curves.csv',index=False)
sel=np.ones(nsec,bool); rows.append(row('all','-',np.nan,sel,*evsel(sel)))
pd.DataFrame(rows).to_csv(B+'data/fig3_curves.csv',index=False); print('done',round(time.time()-t00))
