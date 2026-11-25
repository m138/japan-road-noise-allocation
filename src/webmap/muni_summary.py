# Municipal summary of the main allocations, for the web map and the map in the README (not used in the paper).
#  For each municipality of the 2020 Population Census (wards of the designated cities are separate):
#   adults, excess high sleep disturbance (HSD) before mitigation, cases avoided by the efficient allocation (gamma = 0) and by the
#   equity-weighted allocation (gamma = 16) of low-noise pavement (uniform 3 dB, 10% of the national length of ordinary roads),
#   the same by within-municipality quintile of renter share (lowest and highest), length of ordinary roads and selected length.
#  A receptor belongs to the municipality of its small area. A section belongs to the municipality of the small area that
#  contains its representative point (alloc_prep.py), so the lengths are approximate where a section crosses a boundary.
# output: data/muni_allocation.csv
# usage: muni_summary.py <M dir> <noise data dir>
import os
import sys, numpy as np, pandas as pd, time
B=os.environ.get('F_ROOT','.')+'/'
exec(open(B+'src/targeted/greedy_v6.py',encoding='utf-8').read().split('t00=time.time()')[0])
CEN=ctx_central(); c=1-10**(-0.3); nc=len(city_codes); t00=time.time()
O=pd.read_parquet(B+'data/alloc/v5_opt_orders.parquet'); SK=pd.read_parquet(B+'data/alloc/section_scenarios.parquet',columns=['id','len_m','KEY_CODE'])
assert (O.id.values==SC.id.values).all() and (SK.id.values==SC.id.values).all()
SEL={g:(O['pav_order_g'+g].values>=0)&(O['pav_order_g'+g].values<10) for g in ('0','16')}
# receptors: cases by municipality and quintile (index city*6+q; q = 0 where the renter share is missing)
K=city*6+q; H={k:np.zeros(nc*6) for k in ('0','g0','g16')}
for p,fg,E in zip(P,CEN['fg'],CEN['E']):
    k=K[p['sl']]; H['0']+=np.bincount(k,p['w20']*CEN['exc'](E),nc*6)
    for g in ('0','16'): H['g'+g]+=np.bincount(k,p['w20']*CEN['exc'](E-fg*(p['M']@(c*SEL[g]))),nc*6)
    print(p['pc'],round(time.time()-t00),flush=True)
A=np.bincount(K,w20,nc*6).reshape(nc,6); H={k:v.reshape(nc,6) for k,v in H.items()}
T=pd.DataFrame({'municipality':city_codes,'population':np.bincount(city,pw,nc),'adults':A.sum(1),'excess_hsd':H['0'].sum(1),
                'avoided_efficient':(H['0']-H['g0']).sum(1),'avoided_equity_gamma16':(H['0']-H['g16']).sum(1)})
for j,nm in ((1,'q1'),(5,'q5')):
    T['adults_'+nm]=A[:,j]; T['excess_hsd_'+nm]=H['0'][:,j]; T['avoided_efficient_'+nm]=(H['0']-H['g0'])[:,j]; T['avoided_equity_gamma16_'+nm]=(H['0']-H['g16'])[:,j]
# sections: length by municipality
m=SK.KEY_CODE.str[:5]; ok=m.notna().values; L=pd.DataFrame({'municipality':m[ok].values,'road_m':LEN[ok],'selected_m_efficient':LEN[ok]*SEL['0'][ok],'selected_m_equity_gamma16':LEN[ok]*SEL['16'][ok]}).groupby('municipality').sum()
T=T.merge(L,left_on='municipality',right_index=True,how='left').fillna({'road_m':0,'selected_m_efficient':0,'selected_m_equity_gamma16':0})
t0=T.excess_hsd.sum(); print('excess HSD',round(t0,1),'reduction %',round(T.avoided_efficient.sum()/t0*100,3),round(T.avoided_equity_gamma16.sum()/t0*100,3))
print('length without municipality %',round(LEN[~ok].sum()/LEN.sum()*100,3),'selected % of all length',round(LEN[SEL['0']].sum()/LEN.sum()*100,3),round(LEN[SEL['16']].sum()/LEN.sum()*100,3))
T.round(4).to_csv(B+'data/muni_allocation.csv',index=False); print(len(T),'municipalities')
