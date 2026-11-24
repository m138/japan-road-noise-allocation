# Addition to version 8: how the 10-year average of ASJ is calculated, and the treatment of negative reductions
#  (a) energy mean (asj_avg10 of greedy_v8_asj.py): the attenuation coefficients of sections over 10 years are averaged and sleep disturbance is calculated once
#  (b) mean of the sleep benefit: avoided cases are calculated with the reduction of each year y=0..9 and averaged over the 10 years
#  (c) in (b), reductions for light and heavy vehicles are limited to 0 dB from below (sensitivity analysis with an assumption added to the original equations)
#  Allocations evaluated: all sections, efficient allocation selected with the uniform 3 dB (gamma=0, 10% of length), gamma=16, efficient allocation reselected with the energy mean
# usage: greedy_v8b_asj_years.py <M dir> <noise data dir>
import sys, numpy as np, pandas as pd, time
src=open(__file__.replace('greedy_v8b_asj_years.py','greedy_v7.py'),encoding='utf-8').read().split('# --- 1. pavement effect')[0]
src=src.replace("src=open(__file__.replace('greedy_v7.py','greedy_v6.py')","src=open(__file__.replace('greedy_v8b_asj_years.py','greedy_v6.py')")
exec(src)
def dL(y,floor):
    l,h=5.7-7.3*np.log10(1+y),3.9-3.6*np.log10(1+y)
    return (max(l,0),max(h,0)) if floor else (l,h)
def cy(y,floor=False): l,h=dL(y,floor); return 1-((1-FL)*10**(-l/10)+FL*10**(-h/10))
cavg=1-np.mean([1-cy(y) for y in range(10)],axis=0)
savg=greedy_g(cavg,omR(0),0.10,'len')
SEL={'all':ALL,'uniform_sel_g0':s0,'uniform_sel_g16':s16,'avg10_reopt_g0':savg}
out=[]
for nm,sel in SEL.items():
    ea=ev(sel.astype(float),cavg); out.append(dict(sel=nm,method='energy_avg',red=ea['red'],avoided=ea['avoided'],rent_ratio=ea['rent_ratio'],rent_gap=ea['rent_gap']))
    for floor in (False,True):
        E=[ev(sel.astype(float),cy(y,floor)) for y in range(10)]
        out.append(dict(sel=nm,method='health_avg'+('_floor0' if floor else ''),red=np.mean([e['red'] for e in E]),avoided=np.mean([e['avoided'] for e in E]),
                        rent_ratio=np.mean([e['rent_ratio'] for e in E]),rent_gap=np.mean([e['rent_gap'] for e in E])))
    print(nm,[(o['method'],round(o['red'],3),round(o['rent_ratio'],3)) for o in out if o['sel']==nm],round(time.time()-t00),flush=True)
T=pd.DataFrame(out)
for m in T.method.unique():
    a=T[(T.method==m)&(T.sel=='all')].avoided.iloc[0]; T.loc[T.method==m,'share_all']=T.loc[T.method==m,'avoided']/a*100
print(T.round(3).to_string(index=False)); T.round(5).to_csv(B+'data/v8b_asj_years.csv',index=False)
