# Addition of version 10: length of the sections that remain eligible with continuous treatment units (2 km, 5 km) and the reduction from treating all eligible sections (to check the denominators of 71% and 62%)
# usage: v10_units_denominator.py <M dir> <noise data dir>
import os
import sys, numpy as np, pandas as pd, time
B=os.environ.get('F_ROOT','.')+'/'
S7=open(B+'src/targeted/greedy_v7.py',encoding='utf-8').read()
head=S7.split("cU=np.full(nsec")[0].replace("__file__.replace('greedy_v7.py','greedy_v6.py')","B+'src/targeted/greedy_v6.py'")
exec(head)
blk=S7[S7.index("# --- 3. contiguity"):S7.index("KC5=")]
exec(blk)
cU=np.full(nsec,1-10**(-0.3))
def red(sel):
    x=cU*sel; t=t0=0
    for p,E in zip(P,CEN['E']):
        t0+=(p['w20']*CEN['exc'](E)).sum(); t+=(p['w20']*CEN['exc'](E-(p['M']@x))).sum()
    return (1-t/t0)*100
rows=[dict(name='all_sections',len_pct=100.0,n_sections=nsec,red=red(np.ones(nsec)))]
for Lmin in (2000,5000):
    u=blocks(Lmin); el=u>=0
    rows.append(dict(name=f'eligible_{Lmin//1000}km',len_pct=LEN[el].sum()/LEN.sum()*100,len_km=LEN[el].sum()/1e3,n_sections=int(el.sum()),n_units=int(u.max()+1),red=red(el.astype(float))))
    print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(B+'data/v10_units_denominator.csv',index=False); print(pd.DataFrame(rows))
