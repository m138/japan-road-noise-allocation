# version 6: robustness of the selected sections to estimation errors that are independent between sections (pavement, gamma=0, 10% of length)
# The nighttime energy of each section is multiplied by an independent error 10^(e/10), e~N(0, SD), to make a 'true' state;
#  (a) the allocation reselected in that state and (b) the allocation selected with the central model are evaluated in that state to obtain the regret (loss of benefit, %)
# the preparation is shared with greedy_v6.py
# usage: greedy_v6_secnoise.py <M dir> <noise data dir>
import os
import sys, numpy as np, pandas as pd, time, scipy.sparse as sp
src=open(__file__.replace('greedy_v6_secnoise.py','greedy_v6.py'),encoding='utf-8').read().split('t00=time.time()')[0]
exec(src)
B=os.environ.get('F_ROOT','.')+'/'; CEN=ctx_central(); t0=time.time(); rows=[]
P0=[dict(p) for p in P]
sc,_=greedy('pav',0,0.10,'len',CEN)
rng=np.random.default_rng(707)
for sd in (2.0,3.0):
    for d in range(10):
        f=10**(rng.normal(0,sd,nsec)/10)
        P[:]=[]
        for p in P0:
            q=dict(p); M=p['M'].copy(); M.data=(M.data*f[M.indices]).astype(np.float32); q['M']=M
            q['E']=np.maximum(p['E']-p['Eg'],0)+np.asarray(M.sum(axis=1)).ravel(); P.append(q)
        ctx=dict(fg=CEN['fg'],E=[p['E'] for p in P],exc=CEN['exc'],red=CEN['red'])
        sd_sel,_=greedy('pav',0,0.10,'len',ctx); ec=evaluate(sc.astype(float),'pav',ctx); ed=evaluate(sd_sel.astype(float),'pav',ctx)
        rows.append(dict(sd_db=sd,draw=d,red_central_sel=ec['red'],red_reopt=ed['red'],regret_pct=(1-ec['avoided']/ed['avoided'])*100,
                         overlap_len_pct=LEN[sc&sd_sel].sum()/LEN[sc].sum()*100,ratio_central_sel=ec['ratio'],ratio_reopt=ed['ratio']))
        print(rows[-1],round(time.time()-t0),flush=True)
pd.DataFrame(rows).round(4).to_csv(B+'data/v6_greedy_secnoise.csv',index=False)
