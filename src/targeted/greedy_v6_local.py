# version 6: check how close the greedy solution is (pavement, gamma=0, 10% of length, central model)
#  1. step of selection reduced from 1% to 0.25%
#  2. local search by swaps: selected sections with a small loss per length when removed are swapped with unselected sections with a large benefit per length when added,
#     for the same length, and the swap is accepted if the objective improves (steps 0.5% -> 0.1% -> 0.02%)
# The LP bound from the sum of stand-alone effects is not an upper bound because the effects of selected sections are synergistic (additivity 1.03), so this check is used.
# The preparation (matrix, weights, functions) is shared with greedy_v6.py (the part before execution is read)
# usage: greedy_v6_local.py <M dir> <noise data dir>
import os
import sys, numpy as np, pandas as pd, time
src=open(__file__.replace('greedy_v6_local.py','greedy_v6.py'),encoding='utf-8').read().split('t00=time.time()')[0]
exec(src)
B=os.environ.get('F_ROOT','.')+'/'; CEN=ctx_central(); c=coef('pav',CEN); t0=time.time(); rows=[]
def state(sel):
    x=c*sel; Ec=[p['E']-p['M']@x for p in P]; return Ec
def total(Ec): return sum((p['w20']*CEN['exc'](E)).sum() for p,E in zip(P,Ec))
H0=total([p['E'] for p in P])
def gains(sel,Ec):
    add=np.zeros(nsec); rem=np.zeros(nsec)
    for p,E in zip(P,Ec):
        M=p['M']; r=p['row']; d=c[M.indices]*M.data; ex=CEN['exc'](E)[r]; w=p['w20'][r]
        s=sel[M.indices]
        add+=np.bincount(M.indices,np.where(s,0,w*(ex-CEN['exc'](E[r]-d))),nsec)
        rem+=np.bincount(M.indices,np.where(s,w*(CEN['exc'](E[r]+d)-ex),0),nsec)
    return add,rem
for step in (0.01,0.0025):
    s,_=greedy('pav',0,0.10,'len',CEN,step=step); H=total(state(s)); rows.append(dict(method=f'greedy_step{step}',red=(1-H/H0)*100,len_pct=LEN[s].sum()/LEN.sum()*100))
    print(rows[-1],round(time.time()-t0),flush=True)
    if step==0.01: sel=s.copy()
Ec=state(sel); H=total(Ec); TOT=LEN.sum(); budget=0.10*TOT; it=0
for chunk in (0.005,0.001,0.0002):
    while it<40:
        it+=1; add,rem=gains(sel,Ec)
        ri=np.where(sel)[0]; ri=ri[np.argsort(rem[ri]/np.maximum(LEN[ri],1))]; kr=int(np.searchsorted(np.cumsum(LEN[ri]),chunk*TOT))+1; out=ri[:kr]
        free=budget-(LEN[sel].sum()-LEN[out].sum())
        ai=np.where(~sel)[0]; ai=ai[np.argsort(-add[ai]/np.maximum(LEN[ai],1))]; ka=int(np.searchsorted(np.cumsum(LEN[ai]),free)); inn=ai[:ka]
        new=sel.copy(); new[out]=False; new[inn]=True; Ec2=state(new); H2=total(Ec2)
        rows.append(dict(method=f'swap_chunk{chunk}',iter=it,red=(1-H2/H0)*100,accepted=H2<H,len_pct=LEN[new].sum()/TOT*100,n_out=len(out),n_in=len(inn)))
        print(rows[-1],round(time.time()-t0),flush=True)
        if H2<H: sel,Ec,H=new,Ec2,H2
        else: break
rows.append(dict(method='local_search_final',red=(1-H/H0)*100,len_pct=LEN[sel].sum()/TOT*100))
print(rows[-1],flush=True)
pd.DataFrame(rows).round(5).to_csv(B+'data/v6_greedy_local.csv',index=False)
