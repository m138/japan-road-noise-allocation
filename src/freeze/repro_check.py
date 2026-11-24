# Reproducibility check before freezing version 8:
#  (1) rerun the section selection of the main analysis (uniform 3 dB, 10% of road length, gamma=0 and 16) and compare with the stored orders (data/alloc/v5_opt_orders.parquet)
#  (2) check that the selection does not change with a stable sort (does not depend on ties)
#  (3) summary values of population weights, quintiles and adults (to see that they are the same with another matrix directory)
# usage: repro_check.py <M dir> <noise data dir> <tag>
import os
import sys, numpy as np, pandas as pd, hashlib
src=open(os.environ.get('F_ROOT','.')+'/src/targeted/greedy_v6.py').read().split('t00=time.time()')[0]
exec(src)
TAG=sys.argv[3]
def greedy_k(gamma,kind):
    c=coef('pav',CEN); C=COST['len']; TOT=C.sum(); order=np.full(nsec,-1); Ecur=[E.copy() for E in CEN['E']]; ties=0
    for r in range(10):
        mg=np.zeros(nsec)
        for p,Ec in zip(P,Ecur):
            om=p['w20']*(1+gamma*np.maximum(p['q']-1,0)/4); M=p['M']; ex=CEN['exc'](Ec)
            mg+=np.bincount(M.indices,om[p['row']]*(ex[p['row']]-CEN['exc'](Ec[p['row']]-c[M.indices]*M.data)),nsec)
        eff=np.where(order>=0,-1,mg/np.maximum(C,1)); idx=np.argsort(-eff,kind=kind); k=int(np.searchsorted(np.cumsum(C[idx]),0.01*TOT))+1
        ties+=int(eff[idx[k-1]]==eff[idx[k]])
        new=idx[:k]; new=new[eff[new]>0]; order[new]=r; x=np.zeros(nsec); x[new]=c[new]
        for i,p in enumerate(P): Ecur[i]=Ecur[i]-(p['M']@x)
    return order,ties
CEN=ctx_central()
O=pd.read_parquet(B+'data/alloc/v5_opt_orders.parquet').set_index('id')
rows=[]
for g in (0,16):
    for kind in ('quicksort','stable'):
        o,t=greedy_k(g,kind); sel=(o>=0)
        e=evaluate(sel.astype(float),'pav',CEN)
        saved=O[f'pav_order_g{g}'].reindex(SC.id.values).fillna(-1).values
        ss=(saved>=0)&(saved<10)
        rows.append(dict(tag=TAG,gamma=g,sort=kind,red=e['red'],ratio=e['ratio'],n_sel=int(sel.sum()),
            same_sel_as_saved=bool((sel==ss).all()),same_order_as_saved=bool((np.where(sel,o,-1)==np.where(ss,saved,-1)).all()),
            boundary_ties=t,sel_hash=hashlib.sha256(np.where(sel,o,-1).astype(np.int16).tobytes()).hexdigest()[:16]))
        print(rows[-1],flush=True)
summ=dict(tag=TAG,n_receptor=len(pw),pop=pw.sum(),adult=w20.sum(),W1=W1,W5=W5,
          q_share=[round(float(pw[q==k].sum()/pw.sum()),5) for k in range(6)],
          E_hash=hashlib.sha256(np.concatenate([p['E'] for p in P]).tobytes()).hexdigest()[:16],
          pw_hash=hashlib.sha256(pw.tobytes()).hexdigest()[:16],
          M_hash=hashlib.sha256(b''.join(p['M'].data.tobytes()+p['M'].indices.tobytes() for p in P)).hexdigest()[:16])
print(summ)
pd.DataFrame(rows).to_csv(B+f'data/freeze_repro_{TAG}.csv',index=False)
pd.Series(summ).to_csv(B+f'data/freeze_inputs_{TAG}.csv')
