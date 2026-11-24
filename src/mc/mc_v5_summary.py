# Summary of the Monte Carlo runs (version 5): interval of each allocation, 95% interval of the difference between allocations within the same run, stability between the first and second 500 runs, Japanese exposure-response curve (draw=-1)
import os
import pandas as pd, numpy as np, glob
B=os.environ.get('F_ROOT','.')+'/'
D=pd.concat([pd.read_csv(f) for f in sorted(glob.glob(B+'data/v5_mc_draws_*.csv'))],ignore_index=True)
C=D[D.draw==0].iloc[0]; J=D[D.draw==-1].iloc[0]; X=D[D.draw>0].reset_index(drop=True); print('draws',len(X))
cols=[c for c in X.columns if c.startswith(('red_','ratio_','diff_'))]+['H_base','ratio_base','diff_base']
q=lambda s:s.quantile([.025,.5,.975]).values
rows=[dict(stat=c,central=C[c],jcurve=J[c],lo=q(X[c])[0],med=q(X[c])[1],hi=q(X[c])[2]) for c in cols]
PAIRS=[('pav','opt_g0_b10','old_efficiency_b10'),('pav','opt_g0_b10','solo_b10'),('pav','opt_g1_b10','opt_g0_b10'),('pav','opt_g4_b10','opt_g0_b10'),('pav','opt_g16_b10','opt_g0_b10'),
       ('pav','old_rent_b10','opt_g0_b10'),('pav','old_mix0.25_b10','old_efficiency_b10'),('pav','opt_g0_b10','old_random00_b10'),
       ('hv','opt_g0_b10','old_efficiency_b10'),('hv','opt_g1_b10','opt_g0_b10'),('hv','opt_g4_b10','opt_g0_b10'),('hv','opt_g16_b10','opt_g0_b10'),('hv','old_rent_b10','opt_g0_b10')]
for t,a,b in PAIRS:
    for m in ('red','ratio'):
        ca,cb=f'{m}_{t}_{a}',f'{m}_{t}_{b}'
        if ca not in X or cb not in X: continue
        d=X[ca]-X[cb]; h1=q(d[:len(d)//2]); h2=q(d[len(d)//2:])
        rows.append(dict(stat=f'DIFF {m} {t} {a} − {b}',central=C[ca]-C[cb],jcurve=J[ca]-J[cb],lo=q(d)[0],med=q(d)[1],hi=q(d)[2],p_pos=(d>0).mean(),half1=f'{h1[0]:.3f}..{h1[2]:.3f}',half2=f'{h2[0]:.3f}..{h2[2]:.3f}'))
    if t=='pav' or True:
        ca,cb=f'red_{t}_{a}',f'red_{t}_{b}'
        if ca in X and cb in X: rows.append(dict(stat=f'RATIO_OF_RED {t} {a} / {b}',central=C[ca]/C[cb],jcurve=J[ca]/J[cb],lo=q(X[ca]/X[cb])[0],med=q(X[ca]/X[cb])[1],hi=q(X[ca]/X[cb])[2]))
for t in ('pav','hv'):
    s=X[f'red_{t}_opt_g0_b10']/X[f'red_{t}_all_b100']; rows.append(dict(stat=f'SHARE {t} opt_g0_b10 / all',central=C[f'red_{t}_opt_g0_b10']/C[f'red_{t}_all_b100'],jcurve=J[f'red_{t}_opt_g0_b10']/J[f'red_{t}_all_b100'],lo=q(s)[0],med=q(s)[1],hi=q(s)[2]))
O=pd.DataFrame(rows); pd.set_option('display.width',250); pd.set_option('display.max_rows',300); pd.set_option('display.max_colwidth',60)
print(O.round(3).to_string(index=False)); O.round(4).to_csv(B+'data/v5_mc_summary.csv',index=False)
