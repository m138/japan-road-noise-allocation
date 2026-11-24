# Share of the variance of each indicator over the 1,000 Monte Carlo runs: exposure-response slope b, pavement reduction pav, and the rest (exposure error and interactions)
# Method: R2 of the indicator regressed on b and b^2, R2 on pav and pav^2, and R2 on both (+ b*pav);
#       1 - R2 of the joint model is the unexplained share (exposure error, its interactions with the other factors, and variation not captured by the quadratic form). b and pav are drawn independently.
import os
import pandas as pd, numpy as np
B=os.environ.get('F_ROOT','.')+'/'
d=pd.concat([pd.read_csv(B+'data/v5_mc_draws_a.csv'),pd.read_csv(B+'data/v5_mc_draws_b.csv')]); d=d[d.draw>0]
g0,g16,al='pav_opt_g0_b10','pav_opt_g16_b10','pav_all_b100'
M={'全国の睡眠妨害の人数':d.H_base,'対策前の格差の比':d.ratio_base,
   '削減率（効率配分）':d['red_'+g0],'全区間に対する割合':d['red_'+g0]/d['red_'+al],
   'γ16の相対的な損失':1-d['red_'+g16]/d['red_'+g0],'格差の比の変化（γ16−γ0）':d['ratio_'+g16]-d['ratio_'+g0],
   '大型車：γ16の相対的な損失':1-d['red_hv_opt_g16_b10']/d['red_hv_opt_g0_b10'],
   '大型車：格差の比の変化':d['ratio_hv_opt_g16_b10']-d['ratio_hv_opt_g0_b10']}
def r2(y,X):
    X=np.column_stack([np.ones(len(y))]+X); b=np.linalg.lstsq(X,y,rcond=None)[0]; e=y-X@b; return 1-e.var()/y.var()
b=d.b.values; p=d.pav.values; rows=[]
for k,y in M.items():
    y=y.values.astype(float)
    rb=r2(y,[b,b**2]); rp=r2(y,[p,p**2]); rall=r2(y,[b,b**2,p,p**2,b*p])
    rows.append(dict(metric=k,p2_5=np.percentile(y,2.5),p97_5=np.percentile(y,97.5),slope=rb,pavement=rp,both=rall,unexplained=1-rall))
T=pd.DataFrame(rows); T.to_csv(B+'data/v8_mc_variance_share.csv',index=False)
print(T.round(3).to_string()); print('corr(b,pav)',np.corrcoef(b,p)[0,1])
