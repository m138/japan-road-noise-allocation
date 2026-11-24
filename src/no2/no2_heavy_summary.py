# Fall in deaths attributable to near-road NO2 from a 50% reduction in nighttime heavy-vehicle traffic (no shift; age 30 and over, RR 1.02 per 10 ug/m3, same equation as the main analysis)
# the near-road term is proportional to sqrt(road kernel) in the NO2 model, so in the scenario it is road x sqrt(1-dC/C)
import os
import pandas as pd, numpy as np, glob, sys
B=os.environ.get('F_ROOT','.')+'/'; D=sys.argv[1]
H=pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(D+'/no2h_*.parquet'))],ignore_index=True)
R=pd.read_parquet(B+'data/receptor/receptor_burden_v2.parquet',columns=['rid','KEY_CODE','pop_w','road_2022'])
R=R.merge(H,on='rid',how='left'); print('unmatched',R.C.isna().mean())
cats=pd.Series(pd.Categorical(R.KEY_CODE).categories); code=pd.Categorical(R.KEY_CODE,categories=cats).codes
def mp(s):
    v=cats.map(s); m=v.isna(); v[m]=cats[m].str[:9].map(s); return v.values[code]
A2=pd.read_csv(B+'data/age2_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0)
s30=mp(A2.share30plus); s30=np.where(np.isnan(s30),np.nanmedian(s30),s30); dr=mp(A2.dr_among30); dr=np.where(np.isnan(dr),np.nanmedian(dr),dr)
w30=R.pop_w.values*s30
paf=lambda u,rr:1-1/rr**(np.maximum(u,0)/10)
ur=R.road_2022.values*1.88; D0=(w30*dr*paf(ur,1.02)).sum()
yll=float(pd.read_csv(B+'data/yll_per_death_30plus_lifetable2022.csv').iloc[0,-1])
T=pd.read_csv(B+'data/v5_alloc_opt.csv'); hv=T[T.measure=='hv'].set_index(['method','budget']).HSD_k; base=T[T.method=='base'].HSD_k.iloc[0]
rows=[]
for c in [c for c in R.columns if c.startswith('rel_')]:
    k=c[4:]; rel=np.nan_to_num(R[c].values.astype(float)); Dk=(w30*dr*paf(ur*np.sqrt(np.clip(1-rel,0,1)),1.02)).sum()
    meth,b=k.rsplit('_b',1); dH=(base-hv.get((meth,int(b)),np.nan))*1e3
    rows.append(dict(scenario=k,deaths_avoided=D0-Dk,pct_of_road_no2_deaths=(D0-Dk)/D0*100,YLL_avoided=(D0-Dk)*yll,HSD_avoided=dH,YLD_avoided_dw0010=dH*0.010))
O=pd.DataFrame(rows); O['YLL_over_YLD']=O.YLL_avoided/O.YLD_avoided_dw0010
print('road NO2 deaths base',round(D0),'YLL per death',yll); print(O.round(2).to_string(index=False))
O.round(3).to_csv(B+'data/v5_no2_heavy.csv',index=False)
