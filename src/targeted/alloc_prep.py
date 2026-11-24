# Section-level allocation: census sections of ordinary roads are ranked by the priority of the small area that contains their representative point, and the sections treated up to 10/20/30% of length are decided
# output: data/alloc/section_scenarios.parquet (section id x treatment flag of each scenario), scenario_list.csv
import geopandas as gpd, pandas as pd, numpy as np, glob, sys, os
from pyproj import Geod
S=sys.argv[1]; B=os.environ.get('F_ROOT','.')+'/'; O=B+'data/alloc/'; os.makedirs(O,exist_ok=True)
sys.path.insert(0,B+'src/noise_scen'); import noise_lib as nl
K=nl.load_sections(f'{S}/geo/traffic_census_2021_converted.parquet')
st=(K.motorway=='1')|K.rt.isin(['1','2']); K=K[~st][['id','geometry']].dissolve(by='id').reset_index().to_crs(4326)
g=Geod(ellps='GRS80'); K['len_m']=[g.geometry_length(x) for x in K.geometry.values]
P=K[['id','len_m']].copy(); P=gpd.GeoDataFrame(P,geometry=K.geometry.representative_point().values,crs=4326)
A=pd.concat([gpd.read_file(f)[lambda d:d.HCODE==8101][['KEY_CODE','geometry']].to_crs(4326) for f in sorted(glob.glob(f'{S}/estat_all/shp/r2ka*.shp'))])
j=gpd.sjoin(P,A,how='left',predicate='within').drop_duplicates('id')
P['KEY_CODE']=j.set_index('id').KEY_CODE.reindex(P.id).values
print('sections',len(P),'km',round(P.len_m.sum()/1e3),'unassigned share',round(P.len_m[P.KEY_CODE.isna()].sum()/P.len_m.sum(),4),flush=True)
# small-area scores (same as targeted_pavement.py)
R=pd.read_parquet(B+'data/receptor/receptor_burden_v2.parquet',columns=['KEY_CODE','pop_w','share20plus','Lnight_2021','Ln_S3'])
adi=pd.read_csv(B+'data/adi_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0)
fW=lambda L:(19.4312-0.9336*L+0.0126*L**2)/100
LL=lambda x: np.clip(np.nan_to_num(np.asarray(x,'float64'),nan=40),40,None)
L0=LL(R.Lnight_2021.values); v3=R.Ln_S3.values.astype('float64'); L3=np.where(np.isnan(v3),L0,LL(v3))
w=R.pop_w.values.astype('float64'); w20=w*R.share20plus.fillna(R.share20plus.median()).values
R=pd.DataFrame({'KEY_CODE':R.KEY_CODE.values,'dH':w20*(fW(L0)-fW(L3)),'w':w,'wL':w*L0})
SA=R.groupby('KEY_CODE')[['dH','w','wL']].sum(); del R
SA['Lmean']=SA.wL/SA.w
sl=P.groupby('KEY_CODE').len_m.sum(); SA['road_m']=sl.reindex(SA.index).fillna(0)
pr=adi.p_rent; SA['p_rent']=SA.index.to_series().map(pr).fillna(SA.index.to_series().str[:9].map(pr))
SA['eff']=SA.dH/SA.road_m.replace(0,np.nan)
for c in ['eff','p_rent','Lmean']: SA['r_'+c]=SA[c].rank(pct=True)
P=P.merge(SA[['Lmean','eff','p_rent','r_eff','r_p_rent','r_Lmean']],left_on='KEY_CODE',right_index=True,how='left')
tot=P.len_m.sum(); cols={}
def sel(score):
    o=np.argsort(-np.nan_to_num(score,nan=-1e9),kind='stable'); c=np.cumsum(P.len_m.values[o]); out={}
    for b in (0.1,0.2,0.3):
        m=np.zeros(len(P),bool); m[o[c<=b*tot]]=True; out[b]=m
    return out
strat={'exposure':P.r_Lmean.values,'efficiency':P.r_eff.values,'rent':P.r_p_rent.values}
for lam in (0.25,0.5,0.75): strat[f'mix{lam}']=(1-lam)*P.r_eff.values+lam*P.r_p_rent.values
for i in range(20): strat[f'random{i:02d}']=np.random.default_rng(1000+i).random(len(P))
for k,s in strat.items():
    for b,m in sel(s).items(): cols[f'{k}_b{int(b*100)}']=m
cols['all']=np.ones(len(P),bool)
F=pd.DataFrame(cols); F.insert(0,'id',P.id.values); F.insert(1,'len_m',P.len_m.values); F.insert(2,'KEY_CODE',P.KEY_CODE.values)
F.to_parquet(O+'section_scenarios.parquet'); print('scenarios',len(cols)); print(F.drop(columns=['id','KEY_CODE']).mul(F.len_m,axis=0).sum().div(tot).round(3).head(12).to_string())
