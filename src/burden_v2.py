# v2: burden from noise (receptors and scenarios), NO2 and tropical nights (receptors, urban correction of the tropical-nights paper), and its distribution by socioeconomic indicators
import os
import pandas as pd, numpy as np, pyarrow.parquet as pq, glob, sys
B=os.environ.get('F_ROOT','.')+'/'; NAT=sys.argv[1]
R=pq.read_table(os.environ.get('F_RECEPTORS','receptors_no2_noise_2015_2023.parquet'),
  columns=['rid','KEY_CODE','lon','lat','pop_w','Lnight_2021','NO2_2022','road_2022']).to_pandas()
H=pq.read_table(B+'data/receptor/receptor_tn25_2018_2023.parquet').to_pandas()
R=R.merge(H[['rid','tn25_mean_2018_2023','tn25_2023']],on='rid',how='left')
# scenarios: match by KEY_CODE + rounded coordinates
cols_s=['lon','lat','Ln','Ln_S1','Ln_S2','Ln_S3','Ln_S4']
S=pd.concat([pd.read_parquet(f,columns=cols_s) for f in sorted(glob.glob(NAT+'/rec_*.parquet'))],ignore_index=True)
def kk(d): return (d.lon.values*1e5).round().astype('int64')*100000000+(d.lat.values*1e5).round().astype('int64')
S['k']=kk(S); S=S.drop(columns=['lon','lat']).drop_duplicates('k')
for c in S.columns:
    if c!='k': S[c]=S[c].astype('float32')
R['k']=kk(R)
R=R.merge(S,on='k',how='left',indicator=True); del S
ok=(R._merge=='both').values; R=R.drop(columns='_merge')
print('scenario match share (pop)',round(R.loc[ok,'pop_w'].sum()/R.pop_w.sum(),4))
base=np.clip(np.nan_to_num(R.Ln.values,nan=40),40,None)
print('baseline check |Ln-Lnight_2021|>0.05 dB share',round(float(np.mean(np.abs(base[ok]-R.Lnight_2021.values[ok])>0.05)),4))
adi=pd.read_csv(B+'data/adi_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0)
AG=pd.read_csv(B+'data/age_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0)
def mp(s):
    v=R.KEY_CODE.map(s); mm=v.isna(); v[mm]=R.KEY_CODE[mm].str[:9].map(s); return v
for c in ['ADI_raw','p_rent','p_unemp','p_lone_parent']: R[c]=mp(adi[c])
for c in ['share20plus','death_rate30plus']:
    R[c]=mp(AG[c]); R[c]=R[c].fillna(R[c].median())
R['city']=R.KEY_CODE.str[:5]
w=R.pop_w.values.astype('float64'); wa=w*R.share20plus.values
def hsd(L): L=np.clip(np.nan_to_num(L,nan=40),40,None); return np.where(L>=45,np.clip(19.4312-0.9336*L+0.0126*L**2,0,None)/100,0)
def hsdj(L): L=np.clip(np.nan_to_num(L,nan=40),40,None); return np.where(L>=45,1/(1+np.exp(5.6431-0.0602*L)),0)
L0=R.Lnight_2021.values
R['HSD']=wa*hsd(L0); R['HSDj']=wa*hsdj(L0)
for s in ['S1','S2','S3','S4']:
    Ls=np.where(ok,np.clip(np.nan_to_num(R['Ln_'+s].values,nan=40),40,None),L0); R['HSD_'+s]=wa*hsd(Ls)
dr=R.death_rate30plus.values
ug=R.road_2022.values*1.88
R['D_road']=w*dr*(1-1/1.02**(ug/10))
R['D_road_E50']=w*dr*(1-1/1.02**(ug*np.sqrt(0.5)/10))   # road NOx emissions halved (the near-road term is proportional to sqrt(emissions))
# tropical nights -> insufficient rest (paper E OR per 10 nights), baseline prevalence p0
OR=1.0105; p0=0.35
TN=np.nan_to_num(R.tn25_mean_2018_2023.values)
po=p0*OR**(TN/10)/(1-p0+p0*OR**(TN/10)); R['REST_heat']=wa*(po-p0)
R['wTN']=w*TN
cols=['pop_w','HSD','HSDj','HSD_S1','HSD_S2','HSD_S3','HSD_S4','D_road','D_road_E50','REST_heat','wTN']
tot=R[cols].sum(); print('\nNational totals (thousand persons / deaths):'); print((tot/1e3).round(1).to_string())
def quint(v,within=False):
    v=pd.Series(v)
    if within: rk=v.groupby(R.city.values).rank(method='first',pct=True).values
    else:
        o=np.argsort(v.values,kind='stable'); cw=np.empty(len(v)); cw[o]=np.cumsum(w[o])/w[~np.isnan(v.values)].sum() if False else np.cumsum(w[o])/w.sum(); rk=cw
    q=np.minimum((np.asarray(rk)*5).astype(int),4)+1; q=np.where(v.isna(),0,q); return q
rows=[]
for var in ['ADI_raw','p_rent']:
    for scope in ['national','within']:
        R['q']=quint(R[var].values,scope=='within')
        g=R[R.q>0].groupby('q')[cols].sum()
        for q,s in g.iterrows():
            rows.append(dict(var=var,scope=scope,q=q,pop_M=s.pop_w/1e6,HSD_per1000=s.HSD/s.pop_w*1e3,HSD_S4_per1000=s.HSD_S4/s.pop_w*1e3,
              D_road_per100k=s.D_road/s.pop_w*1e5,REST_heat_per1000=s.REST_heat/s.pop_w*1e3,TN_mean=s.wTN/s.pop_w))
T=pd.DataFrame(rows).round(2); print(T.to_string()); T.to_csv(B+'data/v2_by_quintile.csv',index=False)
R[['rid','KEY_CODE','pop_w','share20plus','Lnight_2021','Ln_S1','Ln_S2','Ln_S3','Ln_S4','tn25_mean_2018_2023','road_2022','HSD','HSD_S4','D_road','REST_heat']].to_parquet(B+'data/receptor/receptor_burden_v2.parquet')
