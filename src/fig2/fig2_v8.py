# Fig. 2 (version 8): the disparity uses the same definition as the allocation analysis
#   Small areas are ranked within each municipality from the lowest renter share (or ADI) and
#   divided into five equal parts of the cumulative population (sum of the population weights of receptors); same rule as greedy_v6.py.
#   The earlier version (burden_v3.py) ranked receptors. Both are written out to check the difference.
# (the preparation part of burden_v3.py is used as it is)
# v3:
#  - high sleep disturbance from the continuous exposure-response curve (counterfactual 40 dB, main). The number at 45 dB or higher (descriptive) and the excess over 45 dB are also given
#  - denominators of rates follow the target ages (HSD: 20 and over; NO2 deaths: 30 and over; tropical nights: 40-74); for NO2 deaths a rate that removes age structure is also given
#  - NO2 deaths: near-road term (main), total NO2 above the WHO guideline, sensitivity with RR 1.04 (HEI 2022), NOx halved
#  - tropical nights: ages 40-74, baseline proportion replaced by the observed values of the NDB Open Data (fiscal years 2018-2023, by area)
#  - YLD: WHO 2024 disability weight (sleep disturbance attributed to noise, 0.010) as main, 0.07 (2011) for comparison
#  - stratified by housing type (detached houses are half or more of main households, or not)
import os
import pandas as pd, numpy as np, pyarrow.parquet as pq, glob
B=os.environ.get('F_ROOT','.')+'/'
R=pd.read_parquet(B+'data/receptor/receptor_burden_v2.parquet',columns=['rid','KEY_CODE','pop_w','share20plus','Lnight_2021','Ln_S1','Ln_S2','Ln_S3','Ln_S4','tn25_mean_2018_2023','road_2022'])
N=pq.read_table(os.environ.get('F_RECEPTORS','receptors_no2_noise_2015_2023.parquet'),columns=['rid','NO2_2022']).to_pandas()
R=R.merge(N,on='rid',how='left'); del N
R['KEY_CODE']=R.KEY_CODE.astype('category'); _cats=pd.Series(R.KEY_CODE.cat.categories); _codes=R.KEY_CODE.cat.codes.values
R['city']=pd.Categorical(_cats.str[:5].values[_codes])
def mp(s):
    v=_cats.map(s); mm=v.isna(); v[mm]=_cats[mm].str[:9].map(s); return v.values[_codes]
adi=pd.read_csv(B+'data/adi_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0)
A2=pd.read_csv(B+'data/age2_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0)
nat30=float(pd.read_csv(B+'data/age2_national.csv',index_col=0).iloc[0,0])
for c in ['ADI_raw','p_rent']: R[c]=mp(adi[c])
for c in ['share30plus','share40_74','dr_among30','dr_na_among30']: R[c]=mp(A2[c]); R[c]=R[c].fillna(R[c].median())
# housing type
fs=[pd.read_csv(f,encoding='cp932',skiprows=5,dtype=str) for f in sorted(glob.glob(B+'census2020_small_area/t08_1/*.csv'))]
d=pd.concat(fs,ignore_index=True); d.columns=[c.strip() for c in d.columns]; d=d[d['町丁字コード']!='-']
d['KEY_CODE']=d['市区町村コード'].str.zfill(5)+d['町丁字コード']; d=d.set_index('KEY_CODE')
num=lambda s: pd.to_numeric(s.str.replace('-','0'),errors='coerce').astype(float)
R['p_detached']=mp((num(d['一戸建'])/num(d['うち主世帯'])).replace([np.inf,-np.inf],np.nan)); del d,fs
import gc; gc.collect()
# area baseline for insufficient rest
ma=pd.read_csv(B+'data/muni_to_area335.csv',dtype=str).set_index('municipality_code').area
pa=pd.read_csv(B+'data/ndb_insufficient_rest_area_2018_2023.csv',dtype={'area_code':str}).set_index('area_code').p_insuff
R['p_obs']=pd.Series(R.city.cat.categories).map(ma).map(pa).values[R.city.cat.codes.values]
print('heat baseline unmatched pop share',round(R.pop_w[R.p_obs.isna()].sum()/R.pop_w.sum(),4)); R['p_obs']=R.p_obs.fillna(0.3496)

w=R.pop_w.values.astype('float64'); w20=w*R.share20plus.values; w30=w*R.share30plus.values; w4074=w*R.share40_74.values
def fW(L): return (19.4312-0.9336*L+0.0126*L**2)/100
def fJ(L): return 1/(1+np.exp(5.6431-0.0602*L))
def LL(x): return np.clip(np.nan_to_num(np.asarray(x,dtype='float64'),nan=40),40,None)
L0=LL(R.Lnight_2021.values)
Ls={'base':L0}
for s in ['S1','S2','S3','S4']:
    v=R['Ln_'+s].values.astype('float64'); Ls[s]=np.where(np.isnan(v),L0,LL(v))
Ls['n2307']=L0+1.1   # sensitivity: about +1.1 dB on ordinary roads when the night is defined as 23:00-07:00 (WHO)
for k,L in Ls.items():
    R['prev45_'+k]=w20*np.where(L>=45,fW(L),0)          # definition of v2 (value of the curve at 45 dB or higher, 0 below)
    R['att40_'+k]=w20*(fW(L)-fW(40))                     # main: counterfactual 40 dB, continuous
    R['att45_'+k]=w20*np.clip(fW(L)-fW(45),0,None)*(L>45) # excess over 45 dB
    R['att40J_'+k]=w20*(fJ(L)-fJ(40))
# NO2
ur=R.road_2022.values*1.88; ut=R.NO2_2022.values*1.88
paf=lambda u,rr: 1-1/rr**(np.clip(u,0,None)/10)
dr=R.dr_among30.values
R['D_road']=w30*dr*paf(ur,1.02); R['D_road_rr104']=w30*dr*paf(ur,1.04); R['D_road_rr104_na']=w30*R.dr_na_among30.values*paf(ur,1.04); R['D_road_E50']=w30*dr*paf(ur*np.sqrt(0.5),1.02)
R['D_road_std']=w30*nat30*paf(ur,1.02); R['D_totWHO']=w30*dr*paf(ut-10,1.02)
# heat: 40-74, observed area prevalence; counterfactual = 0 tropical nights
OR=1.0105; TN=np.nan_to_num(R.tn25_mean_2018_2023.values); po=R.p_obs.values
o0=po/(1-po)/OR**(TN/10); p0=o0/(1+o0); R['REST_heat']=w4074*(po-p0)
R['w']=w; R['w20']=w20; R['w30']=w30; R['w4074']=w4074; R['wTN']=w*TN
hs=[c for c in R.columns if c.startswith(('prev45_','att40_','att45_','att40J_'))]
cols=['w','w20','w30','w4074','wTN']+hs+['D_road','D_road_rr104','D_road_rr104_na','D_road_E50','D_road_std','D_totWHO','REST_heat']

tot=R[cols].sum()
codes=R.KEY_CODE.cat.codes.values; ncat=len(_cats)
def sa_quint(var):
    SA=pd.DataFrame({'KEY_CODE':_cats.values,'w':np.bincount(codes,w,ncat),'v':pd.Series(R[var].values).groupby(codes).first().reindex(range(ncat)).values})
    SA['city']=SA.KEY_CODE.str[:5]
    X=SA.dropna(subset=['v']).sort_values(['city','v'],kind='stable'); cw=X.groupby('city').w.cumsum()/X.groupby('city').w.transform('sum')
    SA['q']=0; SA.loc[X.index,'q']=np.minimum((cw.values-1e-12)*5,4.999).astype(int)+1
    return SA.q.values[codes]
def rec_quint(v):
    v=pd.Series(v); ok=~v.isna().values
    rk=v.groupby(R.city.cat.codes.values).rank(method='first',pct=True).values
    q=np.minimum((np.asarray(rk)*5).astype(int),4)+1; return np.where(ok,q,0)
R['strat']=np.where(R.p_detached>=0.5,'detached_major',np.where(R.p_detached<0.5,'apartment_major','na'))
rows=[]
for var in ['p_rent','ADI_raw']:
    for unit in ['small_area','receptor']:
        R['q']=sa_quint(var) if unit=='small_area' else rec_quint(R[var].values)
        G=R[cols].groupby([R.q.values,R.strat.values]).sum()
        for st in ['all','detached_major','apartment_major']:
            g=G.groupby(level=0).sum() if st=='all' else G.xs(st,level=1)
            g=g[g.index>0]
            for q,s in g.iterrows():
                rows.append(dict(var=var,unit=unit,stratum=st,q=q,pop_M=s.w/1e6,adult_M=s.w20/1e6,
                  HSD_per1000adult=s.att40_base/s.w20*1e3,HSD_J_per1000adult=s.att40J_base/s.w20*1e3,
                  HSD_S3_per1000adult=s.att40_S3/s.w20*1e3,HSD_S4_per1000adult=s.att40_S4/s.w20*1e3,
                  NO2D_per100k30=s.D_road/s.w30*1e5,NO2D_std_per100k30=s.D_road_std/s.w30*1e5,NO2D_rr104_per100k30=s.D_road_rr104/s.w30*1e5,
                  REST_per1000_4074=s.REST_heat/s.w4074*1e3,TN_mean=s.wTN/s.w,
                  HSD_thousand=s.att40_base/1e3,NO2D=s.D_road))
        # quintiles are cut by the population weights of small areas, so each quintile should hold about the same population
T=pd.DataFrame(rows); T.round(4).to_csv(B+'data/fig2_v8_by_quintile.csv',index=False)
S=[]
for (var,unit,st),X in T.groupby(['var','unit','stratum']):
    X=X.set_index('q')
    for c in ['HSD_per1000adult','HSD_J_per1000adult','HSD_S3_per1000adult','HSD_S4_per1000adult','NO2D_per100k30','NO2D_std_per100k30','NO2D_rr104_per100k30','REST_per1000_4074','TN_mean']:
        S.append(dict(var=var,unit=unit,stratum=st,measure=c,Q1=X.loc[1,c],Q5=X.loc[5,c],diff=X.loc[5,c]-X.loc[1,c],ratio=X.loc[5,c]/X.loc[1,c]))
S=pd.DataFrame(S); S.round(3).to_csv(B+'data/fig2_v8_gaps.csv',index=False)
pd.set_option('display.width',250); print(S[S.stratum=='all'].round(3).to_string())
print(T[(T.stratum=='all')&(T['var']=='p_rent')][['unit','q','pop_M','adult_M']].round(3).to_string())
