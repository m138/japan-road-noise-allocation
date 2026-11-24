"""Files of the data deposit: results by small area and by road section.

small_area_burden_2020.csv : health burden before mitigation and within-municipality quintiles for each 2020 census small area
section_allocation.csv     : selection round and stand-alone benefit of each census section of ordinary roads

The definitions are those of the main analysis (src/fig2/fig2_v8.py and src/targeted/greedy_alloc.py).
usage: build_release_data.py <output dir>
"""
import os, sys, glob, numpy as np, pandas as pd
B=os.environ.get('F_ROOT','.')+'/'; OUT=sys.argv[1]; os.makedirs(OUT,exist_ok=True)
R=pd.read_parquet(B+'data/receptor/receptor_burden_v2.parquet',columns=['KEY_CODE','pop_w','share20plus','Lnight_2021','road_2022'])
R['KEY_CODE']=R.KEY_CODE.astype('category'); cats=pd.Series(R.KEY_CODE.cat.categories); code=R.KEY_CODE.cat.codes.values; n=len(cats)
def mp(s):   # value of the small area; the 9-digit town value is used when the 11-digit code is not found
    v=cats.map(s); m=v.isna(); v[m]=cats[m].str[:9].map(s); return v.values.astype(float)
def rec(v):  # receptor-level values, missing values filled with the median over receptors (as in fig2_v8.py)
    x=v[code]; return np.where(np.isnan(x),np.nanmedian(x),x)
adi=pd.read_csv(B+'data/adi_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0)
A2=pd.read_csv(B+'data/age2_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0)
nat30=float(pd.read_csv(B+'data/age2_national.csv',index_col=0).iloc[0,0])
p_rent=mp(adi.p_rent); ADI=mp(adi.ADI_raw)
# housing type (2020 Population Census, small-area tabulation, Table 8-1): detached houses / main households
fs=[pd.read_csv(f,encoding='cp932',skiprows=5,dtype=str) for f in sorted(glob.glob(B+'census2020_small_area/t08_1/*.csv'))]
d=pd.concat(fs,ignore_index=True); d.columns=[c.strip() for c in d.columns]; d=d[d['町丁字コード']!='-']
d['KEY_CODE']=d['市区町村コード'].str.zfill(5)+d['町丁字コード']; d=d.set_index('KEY_CODE')
num=lambda s: pd.to_numeric(s.str.replace('-','0'),errors='coerce').astype(float)
p_det=mp((num(d['一戸建'])/num(d['うち主世帯'])).replace([np.inf,-np.inf],np.nan)); del d,fs
# receptor-level burden
w=R.pop_w.values.astype('float64'); w20=np.nan_to_num(w*R.share20plus.values); w30=w*rec(mp(A2.share30plus)); dr=rec(mp(A2.dr_among30))
fW=lambda L:(19.4312-0.9336*L+0.0126*L**2)/100            # WHO curve, share highly sleep disturbed
L=np.clip(np.nan_to_num(R.Lnight_2021.values.astype('float64'),nan=40),40,None)
hsd=w20*(fW(L)-fW(40))                                      # excess high sleep disturbance, adults, counterfactual 40 dB
ur=R.road_2022.values.astype('float64')*1.88                # near-road NO2, ug/m3
paf=1-1/1.02**(np.clip(ur,0,None)/10)
dno2=w30*dr*paf; dstd=w30*nat30*paf                         # deaths with small-area age structure; with the national crude death rate
bc=lambda x: np.bincount(code,x,n)
SA=pd.DataFrame({'KEY_CODE':cats.values,'municipality':cats.str[:5].values,'n_receptors':np.bincount(code,minlength=n),'population':bc(w),'adults_20plus':bc(w20),'pop_30plus':bc(w30),
                 'renter_share':p_rent,'adi':ADI,'detached_share':p_det,'lnight_mean_db':bc(w*L)/bc(w),'adults_lnight45_share':bc(w20*(L>=45))/np.maximum(bc(w20),1e-12),
                 'near_road_no2_ugm3':bc(w*ur)/bc(w),'excess_hsd_adults':bc(hsd),'no2_deaths_30plus':bc(dno2),'no2_deaths_30plus_national_rate':bc(dstd)})
def quint(v):  # within-municipality quintile: small areas ranked from the lowest value, cut by cumulative population
    S=pd.DataFrame({'city':SA.municipality.values,'v':v,'w':SA.population.values}); X=S.dropna(subset=['v']).sort_values(['city','v'],kind='stable')
    cw=X.groupby('city').w.cumsum()/X.groupby('city').w.transform('sum'); q=np.zeros(n,int); q[X.index.values]=np.minimum((cw.values-1e-12)*5,4.999).astype(int)+1; return q
SA['renter_quintile']=quint(p_rent); SA['adi_quintile']=quint(ADI)
SA['housing']=np.where(p_det>=0.5,'mainly_detached',np.where(p_det<0.5,'mainly_apartment',''))
SA['excess_hsd_per_1000_adults']=SA.excess_hsd_adults/SA.adults_20plus.replace(0,np.nan)*1e3
SA['no2_deaths_national_rate_per_100000']=SA.no2_deaths_30plus_national_rate/SA.pop_30plus.replace(0,np.nan)*1e5
cols=['KEY_CODE','municipality','n_receptors','population','adults_20plus','pop_30plus','renter_share','renter_quintile','adi','adi_quintile','detached_share','housing',
      'lnight_mean_db','adults_lnight45_share','near_road_no2_ugm3','excess_hsd_adults','excess_hsd_per_1000_adults','no2_deaths_30plus','no2_deaths_30plus_national_rate','no2_deaths_national_rate_per_100000']
rd={'population':2,'adults_20plus':2,'pop_30plus':2,'renter_share':4,'adi':3,'detached_share':4,'lnight_mean_db':2,'adults_lnight45_share':4,'near_road_no2_ugm3':3,'excess_hsd_adults':4,
    'excess_hsd_per_1000_adults':3,'no2_deaths_30plus':6,'no2_deaths_30plus_national_rate':6,'no2_deaths_national_rate_per_100000':3}
SA[cols].round(rd).to_csv(f'{OUT}/small_area_burden_2020.csv',index=False)
# checks against the reported values
print('small areas',n,'population (million)',round(SA.population.sum()/1e6,2),'adults',round(SA.adults_20plus.sum()/1e6,2))
print('excess HSD (thousand)',round(SA.excess_hsd_adults.sum()/1e3,1),'NO2 deaths',round(SA.no2_deaths_30plus.sum()))
for k in ('renter_quintile','adi_quintile'):
    g=SA.groupby(k)[['population','adults_20plus','pop_30plus','excess_hsd_adults','no2_deaths_30plus_national_rate']].sum()
    print(k,'HSD per 1,000 adults',(g.excess_hsd_adults/g.adults_20plus*1e3).round(2).to_dict(),'NO2 per 100,000',(g.no2_deaths_30plus_national_rate/g.pop_30plus*1e5).round(2).to_dict())
g=SA[SA.renter_quintile>0].groupby(['housing','renter_quintile'])[['adults_20plus','excess_hsd_adults']].sum(); print((g.excess_hsd_adults/g.adults_20plus*1e3).round(2).unstack().to_string())

# sections: selection round (1 = first 1% of road length) and stand-alone benefit
O=pd.read_parquet(B+'data/alloc/v5_opt_orders.parquet'); SC=pd.read_parquet(B+'data/alloc/section_scenarios.parquet',columns=['id','KEY_CODE','rent_b10'])
assert (O.id.values==SC.id.values).all()
S=pd.DataFrame({'section_id':O.id.values,'length_m':O.len_m.round(1).values,'small_area':SC.KEY_CODE.values})
rnd=lambda c: np.where(O[c].values>=0,O[c].values+1,0)
S['pavement_benefit_alone']=O.pav_solo.round(4).values
for g in ('0','0.5','1','2','4','8','16'): S[f'pavement_round_gamma{g}']=rnd(f'pav_order_g{g}')
S['pavement_selected10_efficient']=((O.pav_order_g0>=0)&(O.pav_order_g0<10)).astype(int).values
S['pavement_selected10_equity_gamma16']=((O.pav_order_g16>=0)&(O.pav_order_g16<10)).astype(int).values
S['pavement_selected10_renter_share_only']=SC.rent_b10.astype(int).values
S['heavy_benefit_alone']=O.hv_solo.round(4).values
for g in ('0','1','2','4','8','16'): S[f'heavy_round_gamma{g}']=rnd(f'hv_order_g{g}')
S.to_csv(f'{OUT}/section_allocation.csv',index=False)
T=S.length_m.sum()
for c in ('pavement_selected10_efficient','pavement_selected10_equity_gamma16','pavement_selected10_renter_share_only'):
    print(c,int(S[c].sum()),'sections',round(S.length_m[S[c]==1].sum()/T*100,2),'% of length')
print('sections',len(S),'km',round(T/1e3))
