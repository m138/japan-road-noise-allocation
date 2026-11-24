# Addition to build_age.py: shares aged 30 and over and 40-74, and the national death rate at age 30 and over (for the comparison that removes age structure)
import os
import pandas as pd, numpy as np, glob
B=os.environ.get('F_ROOT','.')+'/'
fs=[pd.read_csv(f,encoding='cp932',skiprows=4,dtype=str) for f in sorted(glob.glob(B+'census2020_small_area/t03/*.csv'))]
d=pd.concat(fs,ignore_index=True); d.columns=[c.strip() for c in d.columns]
d=d[(d.iloc[:,1]=='総数')&(d['町丁字コード']!='-')].copy()
d['KEY_CODE']=d['市区町村コード'].str.zfill(5)+d['町丁字コード']
ages=['0～4歳','5～9歳','10～14歳','15～19歳','20～24歳','25～29歳','30～34歳','35～39歳','40～44歳','45～49歳','50～54歳','55～59歳','60～64歳','65～69歳','70～74歳','75～79歳','80～84歳','85～89歳','90～94歳','95～99歳','100歳以上']
P=d.set_index('KEY_CODE')[ages].apply(lambda s: pd.to_numeric(s.replace('-','0'),errors='coerce')).astype(float)
lev=d.set_index('KEY_CODE')['地域階層レベル'].astype(int)
v=pd.read_csv(B+'data/vital_2022_deaths_age_statInfId000040067346.csv',encoding='cp932',header=None,dtype=str)
row=v[(v[1]=='年計')&(v[2]=='総数')&(v[3]=='総数')].iloc[0]
deaths=pd.to_numeric(pd.Series(row[list(range(11,32))].values)).values.astype(float)
natpop=P[lev.isin([2,4]).values].sum().values.astype(float)
rate=deaths/natpop
ext=sum(pd.to_numeric(pd.Series(v[(v[1]=='年計')&(v[2]==c)&(v[3]=='総数')].iloc[0][list(range(11,32))].values).replace('-','0')).values.astype(float) for c in ['Se32不慮の事故','Se34自殺'])
rate_na=(deaths-ext)/natpop
print('external share among 30+ deaths',round(ext[6:].sum()/deaths[6:].sum(),4))
nat30=(natpop[6:]*rate[6:]).sum()/natpop[6:].sum()
print('national death rate among 30+ (per 1000):',round(nat30*1000,3),' national share30+',round(natpop[6:].sum()/natpop.sum(),4))
tot=P.sum(axis=1)
out=pd.DataFrame({'share30plus':P[ages[6:]].sum(axis=1)/tot,'share40_74':P[ages[8:15]].sum(axis=1)/tot,
                  'dr_among30':(P[ages[6:]]*rate[6:]).sum(axis=1)/P[ages[6:]].sum(axis=1),'dr_na_among30':(P[ages[6:]]*rate_na[6:]).sum(axis=1)/P[ages[6:]].sum(axis=1)})
out.attrs['nat30']=nat30
out.to_csv(B+'data/age2_small_area_2020.csv'); pd.Series({'nat_dr_among30':nat30}).to_csv(B+'data/age2_national.csv')
print(out.describe().round(4))
