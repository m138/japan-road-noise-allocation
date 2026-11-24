# Age structure of small areas (2020 Population Census, Table 3) and age-specific death rates of 2022 (Vital Statistics, monthly report annual total, Table 6-1, statInfId 000040067346):
# share aged 20 and over and expected death rate at age 30 and over (per person) for each small area
import os
import pandas as pd, numpy as np, glob
B=os.environ.get('F_ROOT','.')+'/'
fs=[pd.read_csv(f,encoding='cp932',skiprows=4,dtype=str) for f in sorted(glob.glob(B+'census2020_small_area/t03/*.csv'))]
d=pd.concat(fs,ignore_index=True); d.columns=[c.strip() for c in d.columns]
d=d[(d.iloc[:,1]=='総数')&(d['町丁字コード']!='-')].copy()
d['KEY_CODE']=d['市区町村コード'].str.zfill(5)+d['町丁字コード']
ages=['0～4歳','5～9歳','10～14歳','15～19歳','20～24歳','25～29歳','30～34歳','35～39歳','40～44歳','45～49歳','50～54歳','55～59歳','60～64歳','65～69歳','70～74歳','75～79歳','80～84歳','85～89歳','90～94歳','95～99歳','100歳以上']
P=d.set_index('KEY_CODE')[ages].apply(lambda s: pd.to_numeric(s.replace('-','0'),errors='coerce')).astype(float)
# national death rates by age 2022 (deaths 2022 / census 2020 population by age; approximation)
v=pd.read_csv(B+'data/vital_2022_deaths_age_statInfId000040067346.csv',encoding='cp932',header=None,dtype=str)
row=v[(v[1]=='年計')&(v[2]=='総数')&(v[3]=='総数')].iloc[0]
dcols=list(range(11,32))  # 0-04 ... 100+
deaths=pd.to_numeric(pd.Series(row[dcols].values)).values.astype(float)
lev=d.set_index('KEY_CODE')['地域階層レベル'].astype(int)
natpop=P[lev.isin([2,4]).values].sum().values.astype(float)  # sum only the lowest level (2 = town without chome, 4 = chome)  # small-area sum (includes some suppressed NaN -> ignored)
rate=deaths/natpop
print(pd.DataFrame({'age':ages,'deaths':deaths,'pop':natpop,'rate_per1000':(rate*1000).round(2)}).to_string())
tot=P.sum(axis=1)
out=pd.DataFrame({'pop_age_total':tot,'share20plus':P[ages[4:]].sum(axis=1)/tot,'share0_14':P[ages[:3]].sum(axis=1)/tot,
    'death_rate30plus':(P[ages[6:]]*rate[6:]).sum(axis=1)/tot,'death_rate_all':(P*rate).sum(axis=1)/tot})
out.to_csv(B+'data/age_small_area_2020.csv'); print(out.describe().round(4))
