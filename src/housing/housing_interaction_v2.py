# v2: adds the adjusted model without interaction and the slopes by housing type (confidence intervals of linear combinations)
# Interaction of renter share and housing type: small areas, municipality fixed effects, standard errors clustered by municipality, population weights
# Outcome: excess high sleep disturbance per 1,000 adults (counterfactual 40 dB, WHO curve). Adjusted for road density of ordinary roads, population density and apartment share
import os
import pandas as pd, numpy as np, glob, geopandas as gpd, sys, statsmodels.api as sm
B=os.environ.get('F_ROOT','.')+'/'; S=sys.argv[1]
R=pd.read_parquet(B+'data/receptor/receptor_burden_v2.parquet',columns=['KEY_CODE','pop_w','share20plus','Lnight_2021'])
fW=lambda L:(19.4312-0.9336*L+0.0126*L**2)/100
L=np.clip(np.nan_to_num(R.Lnight_2021.values.astype('float64'),nan=40),40,None)
R['h']=R.pop_w.astype('float64')*(fW(L)-fW(40)); R['w']=R.pop_w.astype('float64')
A=R.groupby('KEY_CODE')[['w','h']].sum(); del R
adi=pd.read_csv(B+'data/adi_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0)
fs=[pd.read_csv(f,encoding='cp932',skiprows=5,dtype=str) for f in sorted(glob.glob(B+'census2020_small_area/t08_1/*.csv'))]
d=pd.concat(fs,ignore_index=True); d.columns=[c.strip() for c in d.columns]; d=d[d['町丁字コード']!='-']
d['KEY_CODE']=d['市区町村コード'].str.zfill(5)+d['町丁字コード']; d=d.set_index('KEY_CODE')
num=lambda s: pd.to_numeric(s.str.replace('-','0'),errors='coerce').astype(float)
hh=pd.DataFrame({'p_detached':num(d['一戸建'])/num(d['うち主世帯']),'p_apt':num(d['共同住宅'])/num(d['うち主世帯'])}).replace([np.inf,-np.inf],np.nan)
ar=pd.concat([gpd.read_file(f)[lambda x:x.HCODE==8101][['KEY_CODE','AREA']] for f in sorted(glob.glob(f'{S}/estat_all/shp/r2ka*.shp'))]).groupby('KEY_CODE').AREA.sum()
rl=pd.read_csv(B+'data/gen_road_length_by_sa.csv',dtype={'KEY_CODE':str},index_col=0).gen_road_m
AG=pd.read_csv(B+'data/age_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0).share20plus
def mp(s):
    v=A.index.to_series().map(s); m=v.isna(); v[m]=A.index.to_series()[m].str[:9].map(s); return v
A['p_rent']=mp(adi.p_rent); A['p_detached']=mp(hh.p_detached); A['p_apt']=mp(hh.p_apt); A['s20']=mp(AG)
A['area_km2']=A.index.to_series().map(ar)/1e6; A['road_m']=A.index.to_series().map(rl).fillna(0)
A['y']=A.h*1000/(A.w*A.s20)
A['road_dens']=A.road_m/1e3/A.area_km2; A['log_popdens']=np.log(A.w/A.area_km2)
A['apt_major']=(A.p_detached<0.5).astype(float); A['city']=A.index.str[:5]
D=A.replace([np.inf,-np.inf],np.nan).dropna(subset=['y','p_rent','p_detached','road_dens','log_popdens','p_apt']); D=D[(D.w>0)]
D['rent10']=D.p_rent*10; D['rent10_x_apt']=D.rent10*D.apt_major

def fit(X):
    Z=D[['y']+X+['w','city']].copy()
    for c in ['y']+X:
        m=(Z[c]*Z.w).groupby(Z.city).transform('sum')/Z.w.groupby(Z.city).transform('sum'); Z[c]=Z[c]-m
    return sm.WLS(Z.y,Z[X],weights=Z.w).fit(cov_type='cluster',cov_kwds={'groups':pd.factorize(Z.city)[0]})
rows=[]
def add(lab,term,est,lo,hi): rows.append(dict(model=lab,term=term,coef=round(est,3),lo=round(lo,3),hi=round(hi,3)))
M=[('M1 rent only',['rent10'],False),
   ('M2 + road density, log pop density (no interaction)',['rent10','road_dens','log_popdens'],False),
   ('M3 + apartment share (no interaction)',['rent10','road_dens','log_popdens','p_apt'],False),
   ('M4 + apartment-major indicator (no interaction)',['rent10','road_dens','log_popdens','p_apt','apt_major'],False),
   ('M5 rent x apartment-major (unadjusted)',['rent10','apt_major','rent10_x_apt'],True),
   ('M6 rent x apartment-major + road density, log pop density, apartment share',['rent10','apt_major','rent10_x_apt','road_dens','log_popdens','p_apt'],True)]
for lab,X,inter in M:
    r=fit(X); ci=r.conf_int()
    for t in X: add(lab,t,r.params[t],ci.loc[t,0],ci.loc[t,1])
    if inter:
        c=np.zeros(len(X)); c[X.index('rent10')]=1; c[X.index('rent10_x_apt')]=1
        tt=r.t_test(c); cc=tt.conf_int()[0]; add(lab,'slope_in_apartment_major',float(tt.effect[0]),cc[0],cc[1])
        add(lab,'slope_in_detached_major',r.params['rent10'],ci.loc['rent10',0],ci.loc['rent10',1])
T=pd.DataFrame(rows); T['n']=len(D); T.to_csv(B+'data/v8_housing_regression.csv',index=False)
print(T.to_string()); print('n',len(D),'apt_major share of pop',(D.w*D.apt_major).sum()/D.w.sum(),'n apt',int(D.apt_major.sum()))
