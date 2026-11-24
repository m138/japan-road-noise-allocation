# Cross-validation residuals at monitoring sites on ordinary roads (estimate minus measurement, leaving out one prefecture at a time, correction by regional block),
# summarized by housing composition of the municipality (apartment share, renter share), road type and estimated level.
# The coordinates of the sites and their distance from the carriageway edge are not kept in national_validation_points.csv, so the summaries by small area and by distance use the original data (Environmental Observatory) downloaded again.
# usage: cv_residual_strata.py <national_validation_points.csv>
import os
import sys, glob, numpy as np, pandas as pd
B=os.environ.get('F_ROOT','.')+'/'
R=pd.read_csv(sys.argv[1],dtype={'pref':str,'jis':str}); R['pref']=R.pref.str.zfill(2); R['jis']=R.jis.str.zfill(5); R=R[R.grp=='gen'].dropna(subset=['est','meas']).copy()
BL={'北海道・東北':range(1,8),'関東':range(8,15),'中部':range(15,24),'近畿':range(24,31),'中国・四国':range(31,40),'九州・沖縄':range(40,48)}
reg={f'{p:02d}':k for k,v in BL.items() for p in v}; R['blk']=R.pref.map(reg)
out=[]
for p in R.pref.unique():
    tr=R[R.pref!=p]; te=R[R.pref==p].copy(); c=te.blk.map(tr.groupby('blk').err.mean()).fillna(tr.err.mean()); te['pred']=te.est-c; te['res']=te.pred-te.meas; out.append(te)
T=pd.concat(out); print('n',len(T),'bias',T.res.mean().round(3),'rmse',np.sqrt((T.res**2).mean()).round(3))
num=lambda s: pd.to_numeric(s,errors='coerce')
fs=[pd.read_csv(f,encoding='cp932',skiprows=5,dtype=str) for f in sorted(glob.glob(B+'census2020_small_area/t08_1/*.csv'))]
d=pd.concat(fs); d.columns=[c.strip() for c in d.columns]; print(d.columns.tolist()[:14])
d['city']=d['市区町村コード'].str.zfill(5); d=d[d['地域階層レベル'].str.strip()=='1']
H=pd.DataFrame({'city':d.city.values,'main':num(d['うち主世帯']).values,'det':num(d['一戸建']).values,'apt':num(d['共同住宅']).values}).groupby('city').sum()
H['p_apt']=H.apt/H.main
A=pd.read_csv(B+'data/adi_small_area_2020.csv',dtype={'KEY_CODE':str,'city':str},usecols=['KEY_CODE','level','city','H_dwell','rent'])
A=A[A.KEY_CODE.str.len()>=9]; Ar=A.groupby(A.KEY_CODE.str[:5])[['H_dwell','rent']].sum(); H['p_rent']=Ar.rent/Ar.H_dwell
Hs=H[['main','apt']].join(Ar)   # city codes of designated cities (sum of wards)
for j in sorted(set(T.jis)-set(H.index)):
    k=int(j); w=[c for c in Hs.index if k<int(c)<k+(30 if j.endswith('00') else 10)]
    if w: t=Hs.loc[w].sum(); H.loc[j,['p_apt','p_rent']]=[t.apt/t.main,t.rent/t.H_dwell]
T['p_apt']=T.jis.map(H.p_apt); T['p_rent']=T.jis.map(H.p_rent); print('unmatched jis',T.p_apt.isna().sum(),T[T.p_apt.isna()].jis.unique()[:10])
def tab(key,lab):
    g=T.groupby(key,observed=True).apply(lambda x: pd.Series({'n':len(x),'mean_meas':x.meas.mean(),'mean_pred':x.pred.mean(),'mean_residual':x.res.mean(),'se':x.res.std()/np.sqrt(len(x)),'rmse':np.sqrt((x.res**2).mean())}))
    g.insert(0,'stratum',lab); g.index=g.index.astype(str); g.index.name='level'; return g.reset_index()
O=pd.concat([tab(pd.cut(T.p_apt,[0,.2,.3,.4,.5,.6,1.0001],right=False),'municipal apartment share'),
             tab(pd.cut(T.p_rent,[0,.2,.3,.4,.5,1.0001],right=False),'municipal rent share'),
             tab(T.rt,'road type'),
             tab(pd.cut(T.pred,[0,55,60,65,70,100],right=False),'predicted Lnight (dB)'),
             tab(pd.Series('all',index=T.index),'all')])
pd.set_option('display.width',200); print(O.round(2).to_string(index=False)); O.round(4).to_csv(B+'data/v9_cv_residual_by_muni_housing.csv',index=False,encoding='utf-8-sig')
import statsmodels.api as sm
D=T.dropna(subset=['p_apt']); X=sm.add_constant(pd.DataFrame({'p_apt10':D.p_apt*10})); r=sm.OLS(D.res,X).fit(cov_type='cluster',cov_kwds={'groups':D.jis}); print(r.params.round(3).to_dict(),r.conf_int().round(3).values.tolist())
X=sm.add_constant(pd.DataFrame({'p_rent10':D.p_rent*10})); r=sm.OLS(D.res,X).fit(cov_type='cluster',cov_kwds={'groups':D.jis}); print(r.params.round(3).to_dict(),r.conf_int().round(3).values.tolist())

# The coordinates and distance from the carriageway edge are joined from the original data of the Environmental Observatory (souon*.csv in tn2019-2023.zip) by (pref, year, point, jis), and
# the residuals are stratified by (1) distance from the carriageway edge and (2) housing composition of the small area that contains the site (detached houses half or more of main households; within-municipality quintile of renter share).
# usage: cv_residual_strata.py <national_validation_points.csv> <folder of the extracted Environmental Observatory data> <folder of e-Stat small-area shapefiles>
if len(sys.argv)>3:
    import geopandas as gpd
    def rd(f):
        for e in ('utf-8-sig','cp932'):
            try: return pd.read_csv(f,encoding=e,dtype=str)
            except UnicodeDecodeError: pass
    S=pd.concat([rd(f) for f in sorted(glob.glob(sys.argv[2]+'/**/souon*.csv',recursive=True))]); S.columns=[c.strip() for c in S.columns]
    S['pref']=S.pref_code.str.zfill(2); S['year']=S.report_year.astype(int); S['point']=pd.to_numeric(S.noise_measure_point_code,errors='coerce'); S['jis']=S.jis_code.str.zfill(5)
    for c in ['longitude','latitude','distance_from_roadway']: S[c]=pd.to_numeric(S[c],errors='coerce')
    S=S.drop_duplicates(['pref','year','point','jis'])
    T['year']=T.year.astype(int); T['point']=T.point.astype(float); S['point']=S.point.astype(float)
    T=T.merge(S[['pref','year','point','jis','longitude','latitude','distance_from_roadway']],on=['pref','year','point','jis'],how='left')
    print('matched',T.longitude.notna().sum(),'of',len(T))
    # spatial join to small areas (HCODE 8101)
    G=pd.concat([gpd.read_file(f)[lambda d:d.HCODE==8101][['KEY_CODE','geometry']].to_crs(4326) for f in sorted(glob.glob(sys.argv[3]+'/r2ka*.shp'))])
    P=gpd.GeoDataFrame(T[['longitude','latitude']],geometry=gpd.points_from_xy(T.longitude,T.latitude),crs=4326)
    J=gpd.sjoin(P,G,how='left',predicate='within'); J=J[~J.index.duplicated()]; T['KEY_CODE']=J.KEY_CODE.reindex(T.index).values
    print('in small area',T.KEY_CODE.notna().sum())
    # within-municipality quintiles as in the main analysis (cumulative population of receptors; ties in renter share in area-code order) and housing type
    Rr=pd.read_parquet(B+'data/receptor/receptor_burden_v2.parquet',columns=['KEY_CODE','pop_w'])
    SA=Rr.groupby(Rr.KEY_CODE.astype(str)).pop_w.sum().rename('w').reset_index(); del Rr
    adi=pd.read_csv(B+'data/adi_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0).p_rent
    mp=lambda k,s: k.map(s).fillna(k.str[:9].map(s))
    SA['p_rent']=mp(SA.KEY_CODE,adi); SA['city']=SA.KEY_CODE.str[:5]
    X=SA.dropna(subset=['p_rent']).sort_values(['city','p_rent'],kind='stable'); cw=X.groupby('city').w.cumsum()/X.groupby('city').w.transform('sum')
    SA['q']=0; SA.loc[X.index,'q']=np.minimum((cw.values-1e-12)*5,4.999).astype(int)+1
    fs=[pd.read_csv(f,encoding='cp932',skiprows=5,dtype=str) for f in sorted(glob.glob(B+'census2020_small_area/t08_1/*.csv'))]
    d=pd.concat(fs,ignore_index=True); d.columns=[c.strip() for c in d.columns]; d=d[d['町丁字コード']!='-']
    d['KEY_CODE']=d['市区町村コード'].str.zfill(5)+d['町丁字コード']; d=d.drop_duplicates('KEY_CODE').set_index('KEY_CODE')
    n2=lambda s: pd.to_numeric(s.str.replace('-','0'),errors='coerce').astype(float)
    pdet=(n2(d['一戸建'])/n2(d['うち主世帯'])).replace([np.inf,-np.inf],np.nan)
    k=T.KEY_CODE.astype(str)
    T['q']=k.map(SA.set_index('KEY_CODE').q).fillna(mp(k.str[:9],SA.drop_duplicates('KEY_CODE').set_index('KEY_CODE').q) if False else np.nan)
    qq=SA.set_index('KEY_CODE').q; T['q']=k.map(qq); m=T.q.isna()&T.KEY_CODE.notna()
    # small-area codes of receptors have 11 or 9 digits, so the 9-digit code is tried when the 11-digit code is not found
    T.loc[m,'q']=k[m].str[:9].map(qq)
    T['p_det']=mp(k,pdet).where(T.KEY_CODE.notna())
    T['house']=np.where(T.p_det>=0.5,'detached-major',np.where(T.p_det<0.5,'apartment-major','no households'))
    T['rentq']=T.q.fillna(0).astype(int).map(lambda v: f'Q{v}' if v>0 else 'no rent share')
    O2=pd.concat([tab(pd.cut(T.distance_from_roadway,[-0.01,0.5,2,5,10,20,1e9],labels=['0 m','0.5-2 m','2-5 m','5-10 m','10-20 m','20 m+']),'distance from roadway edge'),
                  tab(T.house,'small-area housing type'),
                  tab(T.rentq,'within-municipality rent quintile of small area'),
                  tab(T.house+' / '+T.rentq,'housing type x rent quintile'),
                  tab(pd.Series('all',index=T.index),'all')])
    print(O2.round(2).to_string(index=False))
    O2[O2.stratum=='distance from roadway edge'].round(4).to_csv(B+'data/v9_cv_residual_by_distance.csv',index=False,encoding='utf-8-sig')
    O2[O2.stratum!='distance from roadway edge'].round(4).to_csv(B+'data/v9_cv_residual_by_smallarea_housing.csv',index=False,encoding='utf-8-sig')
    D=T[T.q>0].copy(); X=sm.add_constant(pd.DataFrame({'q':D.q})); r=sm.OLS(D.res,X).fit(cov_type='cluster',cov_kwds={'groups':D.jis}); print('trend per quintile',r.params.round(3).to_dict(),r.conf_int().round(3).values.tolist())
    D=T[T.p_det.notna()].copy(); X=sm.add_constant(pd.DataFrame({'apt_major':(D.p_det<0.5).astype(float)})); r=sm.OLS(D.res,X).fit(cov_type='cluster',cov_kwds={'groups':D.jis}); print('apartment-major diff',r.params.round(3).to_dict(),r.conf_int().round(3).values.tolist())
    print(T.distance_from_roadway.describe().round(2).to_dict())
    # check whether the slope across quintiles remains when road type is added as fixed effects (whether the bias by quintile only repeats the bias by road type)
    D=T[T.q>0].copy(); X=pd.get_dummies(D.rt.str[:1],prefix='rt',drop_first=True,dtype=float); X['q']=D.q; X=sm.add_constant(X)
    r=sm.OLS(D.res,X).fit(cov_type='cluster',cov_kwds={'groups':D.jis}); print('trend per quintile | road type FE',round(r.params['q'],3),r.conf_int().loc['q'].round(3).tolist())
    print(pd.crosstab(D.q,D.rt.str[:1],normalize='index').round(3))
    T[['pref','year','point','jis','rt','distance_from_roadway','KEY_CODE','q','p_det','meas','pred','res']].to_csv(B+'data/v9_cv_points_smallarea.csv',index=False,encoding='utf-8-sig')
