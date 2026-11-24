# Tropical nights with the urban correction of the tropical-nights paper (main specification without a constant, 'u') at each receptor (2018-2023, each year and the mean)
import os
import pandas as pd,numpy as np,pyarrow.parquet as pq,pyarrow as pa,sys
sys.path.insert(0,'hn'); from idw import idw_daily; from lu_lookup import feats
b=pd.read_csv('hn/uhi_coef2.csv',index_col=0).iloc[:,0]
cols=['built_1500m','built_10km','water_4500m']
G=pd.read_csv('hn/popgrid_1km.csv',dtype={'code':str}); G['code']=G.code.str.zfill(4)
dem=pd.read_csv('dem/dem_raw.csv',usecols=['G04a_001','G04a_002']); c=dem.G04a_001.astype(str).str.zfill(8)
dem['iy']=c.str[:2].astype(int)*80+c.str[4].astype(int)*10+c.str[6].astype(int)
dem['ix']=(c.str[2:4].astype(int)+100)*80+c.str[5].astype(int)*10+c.str[7].astype(int)
dem['elev']=pd.to_numeric(dem.G04a_002,errors='coerce'); dem.loc[dem.elev<-100,'elev']=np.nan
G=G.merge(dem[['iy','ix','elev']],on=['iy','ix'],how='left'); G['elev']=G.elev.fillna(0); G['cid']=np.arange(len(G))
R=pq.read_table(os.environ.get('F_RECEPTORS','receptors_no2_noise_2015_2023.parquet'),columns=['rid','lat','lon','pop_w']).to_pandas()
R['iy']=np.floor(R.lat*120).astype(int); R['ix']=np.floor(R.lon*80).astype(int)
R=R.merge(G[['iy','ix','cid']],on=['iy','ix'],how='left'); print('receptors',len(R),'unmatched',R.cid.isna().sum(),flush=True)
R['cid']=R.cid.fillna(-1).astype(int)
fr=feats(R.lat.values,R.lon.values)
for k in cols: R[k]=np.nan_to_num(fr[k])
out={'rid':R.rid.values}
acc=np.zeros(len(R),np.float32)
for y in range(2018,2024):
    T,ja,W,I,meta=idw_daily(G.lat.values,G.lon.values,G.elev.values,y)
    fs=feats(meta.lat.values,meta.lon.values)
    S={k:(W*np.nan_to_num(fs[k])[I]).sum(1) for k in cols}
    cid=np.where(R.cid.values<0,0,R.cid.values)
    d=sum(b['d_'+k]*(R[k].values-S[k][cid]) for k in cols)
    q=np.round(d/0.05).astype(np.int32)
    key=pd.DataFrame({'cid':cid,'q':q}); u,inv=np.unique(key.values,axis=0,return_inverse=True)
    Tc=T[u[:,0]]; tn=((Tc+(u[:,1]*0.05)[:,None])>=25).sum(1).astype(np.float32)
    v=tn[inv.ravel()]; v[R.cid.values<0]=np.nan
    out[f'tn25_{y}']=v; acc+=np.nan_to_num(v)
    print(y,'pop-wtd',round(float(np.nansum(v*R.pop_w)/R.pop_w[~np.isnan(v)].sum()),2),flush=True)
out['tn25_mean_2018_2023']=acc/6
pq.write_table(pa.table(out),'receptor_tn25_2018_2023.parquet'); print('saved')
