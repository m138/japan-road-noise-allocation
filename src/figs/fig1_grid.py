# For Fig. 1: receptors aggregated to grid cells (population-weighted mean Lnight and near-road NO2, number with excess high sleep disturbance)
import os
import numpy as np, pandas as pd, pyarrow.parquet as pq
B=os.environ.get('F_ROOT','.')+'/'
N=pq.read_table(os.environ.get('F_RECEPTORS','receptors_no2_noise_2015_2023.parquet'),columns=['rid','lon','lat','road_2022','NO2_2022']).to_pandas()
R=pd.read_parquet(B+'data/receptor/receptor_burden_v2.parquet',columns=['rid','pop_w','share20plus','Lnight_2021'])
R=R.merge(N,on='rid'); del N
w=R.pop_w.values.astype('float64'); L=np.clip(np.nan_to_num(R.Lnight_2021.values.astype('float64'),nan=40),40,None)
fW=lambda L:(19.4312-0.9336*L+0.0126*L**2)/100
hsd=w*R.share20plus.fillna(R.share20plus.median()).values*(fW(L)-fW(40))
Lraw=R.Lnight_2021.values.astype('float64'); okL=~np.isnan(Lraw)
for name,res,box in [('nat',0.03,(122.9,146.0,24.0,45.6)),('kanto',0.004,(139.2,140.2,35.3,36.0)),('kansai',0.004,(135.05,135.85,34.4,35.0))]:
    x0,x1,y0,y1=box; nx=int(np.ceil((x1-x0)/res)); ny=int(np.ceil((y1-y0)/res))
    ix=((R.lon.values-x0)/res).astype(int); iy=((R.lat.values-y0)/res).astype(int); m=(ix>=0)&(ix<nx)&(iy>=0)&(iy<ny)
    g=iy[m]*nx+ix[m]; n=nx*ny
    W=np.bincount(g,w[m],n); WL=np.bincount(g,(w*np.where(okL,10**(np.nan_to_num(Lraw)/10),0))[m],n); WLw=np.bincount(g,(w*okL)[m],n)
    WN=np.bincount(g,(w*R.road_2022.values)[m],n); H=np.bincount(g,hsd[m],n)
    np.savez_compressed(B+f'data/fig1_grid_{name}.npz',box=np.array(box),res=res,nx=nx,ny=ny,W=W,Lenergy=WL/np.maximum(WLw,1e-9),Lw=WLw,NO2road=WN/np.maximum(W,1e-9),HSD=H)
    print(name,nx,ny,round(W.sum()/1e6,2))
