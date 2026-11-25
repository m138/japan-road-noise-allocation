# version 5: fall in NOx emissions from a 50% reduction in nighttime (22:00-06:00) heavy-vehicle traffic, calculated at each receptor with the road kernel of the NO2 model (1 km, h=20 m)
#  NOx of a section (g/km/day) = light 24h x EF light + heavy 24h x EF heavy (NILIM Technical Note 671, year 2020, daytime speed). Nighttime heavy vehicles = (24h - daytime 12h) x night share (hourly pattern of the noise model)
#  dq = 0.5 x nighttime heavy vehicles x EF heavy. Output: baseline C at each receptor and dC/C of each scenario (the near-road term is proportional to sqrt(C), so it is multiplied by sqrt(1-dC/C))
# usage: no2_heavy.py <noise data dir> <out dir> [prefs]
import pandas as pd, numpy as np, geopandas as gpd, shapely, sys, os, time, pyproj
sys.path.insert(0,os.environ.get('F_NO2_SRC','japan-traffic-no2/src')); import no2_lib as nl2
sys.path.insert(0,os.environ.get('F_ROOT','.')+'/src/noise_scen'); import noise_lib as nln
S=sys.argv[1]; OUT=sys.argv[2]; PREFS=sys.argv[3].split(',') if len(sys.argv)>3 else [f'{i:02d}' for i in range(1,48)]
B=os.environ.get('F_ROOT','.')+'/'; os.makedirs(OUT,exist_ok=True)
GP=f'{S}/geo/traffic_census_2021_converted.parquet'; EFP=os.environ.get('F_NO2_SRC','japan-traffic-no2/src')+'/../inputs/nox_ef_nilim671.csv'
K=nl2.load_sections(GP,EFP,2020); K=K[K.ok].copy()
N=nln.load_sections(GP).drop_duplicates('id').set_index('id'); prof=nln.night_profile(f'{S}/census27.csv')
st=((K.motorway=='1')|K.rt.isin(['1','2'])).values
nshare=np.where(st,prof[('mw','2')].sum(),prof[('gen','2')].sum())
n12=sum((N[f'{p}24l']-N[f'{p}12l']).clip(lower=0).fillna(0) for p in ('u','d')).reindex(K.id.values).fillna(0).values
K['dq_unit']=0.5*n12*nshare*nl2.ef(K.Vl.values,'large',2020,EFP)     # g/km/day removed if selected
SC=pd.read_parquet(B+'data/alloc/section_scenarios.parquet',columns=['id','random00_b10','rent_b10','all']); SC['id']=SC.id.astype(str)
OR=pd.read_parquet(B+'data/alloc/v5_opt_orders.parquet'); OR['id']=OR.id.astype(str)
sel={'opt_g0_b10':OR.set_index('id').hv_order_g0.between(0,9),'opt_g1_b10':OR.set_index('id').hv_order_g1.between(0,9),'opt_g4_b10':OR.set_index('id').hv_order_g4.between(0,9),'opt_g16_b10':OR.set_index('id').hv_order_g16.between(0,9),
     'old_rent_b10':SC.set_index('id').rent_b10.astype(bool),'old_random00_b10':SC.set_index('id').random00_b10.astype(bool),'all_b100':SC.set_index('id')['all'].astype(bool)}
qcols=['nox_gkmd']
for k,s in sel.items(): K['dq_'+k]=K.dq_unit*s.reindex(K.id.values).fillna(False).values.astype(float); qcols.append('dq_'+k)
for k in sel: print(k,'removed share of census NOx %',round(K['dq_'+k].mul(K.geometry.to_crs(6677).length/1000).sum()/K.nox_gkmd.mul(K.geometry.to_crs(6677).length/1000).sum()*100,3),flush=True)
R=pd.read_parquet(os.environ.get('F_RECEPTORS','receptors_no2_noise_2015_2023.parquet'),columns=['rid','KEY_CODE','lon','lat'])
R['pc']=R.KEY_CODE.astype(str).str[:2]
for pc in PREFS:
    o=f'{OUT}/no2h_{pc}.parquet'
    if os.path.exists(o): continue
    t0=time.time(); r=R[R.pc==pc]; z=nl2.utm_epsg(r.lon.mean())
    x,y=pyproj.Transformer.from_crs(4326,z,always_xy=True).transform(r.lon.values,r.lat.values)
    Kz=K.to_crs(z); Kz=Kz[Kz.intersects(shapely.box(x.min(),y.min(),x.max(),y.max()).buffer(1100))]
    f=nl2.road_field_multi(np.asarray(x),np.asarray(y),Kz,qcols,cut=1000,h=20)
    D=pd.DataFrame({'rid':r.rid.values,'C':f[:,0].astype(np.float32)})
    for j,k in enumerate(sel): D['rel_'+k]=np.where(f[:,0]>0,f[:,j+1]/np.maximum(f[:,0],1e-30),0).astype(np.float32)
    D.to_parquet(o); print(pc,len(r),len(Kz),round(time.time()-t0),'s',flush=True)
