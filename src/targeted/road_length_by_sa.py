# Length of ordinary roads (census sections, motorways excluded) by 2020 small area (sections are split into pieces of about 25 m and assigned by their midpoints)
import geopandas as gpd, pandas as pd, numpy as np, shapely, glob, os, sys
S=sys.argv[1]; OUT=sys.argv[2]
sys.path.insert(0,os.environ.get('F_ROOT','.')+'/src/noise_scen'); import noise_lib as nl
K=nl.load_sections(f'{S}/geo/traffic_census_2021_converted.parquet')
st=(K.motorway=='1')|K.rt.isin(['1','2']); K=K[~st][['id','geometry']].dissolve(by='id').reset_index()
print('general sections',len(K),flush=True)
res=[]
for f in sorted(glob.glob(f'{S}/estat_all/shp/r2ka*.shp')):
    A=gpd.read_file(f); A=A[A.HCODE==8101][['KEY_CODE','geometry']].dissolve(by='KEY_CODE').reset_index()
    z=nl.utm_epsg(A.to_crs(4326).union_all().centroid.x); A=A.to_crs(z)
    Kz=K.to_crs(z); Kz=Kz[Kz.intersects(shapely.box(*A.total_bounds))]
    if len(Kz)==0: continue
    seg=Kz.geometry.explode(index_parts=False).segmentize(25)
    pts=[];lens=[]
    for geom in seg.values:
        c=np.asarray(geom.coords)
        if len(c)<2: continue
        pts.append((c[1:]+c[:-1])/2); lens.append(np.hypot(*(c[1:]-c[:-1]).T))
    pts=np.vstack(pts); lens=np.concatenate(lens)
    P=gpd.GeoDataFrame({'len':lens},geometry=gpd.points_from_xy(pts[:,0],pts[:,1]),crs=z)
    j=gpd.sjoin(P,A,how='inner',predicate='within')
    res.append(j.groupby('KEY_CODE').len.sum())
    print(os.path.basename(f),len(A),round(j.len.sum()/1e3),'km',flush=True)
L=pd.concat(res).groupby(level=0).sum(); L.rename('gen_road_m').to_csv(OUT); print('total km',round(L.sum()/1e3))
