# Sensitivity analysis: receptors are placed on building land within each small area (National Land Numerical Information, land use subdivision mesh L03-b 2021, 100 m, building land 0700) and the matrix is built again.
# Small areas with no building-land cell keep the uniform placement. The number of receptors and their population are unchanged. The list of cells is _building_cells.parquet from the sensitivity analysis of the NO2 model.
# Matrix M of noise energy (receptors x sections), stored by prefecture, to obtain the marginal benefit of each section (version 5). All prefectures are run in order to reproduce the random numbers. Based on noise_scenarios.py, which adds receptor-level output and mitigation scenarios (heavy vehicles halved, lower motorway speed, low-noise pavement, combination) to the national calculation of the noise model (v4)
"""National night-time road noise (census roads only) with population-proportional random receptors (1 per ~10 residents)
uniformly placed in each 2020 census small area; no building shielding; region-specific bias correction from national validation.
v4: receptors that fall on a carriageway (distance to a census section centre line < carriageway width/2) are removed and the
population of the small area is divided among the remaining receptors (if all fall on carriageways, they are kept).
usage: noise_national.py <v3 data dir> <noise_lib dir> <out dir> [prefs]"""
import pandas as pd, numpy as np, geopandas as gpd, shapely, sys, glob, os, time, scipy.sparse as sp
sys.path.insert(0,sys.argv[2]); import noise_lib as nl
S=sys.argv[1]; NAT=sys.argv[3]; PREFS=sys.argv[4].split(',') if len(sys.argv)>4 else None
H=4.0; CUT=200.0; PER=10; rng=np.random.default_rng(20261106)
REG={**{p:'北海道・東北' for p in ['01','02','03','04','05','06','07']},**{p:'関東' for p in ['08','09','10','11','12','13','14']},**{p:'中部' for p in ['15','16','17','18','19','20','21','22','23']},
     **{p:'近畿' for p in ['24','25','26','27','28','29','30']},**{p:'中国・四国' for p in ['31','32','33','34','35','36','37','38','39']},**{p:'九州・沖縄' for p in ['40','41','42','43','44','45','46','47']}}
V=pd.read_csv(f'{S}/national_validation_points.csv',dtype={'pref':str})
bias_gen=V[V.grp=='gen'].groupby(V.pref.map(REG)).err.mean().to_dict(); V['rtc']=V.rt.astype(str).str[0].where(V.rt.astype(str).str[0].isin(['1','2']),'other')
bias_mw=V[V.grp=='mw'].groupby('rtc').err.mean().to_dict()   # v2: road-type-specific motorway correction
print('bias gen by region',{k:round(v,2) for k,v in bias_gen.items()},' mw',{k:round(v,2) for k,v in bias_mw.items()},flush=True)
prof=nl.night_profile(f'{S}/census27.csv')
K=nl.load_sections(f'{S}/geo/traffic_census_2021_converted.parquet')
SRC=nl.source_table(K.drop_duplicates('id'),prof).set_index('id')
# source power and speed by vehicle type (for the scenarios)
_K=K.drop_duplicates('id'); _st=(_K.motorway=='1')|_K.rt.isin(['1','2'])
_V=_K.vlim.where(_K.vlim>0); _V=np.where(_st,np.clip(_V.fillna(80),40,120),np.clip(_V.fillna(40),10,60))
_Sl=np.zeros(len(_K))
for p in ('u','d'):
    n12=(_K[f'{p}24l']-_K[f'{p}12l']).clip(lower=0).fillna(0).values
    pf=np.where(_st,prof[('mw','2')].mean(),prof[('gen','2')].mean())
    L=np.where(_st,53.2+30*np.log10(_V),88.8+10*np.log10(_V))
    _Sl+=10**((L-8)/10)*n12*pf/(1000*_V)
SRC['fl']=pd.Series(_Sl,index=_K.id.values).reindex(SRC.index).values/np.where(SRC.S>0,SRC.S,np.nan)
SRC['fl']=SRC.fl.fillna(0).clip(0,1)
SRC['spd']=np.where(_st,((np.maximum(_V-10,30))/_V)**2,1.0)  # motorway speed lowered by 10 km/h (steady flow, S proportional to V^2). Ordinary roads use the non-steady equation, which does not depend on speed
SRC.loc[:,'spd']=pd.Series(SRC.spd.values,index=SRC.index)
K=K[K.id.isin(SRC.index[SRC.ok])].dissolve(by='id').reset_index()
K['pref']=K.muni.str[:2]
rtc=K.rt.astype(str).str[0].where(K.rt.astype(str).str[0].isin(['1','2']),'other')
b=np.where(SRC.loc[K.id,'grp'].values=='mw',rtc.map(bias_mw).values,K.pref.map(REG).map(bias_gen).fillna(np.mean(list(bias_gen.values()))).values)
K['S']=SRC.loc[K.id,'S'].values*10**(-b/10); K['off']=SRC.loc[K.id,'off'].values; K['hw']=K.carw.values/2; K['fl']=SRC.loc[K.id,'fl'].values; K['spd']=SRC.loc[K.id,'spd'].values; K['grp']=SRC.loc[K.id,'grp'].values
K['geometry']=K.geometry.simplify(0.00002)
SC=pd.read_parquet(os.environ.get('F_ROOT','.')+'/data/alloc/section_scenarios.parquet',columns=['id','len_m','KEY_CODE'])
sidx=pd.Series(np.arange(len(SC)),index=SC.id.astype(str).values)
FLv=SRC.fl.reindex(SC.id.astype(str).values).fillna(0).values.astype(np.float32)
os.makedirs(NAT,exist_ok=True); np.save(f'{NAT}/section_fl.npy',FLv)
REC=pd.concat([pd.read_parquet(f,columns=['KEY_CODE','Ln']) for f in sorted(glob.glob(sys.argv[5]+'/rec_*.parquet'))]) if len(sys.argv)>5 else None
fW=lambda L:(19.4312-0.9336*L+0.0126*L**2)/100
def LLc(x): return np.clip(np.where(x>0,10*np.log10(np.maximum(x,1e-30)),40),40,None)
print('sections',len(K),flush=True)
CELLS=pd.read_parquet(os.environ.get('F_BUILDING_CELLS','_building_cells.parquet'))
os.makedirs(NAT,exist_ok=True)
files=sorted(glob.glob(f'{S}/estat_all/shp/r2ka*.shp'))
for f in files:
    pc=os.path.basename(f)[4:6]
    if PREFS and pc not in PREFS: continue
    t0=time.time()
    A=gpd.read_file(f); A=A[A.HCODE==8101][['KEY_CODE','PREF','CITY','JINKO','geometry']]
    A=A.dissolve(by='KEY_CODE',aggfunc={'PREF':'first','CITY':'first','JINKO':'sum'}).reset_index(); A=A[A.JINKO>0]
    z=nl.utm_epsg(A.to_crs(4326).union_all().centroid.x); A=A.to_crs(z)
    # receptors
    n=np.maximum(1,np.round(A.JINKO.values/PER).astype(int)); px=[];py=[];kc=[];pw=[]
    bb=A.to_crs(6668).total_bounds; Cb=CELLS[(CELLS.lon>=bb[0]-.01)&(CELLS.lon<=bb[2]+.01)&(CELLS.lat>=bb[1]-.01)&(CELLS.lat<=bb[3]+.01)]
    Cz=gpd.GeoDataFrame({'lat':Cb.lat.values},geometry=gpd.points_from_xy(Cb.lon,Cb.lat),crs=6668).to_crs(z); Cz=gpd.sjoin(Cz,A[['KEY_CODE','geometry']],predicate='within')
    Cz=pd.DataFrame({'KEY_CODE':Cz.KEY_CODE.values,'ccx':Cz.geometry.x.values,'ccy':Cz.geometry.y.values,'lat':Cz.lat.values}); Cz['dx']=4.5/3600*111320*np.cos(np.radians(Cz.lat)); DY=3/3600*111320
    cells={kk:(g.ccx.values,g.ccy.values,g.dx.values) for kk,g in Cz.groupby('KEY_CODE')}; rb=np.random.default_rng([20261123,int(pc)]); n_moved=0; p_moved=0.0
    for (key,pop,geom),k in zip(A[['KEY_CODE','JINKO','geometry']].itertuples(index=False),n):
        if key in cells:
            cx,cy,dx=cells[key]; got=np.empty((0,2)); tries=0
            while len(got)<k and tries<30:
                m=max(4*(k-len(got)),16); j=rb.integers(0,len(cx),m); c=np.c_[cx[j]+rb.uniform(-.5,.5,m)*dx[j],cy[j]+rb.uniform(-.5,.5,m)*DY]
                ins=shapely.contains_xy(geom,c[:,0],c[:,1]); got=np.vstack([got,c[ins]]); tries+=1
            if len(got)<k: j=rb.integers(0,len(cx),k-len(got)); got=np.vstack([got,np.c_[cx[j],cy[j]]])
            got=got[:k]; px.append(got[:,0]); py.append(got[:,1]); kc+= [key]*len(got); pw.append(np.full(len(got),pop/len(got))); n_moved+=len(got); p_moved+=pop; continue
        x0,y0,x1,y1=geom.bounds; got=np.empty((0,2)); tries=0
        prep=shapely.prepared.prep(geom) if hasattr(shapely,'prepared') else None
        while len(got)<k and tries<50:
            m=max(4*(k-len(got)),16); c=np.c_[rng.uniform(x0,x1,m),rng.uniform(y0,y1,m)]
            ins=shapely.contains_xy(geom,c[:,0],c[:,1]); got=np.vstack([got,c[ins]]); tries+=1
        if len(got)==0: rp=geom.representative_point(); got=np.array([[rp.x,rp.y]])
        got=got[:k]; px.append(got[:,0]); py.append(got[:,1]); kc+= [key]*len(got); pw.append(np.full(len(got),pop/len(got)))
    px=np.concatenate(px); py=np.concatenate(py); pw=np.concatenate(pw); kc=np.array(kc)
    # sections near this prefecture
    box=shapely.box(*A.total_bounds).buffer(CUT+100)
    Kz=K.to_crs(z); Kz=Kz[Kz.intersects(box)]
    # v4: remove receptors on a carriageway and re-divide the small-area population among the rest
    pts=shapely.points(px,py); ktree=shapely.STRtree(Kz.geometry.values)
    pi,ki=ktree.query(pts,predicate='dwithin',distance=float(np.nanmax(Kz.hw.values)))
    on=np.zeros(len(pts),bool); hit=shapely.distance(pts[pi],Kz.geometry.values[ki])<Kz.hw.values[ki]; on[pi[hit]]=True
    R0=pd.DataFrame({'k':kc,'on':on}); allon=R0.groupby('k').on.transform('all').values; keep=~on|allon
    n_on=int(on.sum()); p_on=pw[on].sum()
    px,py,kc=px[keep],py[keep],kc[keep]; tot=A.set_index('KEY_CODE').JINKO; cnt=pd.Series(kc).value_counts()
    pw=(pd.Series(kc).map(tot)/pd.Series(kc).map(cnt)).values
    pts=shapely.points(px,py); rtree=shapely.STRtree(pts)
    SG=nl.segments(Kz.geometry.values,Kz.id.values); src=Kz.set_index('id')
    E=np.zeros(len(pts)); rows=[];cols=[];vals=[]
    for i0 in range(0,len(SG),20000):
        sg=SG.iloc[i0:i0+20000]
        lines=shapely.linestrings(np.stack([np.c_[sg.ax,sg.ay],np.c_[sg.bx,sg.by]],axis=1))
        si,ri=rtree.query(lines,predicate='dwithin',distance=CUT)
        if len(si)==0: continue
        e=nl.seg_energy(sg.ax.values[si],sg.ay.values[si],sg.bx.values[si],sg.by.values[si],px[ri],py[ri],
                        src.S.loc[sg.id.values[si]].values,src.off.loc[sg.id.values[si]].values,H)
        np.add.at(E,ri,e)
        ci=sidx.reindex(sg.id.values[si].astype(str)).values; ok=~np.isnan(ci)
        rows.append(ri[ok]); cols.append(ci[ok].astype(np.int64)); vals.append(e[ok])
    M=sp.csr_matrix((np.concatenate(vals).astype(np.float32),(np.concatenate(rows),np.concatenate(cols))),shape=(len(pts),len(SC)))
    M.sum_duplicates()
    np.savez(f'{NAT}/M_{pc}.npz',data=M.data,indices=M.indices.astype(np.int32),indptr=M.indptr.astype(np.int64),E=E.astype(np.float64),pw=pw.astype(np.float32))
    pd.DataFrame({'KEY_CODE':kc}).to_parquet(f'{NAT}/kc_{pc}.parquet')
    print(pc,'receptors',len(pts),'nnz',M.nnz,'E>0',int((E>0).sum()),'moved_pop_pct',round(p_moved/A.JINKO.sum()*100,1),round(time.time()-t0),'s',flush=True)
