# Monte Carlo: with the treated sections fixed, propagate the uncertainty of the exposure error common to an area, the exposure-response relationship and the pavement effect
# exposure error (estimated from the national comparison with measurements, dB): ordinary roads prefecture SD 0.97 + municipality SD 1.37; motorways prefecture SD 2.87
# exposure-response: the excess of the WHO curve is scaled by the ratio of slopes of a logistic curve with OR 2.13 (1.82-2.48) per 10 dB
# pavement effect: triangular distribution (1.5, 3, 4.5 dB), assumed. Heavy vehicles are fixed at a 50% reduction at night (definition of the scenario)
# version 5: separate allocations for pavement and heavy vehicles (columns Ep_* and Eh_*); random seed and output tag are arguments (for parallel runs). draw=-1 is the Japanese exposure-response curve (central values)
# usage: mc_alloc_v5.py <mcrec dir> <ND> <seed> <tag>
import os
import pandas as pd, numpy as np, glob, sys
B=os.environ.get('F_ROOT','.')+'/'; NAT=sys.argv[1]; ND=int(sys.argv[2]); SEED=int(sys.argv[3]); TAG=sys.argv[4]
R=pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(NAT+'/mcrec_*.parquet'))],ignore_index=True)
PK=[c[3:] for c in R.columns if c.startswith('Ep_')]; HK=[c[3:] for c in R.columns if c.startswith('Eh_')]
AG=pd.read_csv(B+'data/age_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0).share20plus
adi=pd.read_csv(B+'data/adi_small_area_2020.csv',dtype={'KEY_CODE':str},index_col=0).p_rent
R['KEY_CODE']=R.KEY_CODE.astype('category'); cats=pd.Series(R.KEY_CODE.cat.categories); code=R.KEY_CODE.cat.codes.values
def mp(s):
    v=cats.map(s); m=v.isna(); v[m]=cats[m].str[:9].map(s); return v.values
s20=mp(AG); s20=np.where(np.isnan(s20),np.nanmedian(s20),s20)
w20=(R.w.values.astype('float64')*s20[code])
# within-city rent quintile by small area (population weighted)
SA=pd.DataFrame({'KEY_CODE':cats.values,'w':np.bincount(code,R.w.values.astype('float64'),len(cats)),'p_rent':mp(adi)})
SA['city']=SA.KEY_CODE.str[:5]; X=SA.dropna(subset=['p_rent']).sort_values(['city','p_rent'])
cw=X.groupby('city').w.cumsum()/X.groupby('city').w.transform('sum'); SA['q']=0; SA.loc[X.index,'q']=np.minimum((cw.values-1e-12)*5,4.999).astype(int)+1
q=SA.q.values[code]; q1=q==1; q5=q==5
pref=cats.str[:2].values[code].astype(int); city_codes,city=np.unique(cats.str[:5].values[code],return_inverse=True)
E=R.E.values.astype('float64'); Eg=R.Egen.values.astype('float64'); Em=np.maximum(E-Eg,0)
# for speed, keep only the receptors that change in each scenario (the others are the same as before mitigation)
SCN=[('pav',k) for k in PK]+[('hv',k) for k in HK]; IX={}; DE={}
for t,k in SCN:
    a=R[('Ep_' if t=='pav' else 'Eh_')+k].values; ix=np.nonzero(a>0)[0]; IX[(t,k)]=ix; DE[(t,k)]=a[ix].astype('float64')
del R
fW=lambda L:(19.4312-0.9336*L+0.0126*L**2)/100
def LLc(x): return np.clip(np.where(x>0,10*np.log10(np.maximum(x,1e-30)),40),40,None)
lg=lambda L,b: 1/(1+np.exp(-(np.log(fW(40)/(1-fW(40)))+b*(L-40)/10)))
bc=np.log(2.13); bse=(np.log(2.48)-np.log(1.82))/(2*1.96)
fJ=lambda L:1/(1+np.exp(5.6431-0.0602*L))
def exc(L,b):
    if b is None: return fJ(L)-fJ(40)
    base=fW(L)-fW(40); den=lg(L,bc)-fW(40); num=lg(L,b)-fW(40)
    return base*np.where(den>1e-9,num/np.maximum(den,1e-9),1.0)
def summ(h):
    t=(w20*h).sum(); r1=(w20*h)[q1].sum()/w20[q1].sum(); r5=(w20*h)[q5].sum()/w20[q5].sum(); return t,r5-r1,r5/r1
rng=np.random.default_rng(SEED); rows=[]
for d in ([-1,0] if SEED==20261122 else [])+list(range(1,ND+1)):
    if d<=0: dg=np.zeros(len(E)); dm=np.zeros(len(E)); b=bc if d==0 else None; pav=3.0
    else:
        gp=rng.normal(0,0.97,48); gc=rng.normal(0,1.37,len(city_codes)); mpf=rng.normal(0,2.87,48)
        dg=gp[pref]+gc[city]; dm=mpf[pref]; b=rng.normal(bc,bse); pav=rng.triangular(1.5,3.0,4.5)
    fg=10**(dg/10); fm=10**(dm/10); Eb=Eg*fg+Em*fm
    h0=exc(LLc(Eb),b); t0,d0,r0=summ(h0); row=dict(draw=d,seed=SEED,b=b,pav=pav,H_base=t0,diff_base=d0,ratio_base=r0)
    red=1-10**(-pav/10)
    S0=w20*h0; T1=S0[q1].sum(); T5=S0[q5].sum(); W1=w20[q1].sum(); W5=w20[q5].sum()
    for t,k in SCN:
        ix=IX[(t,k)]; Ebk=Eb[ix]; fgk=fg[ix]; wk=w20[ix]; q1k=q1[ix]; q5k=q5[ix]; h0k=h0[ix]
        Ek=Ebk-(red if t=='pav' else 0.5)*DE[(t,k)]*fgk
        dlt=wk*(exc(LLc(Ek),b)-h0k); tt=t0+dlt.sum(); r1=(T1+dlt[q1k].sum())/W1; r5=(T5+dlt[q5k].sum())/W5
        row.update({f'red_{t}_{k}':(1-tt/t0)*100,f'diff_{t}_{k}':r5-r1,f'ratio_{t}_{k}':r5/r1})
    rows.append(row)
    if d%25==0: print(d,round(t0/1e3),round(row['red_pav_opt_g0_b10'],2),flush=True); pd.DataFrame(rows).to_csv(B+f'data/v5_mc_draws_{TAG}.csv',index=False)
pd.DataFrame(rows).to_csv(B+f'data/v5_mc_draws_{TAG}.csv',index=False); print('done')
