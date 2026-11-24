# Additions of version 9: the definitions of the analysis are those of version 8; only outputs are added
#  main: rates by renter-share quintile for each allocation (before and after mitigation), adults and avoided cases (Supplementary Table 1); main values with the curve capped at 65 dB; breakdown of adults assigned 40 dB
#  bldg: with the matrix of receptors placed on building land (noise_matrix_bldg.py), the disparity before mitigation and the benefit and disparity of the sections selected in the main analysis (fixed). Reselection is also calculated
# usage: v9_quintiles.py <M dir> <noise data dir> <main|bldg>
import os
import sys, numpy as np, pandas as pd, time
B=os.environ.get('F_ROOT','.')+'/'; TAG=sys.argv[3]
src=open(B+'src/targeted/greedy_v7.py',encoding='utf-8').read().split("cU=np.full(nsec")[0]
src=src.replace("__file__.replace('greedy_v7.py','greedy_v6.py')","B+'src/targeted/greedy_v6.py'")
exec(src)
cU=np.full(nsec,1-10**(-0.3)); cH=0.5*FL; ALL=np.ones(nsec,bool)
O=pd.read_parquet(B+'data/alloc/v5_opt_orders.parquet'); SCN=pd.read_parquet(B+'data/alloc/section_scenarios.parquet',columns=['id','rent_b10','random00_b10','random01_b10','random02_b10'])
assert (O.id.values==SC.id.values).all() and (SCN.id.values==SC.id.values).all()
selo=lambda col: (O[col].values>=0)&(O[col].values<10)
Wq=np.bincount(q,w20,6); exc0=CEN['exc']; excC=lambda E: fW(np.minimum(LLc(E),65))-fW(40)
def evq(sel,c,exc=exc0):
    x=c*np.asarray(sel,float); h0q=np.zeros(6); hq=np.zeros(6)
    for p,fg,E in zip(P,CEN['fg'],CEN['E']):
        h0=p['w20']*exc(E); h=p['w20']*exc(E-fg*(p['M']@x)); h0q+=np.bincount(p['q'],h0,6); hq+=np.bincount(p['q'],h,6)
    return h0q,hq
def summ(name,meas,sel,c,exc=exc0,extra=None):
    h0q,hq=evq(sel,c,exc); t0,t=h0q.sum(),hq.sum(); r0=h0q/Wq*1e3; r=hq/Wq*1e3
    d=dict(name=name,measure=meas,len_pct=LEN[np.asarray(sel)>0].sum()/LEN.sum()*100,n_sections=int(np.asarray(sel).sum()),HSD0=t0,HSD=t,avoided=t0-t,red_pct=(1-t/t0)*100,ratio=r[5]/r[1],gap=r[5]-r[1],ratio0=r0[5]/r0[1],gap0=r0[5]-r0[1])
    for k in range(1,6): d.update({f'adults_q{k}':Wq[k],f'rate0_q{k}':r0[k],f'rate_q{k}':r[k],f'avoided_q{k}':h0q[k]-hq[k]})
    if extra: d.update(extra)
    print(name,meas,{k:round(v,3) for k,v in d.items() if k in('len_pct','red_pct','ratio','gap','rate_q1','rate_q5')},round(time.time()-t00),flush=True); return d
rows=[]
if TAG=='main':
    Etot=np.concatenate([p['E'] for p in P]); L=np.where(Etot>0,10*np.log10(np.maximum(Etot,1e-30)),np.nan); A=w20.sum()
    fl=dict(adults=A,share_no_road_200m=w20[Etot<=0].sum()/A*100,share_below40_with_road=w20[(Etot>0)&(L<40)].sum()/A*100,share_floor40=w20[~(L>=40)].sum()/A*100,share_over65=w20[L>65].sum()/A*100,share_45plus=w20[L>=45].sum()/A*100)
    print(fl); pd.Series(fl).to_csv(B+'data/v9_floor40_shares.csv')
    rows.append(summ('none','pav',np.zeros(nsec),cU))
    for g in ('0','0.5','1','2','4','8','16'): rows.append(summ('gamma'+g,'pav',selo('pav_order_g'+g),cU))
    rows.append(summ('rent_only','pav',SCN.rent_b10.values,cU))
    for i in range(3): rows.append(summ(f'random0{i}','pav',SCN[f'random0{i}_b10'].values,cU))
    sgap=greedy_g(cU,omGap,0.10,'len'); rows.append(summ('gap_only','pav',sgap,cU)); rows.append(summ('all','pav',ALL,cU))
    for g in ('0','1','4','16'): rows.append(summ('gamma'+g,'hv',selo('hv_order_g'+g),cH))
    rows.append(summ('all','hv',ALL,cH))
    T=pd.DataFrame(rows); T.to_csv(B+'data/v9_quintile_rates.csv',index=False)
    rc=[summ(n,'pav',s,cU,excC) for n,s in (('none',np.zeros(nsec)),('gamma0',selo('pav_order_g0')),('gamma16',selo('pav_order_g16')),('rent_only',SCN.rent_b10.values),('all',ALL))]
    pd.DataFrame(rc).to_csv(B+'data/v9_cap65.csv',index=False)
else:
    s0c,s16c=selo('pav_order_g0'),selo('pav_order_g16')
    for n,s in (('none',np.zeros(nsec)),('fixed_gamma0',s0c),('fixed_gamma16',s16c),('fixed_rent_only',SCN.rent_b10.values),('all',ALL)): rows.append(summ(n,'pav',s,cU))
    r0=greedy_g(cU,omR(0),0.10,'len'); r16=greedy_g(cU,omR(16),0.10,'len')
    ov=lambda a,b: LEN[a&b].sum()/LEN[a].sum()*100
    rows.append(summ('reopt_gamma0','pav',r0,cU,extra=dict(overlap_with_main=ov(r0,s0c)))); rows.append(summ('reopt_gamma16','pav',r16,cU,extra=dict(overlap_with_main=ov(r16,s16c))))
    Etot=np.concatenate([p['E'] for p in P]); L=np.where(Etot>0,10*np.log10(np.maximum(Etot,1e-30)),np.nan); A=w20.sum()
    print(dict(share_floor40=w20[~(L>=40)].sum()/A*100,share_45plus=w20[L>=45].sum()/A*100,share_over65=w20[L>65].sum()/A*100))
    pd.DataFrame(rows).to_csv(B+'data/v9_bldg_placement.csv',index=False)
