# Additions of version 10: the definitions of the analysis are those of version 8 (freeze_v8); only sensitivity analyses and summaries are added
#  muni : comparison of the disparity within municipalities. For each municipality that has both Q1 and Q5 of renter share, the rates and their difference (Q5-Q1) are obtained and
#         combined with the number of adults of the municipality (all quintiles, common to all allocations) as the weight. Before mitigation, efficient allocation (gamma=0), gamma=16, renter-share-only, all sections.
#         The difference gamma=16 minus gamma=0 gets a 95% interval from a bootstrap of municipalities (2,000 samples, seed 20261123).
#  road : bias by road type. The mean cross-validation residual (estimate minus measurement; leave-one-prefecture-out, 11,978 measurement records on ordinary roads) of each road type
#         is used as a correction of the contribution of sections of that type (contribution x 10^(-residual/10)); a diagnostic sensitivity analysis.
#         Evaluation of the sections selected in the main analysis (fixed), reselection with the corrected model, and overlap of the selected length.
#  rand : repeated random allocation. Sections are put in random order (seeds 1000-1099) and selected up to 10% of road length (same method as Fig. 3; seeds 1000-1002 are the three orders of the earlier Fig. 3).
#         Low-noise pavement (uniform 3 dB) and a 50% reduction in nighttime heavy-vehicle traffic.
#  rentq: the mean residual by renter-share quintile of the small area that contains the monitoring site (data/v9_cv_residual_by_smallarea_housing.csv) is used, according to the quintile of the receptor,
#         as a correction of the contribution of ordinary roads (contribution x 10^(-residual/10)); a diagnostic sensitivity analysis. Receptors without a quintile (renter share missing) are not corrected. Fixed and reselected.
#  popden: comparison with a simple priority rule. The adult population within 200 m of a section (receptors with a contribution in matrix M) is divided by the section length, and
#         sections are selected in descending order of this value up to 10% of road length. The rule uses neither noise contributions nor the exposure-response relationship; difference in benefit from the efficient allocation and overlap.
#  curve: for the curves of Fig. 3, the same 100 random orders evaluated in steps of 5% of road length and lane length (arguments i0 i1 split the range of seeds)
# usage: v10_review.py <M dir> <noise data dir> <muni|road|rand|curve> [i0 i1]
import os
import sys, numpy as np, pandas as pd, time
B=os.environ.get('F_ROOT','.')+'/'; TAG=sys.argv[3]
src=open(B+'src/targeted/greedy_v7.py',encoding='utf-8').read().split("cU=np.full(nsec")[0]
src=src.replace("__file__.replace('greedy_v7.py','greedy_v6.py')","B+'src/targeted/greedy_v6.py'")
exec(src)
cU=np.full(nsec,1-10**(-0.3)); cH=0.5*FL; ALL=np.ones(nsec,bool); NONE=np.zeros(nsec,bool)
O=pd.read_parquet(B+'data/alloc/v5_opt_orders.parquet'); SCN=pd.read_parquet(B+'data/alloc/section_scenarios.parquet',columns=['id','rent_b10'])
assert (O.id.values==SC.id.values).all() and (SCN.id.values==SC.id.values).all()
selo=lambda col: (O[col].values>=0)&(O[col].values<10)
s0,s16,srent=selo('pav_order_g0'),selo('pav_order_g16'),SCN.rent_b10.values.astype(bool)
def hafter(sel,c,ctx=CEN):
    x=c*np.asarray(sel,float); return np.concatenate([p['w20']*ctx['exc'](E-fg*(p['M']@x)) for p,fg,E in zip(P,ctx['fg'],ctx['E'])])
def ev(sel,c,ctx=CEN):
    h=hafter(sel,c,ctx); h0=hafter(NONE,c,ctx); r=lambda hh,k: hh[q==k].sum()/w20[q==k].sum()*1e3
    return dict(len_pct=LEN[np.asarray(sel)>0].sum()/LEN.sum()*100,HSD0=h0.sum(),HSD=h.sum(),red=(1-h.sum()/h0.sum())*100,
                ratio=r(h,5)/r(h,1),gap=r(h,5)-r(h,1),ratio0=r(h0,5)/r(h0,1),gap0=r(h0,5)-r(h0,1))
t00=time.time()
if TAG=='muni':
    nc=len(city_codes); W=lambda k: np.bincount(city,w20*(q==k),nc); Wc=np.bincount(city,w20,nc); W1c,W5c=W(1),W(5)
    inc=(W1c>0)&(W5c>0)
    print('cities',nc,'with Q1&Q5',inc.sum(),'adult share',Wc[inc].sum()/Wc.sum()*100,'Q1 share in inc',W1c[inc].sum()/W1c.sum()*100,'Q5',W5c[inc].sum()/W5c.sum()*100,flush=True)
    R={}
    for nm,sel,c in (('none',NONE,cU),('gamma0',s0,cU),('gamma16',s16,cU),('rent_only',srent,cU),('all',ALL,cU)):
        h=hafter(sel,c); R[nm]=(np.bincount(city,h*(q==1),nc),np.bincount(city,h*(q==5),nc)); print(nm,round(time.time()-t00),flush=True)
    def agg(idx,wt):
        out={}
        for nm,(h1,h5) in R.items():
            r1=h1[idx]/W1c[idx]*1e3; r5=h5[idx]/W5c[idx]*1e3; w=wt[idx]
            out[nm]=dict(gap=(w*(r5-r1)).sum()/w.sum(),ratio=(w*r5).sum()/(w*r1).sum(),r1=(w*r1).sum()/w.sum(),r5=(w*r5).sum()/w.sum())
        return out
    rows=[]; I=np.where(inc)[0]
    for wn,wt in (('adults_city',Wc),('adults_q1q5',W1c+W5c)):
        A=agg(I,wt)
        for nm,d in A.items(): rows.append(dict(weight=wn,scenario=nm,**d))
    # pooled (national aggregation) restricted to the same cities, for comparison
    for nm,(h1,h5) in R.items():
        r1=h1[I].sum()/W1c[I].sum()*1e3; r5=h5[I].sum()/W5c[I].sum()*1e3; rows.append(dict(weight='pooled_same_cities',scenario=nm,gap=r5-r1,ratio=r5/r1,r1=r1,r5=r5))
        r1=h1.sum()/W1c.sum()*1e3; r5=h5.sum()/W5c.sum()*1e3; rows.append(dict(weight='pooled_all(main)',scenario=nm,gap=r5-r1,ratio=r5/r1,r1=r1,r5=r5))
    # per-city direction
    g=lambda nm: R[nm][1][I]/W5c[I]*1e3-R[nm][0][I]/W1c[I]*1e3
    d16=g('gamma16')-g('gamma0'); d0=g('gamma0')-g('none')
    rows.append(dict(weight='city_count',scenario='gamma16_lt_gamma0_gap',gap=(d16<-1e-9).sum(),ratio=(np.abs(d16)<=1e-9).sum(),r1=(d16>1e-9).sum(),r5=Wc[I][d16<-1e-9].sum()/Wc[I].sum()*100))
    # among municipalities with the same absolute gap (|d16|<=1e-9), the number in which the Q1 and Q5 rates are each identical
    e1=np.abs(R['gamma16'][0][I]-R['gamma0'][0][I])/W1c[I]*1e3<=1e-9; e5=np.abs(R['gamma16'][1][I]-R['gamma0'][1][I])/W5c[I]*1e3<=1e-9; same=np.abs(d16)<=1e-9
    rows.append(dict(weight='city_count',scenario='same_gap_and_both_rates_equal',gap=(same&e1&e5).sum(),ratio=same.sum(),r1=(same&~(e1&e5)).sum(),r5=Wc[I][same].sum()/Wc[I].sum()*100))
    rows.append(dict(weight='city_count',scenario='gamma0_lt_none_gap',gap=(d0<-1e-9).sum(),ratio=(np.abs(d0)<=1e-9).sum(),r1=(d0>1e-9).sum(),r5=Wc[I][d0<-1e-9].sum()/Wc[I].sum()*100))
    rng=np.random.default_rng(20261123); bs=[]
    for b in range(2000):
        J=rng.choice(I,len(I),replace=True); A=agg(J,Wc)
        bs.append([A['none']['gap'],A['gamma0']['gap'],A['gamma16']['gap'],A['rent_only']['gap'],A['gamma16']['gap']-A['gamma0']['gap'],A['gamma0']['gap']-A['none']['gap'],
                   A['none']['ratio'],A['gamma0']['ratio'],A['gamma16']['ratio'],A['gamma16']['ratio']-A['gamma0']['ratio'],A['gamma16']['gap']-A['rent_only']['gap']])
    bs=np.array(bs); names=['gap_none','gap_g0','gap_g16','gap_rent','dgap_g16_g0','dgap_g0_none','ratio_none','ratio_g0','ratio_g16','dratio_g16_g0','dgap_g16_rent']
    for i,n in enumerate(names): rows.append(dict(weight='bootstrap_adults_city',scenario=n,gap=np.percentile(bs[:,i],2.5),ratio=np.percentile(bs[:,i],97.5),r1=np.mean(bs[:,i]<0)*100,r5=np.nan))
    rows.append(dict(weight='meta',scenario='n_cities_total',gap=nc,ratio=inc.sum(),r1=Wc[inc].sum()/Wc.sum()*100,r5=np.nan))
    T=pd.DataFrame(rows); T.to_csv(B+'data/v10_within_muni_gap.csv',index=False); print(T.to_string())
elif TAG=='road':
    G2=pd.read_parquet(f'{S}/geo/traffic_census_2021_converted.parquet',columns=['census','道路種別']).drop_duplicates('census')
    rt=G2.set_index(G2.census.astype(str))['道路種別'].reindex(SC.id.astype(str).values).astype(str).str[:1].values
    CV=pd.read_csv(B+'data/v9_cv_residual_by_muni_housing.csv'); CV=CV[CV.stratum=='road type']; res=dict(zip(CV.level.str[:1],CV.mean_residual))
    print('road types',pd.Series(rt).value_counts().to_dict(),'residual',res,flush=True)
    f=np.array([10**(-res.get(t,0.0)/10) for t in rt])
    ctxR=dict(fg=CEN['fg'],exc=CEN['exc'],E=[E-p['Eg']+p['M']@f for p,E in zip(P,CEN['E'])])
    cR=cU*f; rows=[]
    def rec(name,sel,c,ctx,**kw): e=ev(sel,c,ctx); e.update(name=name,**kw); rows.append(e); print(name,{k:round(v,4) for k,v in e.items() if isinstance(v,float)},round(time.time()-t00),flush=True); return e
    for model,c,ctx in (('main',cU,CEN),('roadtype_corrected',cR,ctxR)):
        ea=rec(f'{model}_all',ALL,c,ctx,model=model)
        e0=rec(f'{model}_fixed_g0',s0,c,ctx,model=model); e16=rec(f'{model}_fixed_g16',s16,c,ctx,model=model); rec(f'{model}_fixed_rent',srent,c,ctx,model=model)
        if model!='main':
            r0=greedy_g(c,omR(0),0.10,'len',ctx); r16=greedy_g(c,omR(16),0.10,'len',ctx)
            f0=rec(f'{model}_reopt_g0',r0,c,ctx,model=model); f16=rec(f'{model}_reopt_g16',r16,c,ctx,model=model)
            ov=lambda a,b: LEN[a&b].sum()/LEN[a].sum()*100
            rows.append(dict(name=f'{model}_summary',model=model,share_all_fixed_g0=e0['red']/ea['red']*100,share_all_reopt_g0=f0['red']/ea['red']*100,
                 regret_fixed_g0=(1-e0['red']/f0['red'])*100,loss_g16_fixed=(1-e16['red']/e0['red'])*100,loss_g16_reopt=(1-f16['red']/f0['red'])*100,
                 overlap_g0=ov(r0,s0),overlap_g16=ov(r16,s16),len_share_corrected_up=LEN[f>1].sum()/LEN.sum()*100,len_share_corrected_down=LEN[f<1].sum()/LEN.sum()*100))
            print(rows[-1],flush=True)
            rtsel=pd.DataFrame({'rt':rt,'len':LEN,'s0':s0,'r0':r0}); print(rtsel.groupby('rt').apply(lambda d: pd.Series({'len_all':d.len.sum()/LEN.sum()*100,'main_g0':d.len[d.s0].sum()/LEN[s0].sum()*100,'corr_g0':d.len[d.r0].sum()/LEN[r0].sum()*100})).round(1),flush=True)
    pd.DataFrame(rows).to_csv(B+'data/v10_roadtype_bias.csv',index=False)
elif TAG=='rentq':
    CV=pd.read_csv(B+'data/v9_cv_residual_by_smallarea_housing.csv'); CV=CV[CV.stratum=='within-municipality rent quintile of small area']
    res={int(l[1]):v for l,v in zip(CV.level,CV.mean_residual) if l.startswith('Q')}; print('residual by quintile',res,flush=True)
    fq=lambda qq: np.array([10**(-res.get(int(k),0.0)/10) for k in range(6)])[qq]
    ctxQ=dict(exc=CEN['exc'],fg=[fg*fq(p['q']) for p,fg in zip(P,CEN['fg'])],E=[E+(fq(p['q'])-1)*p['Eg']*fg for p,E,fg in zip(P,CEN['E'],CEN['fg'])])
    rows=[]
    def rec(name,sel,c,ctx,**kw): e=ev(sel,c,ctx); e.update(name=name,**kw); rows.append(e); print(name,{k:round(v,4) for k,v in e.items() if isinstance(v,float)},round(time.time()-t00),flush=True); return e
    for model,ctx in (('main',CEN),('rentq_corrected',ctxQ)):
        ea=rec(f'{model}_all',ALL,cU,ctx,model=model)
        e0=rec(f'{model}_fixed_g0',s0,cU,ctx,model=model); e16=rec(f'{model}_fixed_g16',s16,cU,ctx,model=model); rec(f'{model}_fixed_rent',srent,cU,ctx,model=model)
        if model!='main':
            r0=greedy_g(cU,omR(0),0.10,'len',ctx); r16=greedy_g(cU,omR(16),0.10,'len',ctx)
            f0=rec(f'{model}_reopt_g0',r0,cU,ctx,model=model); f16=rec(f'{model}_reopt_g16',r16,cU,ctx,model=model)
            ov=lambda a,b: LEN[a&b].sum()/LEN[a].sum()*100
            rows.append(dict(name=f'{model}_summary',model=model,share_all_fixed_g0=e0['red']/ea['red']*100,share_all_reopt_g0=f0['red']/ea['red']*100,
                 loss_g16_fixed=(1-e16['red']/e0['red'])*100,loss_g16_reopt=(1-f16['red']/f0['red'])*100,overlap_g0=ov(r0,s0),overlap_g16=ov(r16,s16)))
            print(rows[-1],flush=True)
    pd.DataFrame(rows).to_csv(B+'data/v10_rentq_bias.csv',index=False)
elif TAG=='popden':
    A=np.zeros(nsec)
    for p in P: A+=np.bincount(p['M'].indices,p['w20'][p['row']],nsec)
    Cl=COST['len']; TOT=Cl.sum(); rows=[]
    def pick(score):
        o=np.argsort(-score,kind='stable'); sel=np.zeros(nsec,bool); sel[o[np.cumsum(Cl[o])<=0.10*TOT]]=True; return sel
    ea=ev(ALL,cU); e0=ev(s0,cU); e16=ev(s16,cU)
    for nm,sel in (('all',ALL),('gamma0',s0),('gamma16',s16),('rent_only',srent),('adults200m_per_length',pick(A/np.maximum(Cl,1e-9))),('adults200m_total',pick(A))):
        e=ev(sel,cU); e.update(name=nm,share_all=e['red']/ea['red']*100,loss_vs_g0=(1-e['red']/e0['red'])*100,overlap_g0_len=LEN[sel&s0].sum()/LEN[s0].sum()*100,
                 avoided=e['HSD0']-e['HSD'],n_sections=int(np.asarray(sel).sum())); rows.append(e); print(nm,{k:(round(v,3) if isinstance(v,float) else v) for k,v in e.items()},round(time.time()-t00),flush=True)
    print('adults within 200 m of any section (share)',sum(p['w20'][np.diff(p['M'].indptr)>0].sum() for p in P)/w20.sum()*100)
    pd.DataFrame(rows).to_csv(B+'data/v10_popdensity_rule.csv',index=False)
elif TAG=='rand':
    rows=[]; Cl=COST['len']; TOT=Cl.sum()
    for i in range(100):
        o=np.argsort(-np.random.default_rng(1000+i).random(nsec),kind='stable'); sel=np.zeros(nsec,bool); sel[o[np.cumsum(Cl[o])<=0.10*TOT]]=True
        for meas,c in (('pav',cU),('hv',cH)):
            e=ev(sel,c); e.update(seed=1000+i,measure=meas,lane_pct=COST['lane'][sel].sum()/COST['lane'].sum()*100); rows.append(e)
        if i%10==9: print(i,round(time.time()-t00),flush=True); pd.DataFrame(rows).to_csv(B+'data/v10_random100.csv',index=False)
    pd.DataFrame(rows).to_csv(B+'data/v10_random100.csv',index=False)
elif TAG=='curve':
    # curves of Fig. 3: random orders (seeds i0 to i1-1) evaluated in steps of 5% of road length and lane length. Run in three parallel parts; v10_random_curves_*.csv are then combined
    i0,i1=int(sys.argv[4]),int(sys.argv[5]); new=[]; h0=hafter(NONE,cU); T0=h0.sum()
    for i in range(i0,i1):
        o=np.argsort(-np.random.default_rng(1000+i).random(nsec),kind='stable')
        for cost in ('len','lane'):
            C=COST[cost]; TC=C.sum(); cl=np.cumsum(C[o])
            for b in range(5,101,5):
                sel=np.zeros(nsec,bool); sel[o[cl<=b/100*TC]]=True; h=hafter(sel,cU); a1,a5=h[q==1].sum(),h[q==5].sum()
                new.append(dict(method=f'random{i}',cost=cost,gamma=np.nan,len_pct=LEN[sel].sum()/LEN.sum()*100,lane_pct=COST['lane'][sel].sum()/COST['lane'].sum()*100,
                                red=(1-h.sum()/T0)*100,ratio=(a5/W5)/(a1/W1),gap=(a5/W5-a1/W1)*1e3))
        print(i,round(time.time()-t00),flush=True); pd.DataFrame(new).to_csv(B+f'data/v10_random_curves_{i0}_{i1}.csv',index=False)
print('done',round(time.time()-t00))
