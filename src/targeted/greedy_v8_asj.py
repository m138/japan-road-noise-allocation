# version 8: sensitivity analysis with the effect of low-noise (porous asphalt) pavement given by the equations of ASJ RTN-Model 2018 for ordinary roads
#  Source: Sakamoto S. Road traffic noise prediction model "ASJ RTN-Model 2018". Acoust Sci Tech 41(3), 2020 (English version of J Acoust Soc Jpn 75(4), 2019)
#   Dense asphalt (Table 2.3, non-steady flow 10<=V<=60 km/h, two vehicle categories): L_WA = a + 10 lg V; light vehicles a=82.3, heavy vehicles a=88.8
#   Porous asphalt, ordinary roads (Appendix A4, Table A4.1, non-steady flow): L_WA = a + 10 lg V + c lg(1+y); light vehicles a=76.6 c=7.3, heavy vehicles a=84.9 c=3.6 (y = years after paving)
#   -> reduction dL_light = 5.7 - 7.3 lg(1+y), dL_heavy = 3.9 - 3.6 lg(1+y) (independent of speed). The noise model uses the non-steady equation for all ordinary roads, so the same treatment is used.
#  Reduction coefficient of a section c = 1 - [(1-fl)*10^(-dL_light/10) + fl*10^(-dL_heavy/10)]; fl is the share of heavy vehicles in the nighttime energy.
#  Conditions: immediately after paving (y=0), after 3 years (y=3), and the energy mean over 10 years after paving (y=0-9). As the equation gives, the reduction for light vehicles is negative for y>5 (used as it is).
#  The preparation and functions are read from the first half of greedy_v7.py (before part 1).
# usage: greedy_v8_asj.py <M dir> <noise data dir>
import sys, numpy as np, pandas as pd, time
src=open(__file__.replace('greedy_v8_asj.py','greedy_v7.py'),encoding='utf-8').read().split('# --- 1. pavement effect')[0]
src=src.replace("src=open(__file__.replace('greedy_v7.py','greedy_v6.py')","src=open(__file__.replace('greedy_v8_asj.py','greedy_v6.py')")
exec(src)
def dL(y): return 5.7-7.3*np.log10(1+y), 3.9-3.6*np.log10(1+y)
def rfac(y): l,h=dL(y); return (1-FL)*10**(-l/10)+FL*10**(-h/10)
for y in range(10): print('y',y,'dL light %.2f heavy %.2f'%dL(y),flush=True)
SCEN={'asj_y0':1-rfac(0),'asj_y3':1-rfac(3),'asj_avg10':1-np.mean([rfac(y) for y in range(10)],axis=0)}
for scen,c in SCEN.items():
    db=-10*np.log10(1-c); print(scen,'section reduction dB pct 5/25/50/75/95',np.round(np.percentile(db,[5,25,50,75,95]),2),flush=True)
    ea=rec('asj',scen+'_all',ALL,c,scen=scen)
    for g in (0,16):
        sg=greedy_g(c,omR(g),0.10,'len'); e=rec('asj',f'{scen}_reopt_g{g}_len',sg,c,scen=scen,gamma=g); rows[-1]['share_all']=e['red']/ea['red']*100
        if g==0:
            eu=rec('asj',f'{scen}_uniformsel_g0_len',s0,c,scen=scen,gamma=0); rows[-1]['regret_pct']=(1-eu['avoided']/e['avoided'])*100; rows[-1]['overlap_len_pct']=LEN[s0&sg].sum()/LEN[s0].sum()*100; rows[-1]['share_all']=eu['red']/ea['red']*100
    sl=greedy_g(c,omR(0),0.10,'lane'); e=rec('asj',f'{scen}_reopt_g0_lane',sl,c,scen=scen,gamma=0); rows[-1]['share_all']=e['red']/ea['red']*100
pd.DataFrame(rows).round(5).to_csv(B+'data/v8_asj_checks.csv',index=False); print('done',round(time.time()-t00),flush=True)
