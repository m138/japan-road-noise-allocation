# Areal deprivation index of small areas (weights of Nakaya 2011) from the small-area tabulation of the 2020 Population Census
# input: census2020_small_area/{t06_1,t07_1,t09,t12}/PP.csv  output: data/adi_small_area_2020.csv
import pandas as pd, numpy as np, glob, os
B=os.environ.get('F_ROOT','.')+'/'
def rd(tab, has_type):
    fs=[]
    for f in sorted(glob.glob(f'{B}census2020_small_area/{tab}/*.csv')):
        d=pd.read_csv(f,encoding='cp932',skiprows=3,dtype=str)
        fs.append(d)
    d=pd.concat(fs,ignore_index=True)
    d.columns=[c.strip() for c in d.columns]
    d=d[d['町丁字コード']!='-'].copy()
    d['KEY_CODE']=d['市区町村コード'].str.zfill(5)+d['町丁字コード']
    return d
def num(s): return pd.to_numeric(s.replace({'-':'0'}),errors='coerce')
t6=rd('t06_1',True); t7=rd('t07_1',False); t9=rd('t09',True); t12=rd('t12',True)
c6=t6.columns
a=t6[t6.iloc[:,1]=='総数'].set_index('KEY_CODE')
k=t6[t6.iloc[:,1]=='うち18歳未満世帯員のいる一般世帯'].set_index('KEY_CODE')
out=pd.DataFrame(index=a.index)
out['level']=a['地域階層レベル'].astype(int); out['secret']=a['秘匿処理'].fillna(''); out['secret_to']=a['秘匿先情報'].fillna('')
out['city']=a['市区町村コード']
out['H']=num(a['総数'])
out['eld_couple']=num(a['（再掲）夫65歳以上，妻60歳以上の夫婦のみの世帯'])
out['eld_single']=num(a['（再掲）65歳以上の単独世帯'])
lone=num(k['核家族世帯'])-num(k['うち夫婦のみの世帯'])-num(k['うち夫婦と子供から成る世帯'])
out['lone_parent']=lone.reindex(out.index)
h7=t7.set_index('KEY_CODE')
out['H_dwell']=num(h7['住宅に住む一般世帯']).reindex(out.index)
out['rent']=(num(h7['公営・都市再生機構・公社の借家'])+num(h7['民営の借家'])+num(h7['給与住宅'])+num(h7['間借り'])).reindex(out.index)
l9=t9[t9.iloc[:,1]=='総数'].set_index('KEY_CODE')
out['LF']=num(l9['労働力人口']).reindex(out.index)
o=t12[t12.iloc[:,1]=='総数'].set_index('KEY_CODE')
out['EMP']=num(o['0_総数']).reindex(out.index)
out['svc_sales']=(num(o['D_販売従事者'])+num(o['E_サービス職業従事者'])).reindex(out.index)
out['agri']=num(o['G_農林漁業従事者']).reindex(out.index)
out['blue']=(num(o['H_生産工程従事者'])+num(o['I_輸送・機械運転従事者'])+num(o['J_建設・採掘従事者'])+num(o['K_運搬・清掃・包装等従事者'])).reindex(out.index)
out['unemp']=(out.LF-out.EMP).clip(lower=0)
W={'eld_couple':3.66,'eld_single':8.73,'lone_parent':17.31,'rent':2.37,'svc_sales':3.40,'agri':4.57,'blue':4.79,'unemp':14.97}
den={'eld_couple':'H','eld_single':'H','lone_parent':'H','rent':'H_dwell','svc_sales':'LF','agri':'LF','blue':'LF','unemp':'LF'}
for v in W: out['p_'+v]=out[v]/out[den[v]].where(out[den[v]]>0)
out['ADI_raw']=sum(W[v]*out['p_'+v] for v in W)
ok=out[[ 'p_'+v for v in W]].notna().all(axis=1) & (out.H>=10) & (out.LF>=10)
out.loc[~ok,'ADI_raw']=np.nan
out['ADI_raw_own']=out['ADI_raw']
# missing values (suppressed or very small areas) are filled from the larger unit (4-digit town code), then from the municipality
par=out.index.str[:9]
pv=out['ADI_raw'].groupby(out.index).first()
fill_par=pd.Series(par,index=out.index).map(pv)
out['ADI_src']=np.where(out.ADI_raw.notna(),'own',np.where(fill_par.notna(),'parent',''))
out['ADI_raw']=out['ADI_raw'].fillna(fill_par)
# municipality values (from the row for the whole municipality)
def cityrow(t,typ=None):
    d=t if typ is None else t[t.iloc[:,1]==typ]
    return d
out.to_csv(B+'data/adi_small_area_2020.csv')
print(out.level.value_counts().sort_index()); print(out.ADI_src.value_counts()); print(out[['p_'+v for v in W]+['ADI_raw']].describe().T.round(3))
