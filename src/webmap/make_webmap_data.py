"""Make docs/muni.geojson for the web map: excess high sleep disturbance and the main allocations by municipality.

Usage: python src/webmap/make_webmap_data.py data/muni_allocation.csv BOUNDARIES.geojson docs/muni.geojson
BOUNDARIES.geojson: simplified municipal polygons with the properties c (5-digit code), n (name) and pop, for example
docs/muni.geojson of https://github.com/m138/japan-night-road-noise (2020 census municipalities).

Properties of each municipality (lists hold the efficient allocation and the equity-weighted allocation with gamma = 16):
  ad    adults (20 years and over)                    n0   adults with excess high sleep disturbance (HSD) before mitigation
  hsd   excess HSD per 1,000 adults                    gap  rate in the highest renter-share quintile minus the lowest (per 1,000 adults)
  sel   share of the length of ordinary roads selected (%)      red  reduction in excess HSD (%)
  gapa  gap after mitigation (per 1,000 adults)        av   avoided cases       km, skm   length of ordinary roads and selected length (km)
A value is null where it is not defined (no adults in one of the two quintiles, no excess HSD, no section).
"""
import sys, json
import numpy as np, pandas as pd

csv, bnd, out = sys.argv[1:4]
T = pd.read_csv(csv, dtype={'municipality': str}).set_index('municipality')
AL = ('efficient', 'equity_gamma16')


def r(v, d=2):
    return None if not np.isfinite(v) else round(float(v), d)


def props(s):
    div = lambda a, b: a / b if b > 0 else np.nan
    both = s.adults_q1 > 0 and s.adults_q5 > 0
    rate = lambda k, q: div(s[f'excess_hsd_{q}'] - (s[f'avoided_{k}_{q}'] if k else 0), s[f'adults_{q}']) * 1e3
    gap = lambda k: rate(k, 'q5') - rate(k, 'q1') if both else np.nan
    return {'ad': int(round(s.adults)), 'n0': r(s.excess_hsd, 1), 'hsd': r(div(s.excess_hsd, s.adults) * 1e3), 'gap': r(gap(None)),
            'sel': [r(100 * div(s[f'selected_m_{k}'], s.road_m), 1) for k in AL], 'red': [r(100 * div(s[f'avoided_{k}'], s.excess_hsd), 1) for k in AL],
            'gapa': [r(gap(k)) for k in AL], 'av': [r(s[f'avoided_{k}'], 1) for k in AL],
            'km': r(s.road_m / 1e3, 1), 'skm': [r(s[f'selected_m_{k}'] / 1e3, 1) for k in AL]}


g = json.load(open(bnd, encoding='utf8'))
keep = []
for f in g['features']:
    p = f['properties']
    if p['c'] not in T.index:
        continue
    f['properties'] = {'c': p['c'], 'n': p['n'], 'pop': p['pop'], **props(T.loc[p['c']])}
    keep.append(f)
g['features'] = keep
g['national'] = props(T.sum())
open(out, 'w', encoding='utf8').write(json.dumps(g, ensure_ascii=False, separators=(',', ':')))
print(len(keep), 'municipalities;', g['national'])
