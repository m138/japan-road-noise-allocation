"""Map for the README: share of the length of ordinary roads selected for low-noise pavement, by municipality.

Efficient allocation (gamma = 0), uniform 3 dB, 10% of the national length of ordinary roads.

Usage: python src/webmap/make_readme_map.py data/muni_allocation.csv BOUNDARIES.geojson docs/images/map_selected_share.png
BOUNDARIES.geojson: simplified municipal polygons with the property c (5-digit code), for example docs/muni.geojson of
https://github.com/m138/japan-night-road-noise (2020 census municipalities).
"""
import sys, pandas as pd, geopandas as gpd, shapely
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

CSV, BND, OUT = sys.argv[1:4]
have = {f.name for f in font_manager.fontManager.ttflist}
FONT = next(f for f in ('Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans') if f in have)
plt.rcParams.update({'font.family': FONT, 'font.size': 9})
T = pd.read_csv(CSV, dtype={'municipality': str}).set_index('municipality')
T['share'] = 100 * T.selected_m_efficient / T.road_m
G = gpd.read_file(BND, columns=['c', 'geometry']).to_crs(6677).merge(T[['share']], left_on='c', right_index=True, how='left')

bins = [-1, 0, 10, 25, 50, 75, 101]          # none, >0-10, 10-25, 25-50, 50-75, 75-100
cols = ['#eeede8', '#cde2fb', '#86b6ef', '#3987e5', '#1f5fae', '#0d366b']
labs = ['None', '<10%', '10–25%', '25–50%', '50–75%', '≥75%']
G['cls'] = pd.cut(G.share.fillna(0), bins, labels=False)      # intervals closed on the right: 0 -> "None"
isl = G.c.str.startswith('47') | G.c.isin([f'13{c}' for c in ['361', '362', '363', '364', '381', '382', '401', '402', '421']])
fig, ax = plt.subplots(figsize=(6.5, 7))


def draw(a, D):
    for i, c in enumerate(cols):
        D[D.cls == i].plot(ax=a, color=c, lw=0)
    a.set_axis_off(); a.set_aspect('equal')


draw(ax, G[~isl])
xmin, ymin, xmax, ymax = G[~isl].total_bounds
ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
ia = ax.inset_axes([0.04, 0.66, 0.2, 0.26]); draw(ia, G[G.c.str.startswith('47')])
bb = gpd.GeoSeries([shapely.box(127.6, 26.05, 128.35, 26.9)], crs=4326).to_crs(6677).total_bounds
ia.set_xlim(bb[0], bb[2]); ia.set_ylim(bb[1], bb[3]); ia.set_axis_on(); ia.set_xticks([]); ia.set_yticks([]); ia.set_xlabel(''); ia.set_ylabel('')
for sp in ia.spines.values():
    sp.set_color('#c3c2b7'); sp.set_visible(True)
ia.set_title('Okinawa Island', fontsize=7, loc='left')
h = [plt.Rectangle((0, 0), 1, 1, color=c) for c in cols]
ax.legend(h, labs, title='Share of the length of ordinary roads\nselected for low-noise pavement\n(efficient allocation of 10% of the\nnational length, municipality)',
          loc='lower right', frameon=False, fontsize=7, title_fontsize=8)
fig.tight_layout(); fig.savefig(OUT, dpi=300, metadata={'Software': None}); plt.close()
n = (T.selected_m_efficient > 0).sum()
print('font:', FONT, '| municipalities with selected sections', n, 'of', len(T), '| national share %.2f%%' % (100 * T.selected_m_efficient.sum() / T.road_m.sum()))
