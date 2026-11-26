"""Make the data files of the fine-grid web map (docs/grid.html): nighttime road traffic noise and near-road NO2
as population-weighted means of the receptors in grid cells of about 0.4 km, with two coarser levels for small scales.

Usage: python src/webmap/make_grid_data.py RECEPTORS_NO2.parquet RECEPTOR_BURDEN.parquet OUT_DIR
  RECEPTORS_NO2.parquet     receptors_no2_noise_2015_2023.parquet of the NO2 study (columns rid, lon, lat, road_2022)
  RECEPTOR_BURDEN.parquet   receptor_burden_v2.parquet (columns rid, pop_w, Lnight_2021)
Writes OUT_DIR/grid_lnight.js and OUT_DIR/grid_no2.js. Each file holds, for every level, the non-empty blocks of
256 x 256 cells as 8-bit PNG images in base64. Codes: 0 no residents, 1 few residents, 3-255 the value on a linear
scale (Lnight 40-75 dB, near-road NO2 0-20 ug/m3). The receptor file gives NO2 in ppb; it is converted here (1 ppb = 1.88 ug/m3). The cells are the same as in Fig. 1c-f of the paper (0.004 degrees).
"""
import sys, io, json, base64
import numpy as np, pandas as pd, pyarrow.parquet as pq
from PIL import Image

f_no2, f_burden, out = sys.argv[1:4]
BOX = (122.9, 146.0, 24.0, 45.6)                     # lon0, lon1, lat0, lat1 (same as the national map of Fig. 1)
LEVELS = [(0.004, 3), (0.016, 10), (0.064, 20)]      # cell size in degrees, residents below which a cell is shown in gray
BLOCK = 256
SCALE = {'lnight': (40.0, 75.0), 'no2': (0.0, 20.0)}
UG_PER_PPB = 1.88                                   # NO2 at 25 degrees C

N = pq.read_table(f_no2, columns=['rid', 'lon', 'lat', 'road_2022']).to_pandas()
R = pd.read_parquet(f_burden, columns=['rid', 'pop_w', 'Lnight_2021']).merge(N, on='rid'); del N
w = R.pop_w.values.astype('float64'); L = R.Lnight_2021.values.astype('float64'); okL = ~np.isnan(L)
E = np.where(okL, 10 ** (np.nan_to_num(L) / 10), 0.0); road = R.road_2022.values.astype('float64') * UG_PER_PPB
x0, x1, y0, y1 = BOX


def code(v, lo, hi):
    return (np.clip(np.round((v - lo) / (hi - lo) * 252), 0, 252) + 3).astype(np.uint8)


data = {k: {'box': BOX, 'block': BLOCK, 'scale': SCALE[k], 'levels': []} for k in SCALE}
for res, minpop in LEVELS:
    nx = int(np.ceil((x1 - x0) / res)); ny = int(np.ceil((y1 - y0) / res))
    ix = np.floor((R.lon.values - x0) / res).astype(np.int64); iy = np.floor((R.lat.values - y0) / res).astype(np.int64)
    m = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    cell, inv = np.unique(iy[m] * nx + ix[m], return_inverse=True); n = len(cell)
    W = np.bincount(inv, w[m], n); Lw = np.bincount(inv, (w * okL)[m], n); Es = np.bincount(inv, (w * E)[m], n); Ns = np.bincount(inv, (w * road)[m], n)
    val = {'lnight': np.where(Lw >= minpop, code(10 * np.log10(np.maximum(Es, 1e-9) / np.maximum(Lw, 1e-9)), *SCALE['lnight']), 1).astype(np.uint8),
           'no2': np.where(W >= minpop, code(Ns / np.maximum(W, 1e-9), *SCALE['no2']), 1).astype(np.uint8)}
    cx, cy = (cell % nx) // BLOCK, (cell // nx) // BLOCK; ncx = int(np.ceil(nx / BLOCK)); key = cy * ncx + cx
    order = np.argsort(key, kind='stable'); ks, start = np.unique(key[order], return_index=True); end = np.append(start[1:], len(order))
    for k in SCALE:
        blocks = {}
        for kk, a, b in zip(ks, start, end):
            j = order[a:b]; img = np.zeros((BLOCK, BLOCK), np.uint8)
            img[(cell[j] // nx) % BLOCK, (cell[j] % nx) % BLOCK] = val[k][j]            # row 0 is the southern edge of the block
            buf = io.BytesIO(); Image.fromarray(img, mode='L').save(buf, 'PNG', optimize=True)
            blocks[int(kk)] = base64.b64encode(buf.getvalue()).decode()
        data[k]['levels'].append({'res': res, 'nx': nx, 'ny': ny, 'ncx': ncx, 'minpop': minpop, 'blocks': blocks})
    print(f'res {res}: {nx} x {ny}, cells with residents {n:,}, shown in color: Lnight {int((val["lnight"] > 1).sum()):,}, NO2 {int((val["no2"] > 1).sum()):,}, blocks {len(ks)}')
for k in SCALE:
    s = 'window.GRID_DATA=window.GRID_DATA||{};window.GRID_DATA.%s=%s;\n' % (k, json.dumps(data[k], separators=(',', ':')))
    open(f'{out}/grid_{k}.js', 'w').write(s); print(f'grid_{k}.js {len(s) / 1e6:.2f} MB')
