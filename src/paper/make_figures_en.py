"""English Figures 1-4 (EPS and PNG) of the paper.

Inputs: data/ under F_ROOT (default: current directory).
Outputs: figures/Fig{1..4}.{eps,png}
Font: Arial if it is installed, otherwise Liberation Sans (metric-compatible with Arial; the submitted figures use it). Set F_FONT to choose another font. Text stays as text
(EPS Type 42, no outlines), so labels can be edited. EPS has no transparency, so maps are composed as
opaque RGB images.
usage: python3 make_figures_en.py [1 2 3 4]
"""
import os, sys, numpy as np, pandas as pd
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

B = os.environ.get('F_ROOT','.') + '/'
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(B, 'figures'); os.makedirs(OUT, exist_ok=True)
def pick_font():
    from matplotlib import font_manager
    for f in filter(None, [os.environ.get('F_FONT'), 'Arial', 'Liberation Sans', 'Helvetica']):
        try: font_manager.findfont(f, fallback_to_default=False); return f
        except ValueError: pass
    return 'DejaVu Sans'
FONT = pick_font(); print('font:', FONT)
mpl.rcParams.update({'font.family': FONT, 'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'legend.fontsize': 6.5, 'axes.edgecolor': '#c3c2b7', 'axes.linewidth': 0.6, 'xtick.color': '#52514e', 'ytick.color': '#52514e', 'axes.labelcolor': '#0b0b0b',
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False, 'axes.spines.right': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.dpi': 600, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.unicode_minus': True, 'legend.frameon': False,
    'mathtext.fontset': 'custom', 'mathtext.rm': FONT, 'mathtext.it': FONT + ':italic', 'mathtext.bf': FONT + ':bold', 'mathtext.default': 'regular'})
INK = '#0b0b0b'; INK2 = '#52514e'; MUTED = '#898781'; GRID = '#e1e0d9'
BLUE = '#2a78d6'; ORANGE = '#eb6834'; AQUA = '#1baf7a'; VIOLET = '#6250d6'; GRAY = '#898781'; LBLUE = '#86b6ef'
BLUES = LinearSegmentedColormap.from_list('b', ['#cde2fb', '#86b6ef', '#3987e5', '#256abf', '#184f95', '#0d366b'])
ORANGES = LinearSegmentedColormap.from_list('o', ['#fde6da', '#f6b393', '#ef8a5d', '#eb6834', '#b9461b', '#7a2c0f'])
MM = 1 / 25.4
META_PNG = {'Software': None}


def save(fig, name):
    fig.savefig(f'{OUT}/{name}.png', metadata=META_PNG)
    fig.savefig(f'{OUT}/{name}.eps', metadata={'Creator': ''})
    plt.close(fig)


# ---------------------------------------------------------------- Fig. 1
def fig1():
    def load(n): return dict(np.load(B + f'data/fig1_grid_{n}.npz'))

    def show(ax, d, var, cmap, vmin, vmax, minpop, lim=None):
        x0, x1, y0, y1 = d['box']; nx, ny = int(d['nx']), int(d['ny'])
        if var == 'L': v = 10 * np.log10(np.maximum(d['Lenergy'], 1e-9)); ok = d['Lw'] >= minpop
        else: v = d['NO2road']; ok = d['W'] >= minpop
        rgb = np.ones((ny * nx, 3))                                  # white: no population
        rgb[d['W'] > 0] = mpl.colors.to_rgb('#e1e0d9')               # light gray: few residents
        rgb[ok] = cmap(Normalize(vmin, vmax, clip=True)(v[ok]))[:, :3]
        rgb = rgb.reshape(ny, nx, 3); res = float(d['res']); asp = 1 / np.cos(np.deg2rad((y0 + y1) / 2))
        if lim is not None:                                          # crop to the shown window (keeps the EPS small)
            i0 = max(int((lim[0] - x0) / res), 0); i1 = min(int(np.ceil((lim[1] - x0) / res)), nx); j0 = max(int((lim[2] - y0) / res), 0); j1 = min(int(np.ceil((lim[3] - y0) / res)), ny)
            rgb = rgb[j0:j1, i0:i1]; x0, x1, y0, y1 = x0 + i0 * res, x0 + i1 * res, y0 + j0 * res, y0 + j1 * res
        ax.imshow(rgb, origin='lower', extent=(x0, x1, y0, y1), interpolation='none')
        if lim is not None: ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3])
        ax.set_aspect(asp); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
        return ScalarMappable(Normalize(vmin, vmax), cmap)

    fig = plt.figure(figsize=(180 * MM, 205 * MM))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.25, 1, 1], width_ratios=[1, 1, 0.95], hspace=0.12, wspace=0.08, left=0.01, right=0.99, top=0.97, bottom=0.05)
    nat = load('nat'); ka = load('kanto'); ks = load('kansai')
    UG = 'µg m$^{-3}$'
    specs = [('L', BLUES, 40, 60, 'Nighttime noise $\\mathit{L}_{night}$ (dB)', 'nighttime noise'), ('N', ORANGES, 0, 4, 'Near-road NO$_2$ (' + UG + ', FY2022)', 'near-road NO$_2$')]
    for j, (var, cm, a, b, lab, short) in enumerate(specs):
        ax = fig.add_subplot(gs[0, j]); im = show(ax, nat, var, cm, a, b, 20, (127.5, 146, 30.5, 45.6))
        ax.set_title(r'$\bf{' + 'ab'[j] + '}$   ' + lab, loc='left', color=INK, pad=2, fontsize=7.5)
        ia = ax.inset_axes([0.0, 0.62, 0.30, 0.22]); show(ia, nat, var, cm, a, b, 20, (122.9, 131.5, 24.0, 28.0))
        for s in ia.spines.values(): s.set_visible(True); s.set_color('#c3c2b7'); s.set_linewidth(0.5)
        ax.add_patch(plt.Rectangle((139.2, 35.3), 1.0, 0.7, fill=False, ec=INK, lw=0.6)); ax.add_patch(plt.Rectangle((135.05, 34.4), 0.8, 0.6, fill=False, ec=INK, lw=0.6))
        cax = ax.inset_axes([0.5, 0.09, 0.45, 0.03]); cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='max'); cb.outline.set_visible(False)
        cb.ax.tick_params(labelsize=6.5, length=2); cb.set_label('dB' if var == 'L' else UG, fontsize=6.5, labelpad=1)
        ax2 = fig.add_subplot(gs[1, j]); show(ax2, ka, var, cm, a, b, 3)
        ax2.set_title(r'$\bf{' + 'cd'[j] + '}$   Tokyo area (' + short + ')', loc='left', color=INK, pad=2, fontsize=7.5)
        ax3 = fig.add_subplot(gs[2, j]); show(ax3, ks, var, cm, a, b, 3)
        ax3.set_title(r'$\bf{' + 'ef'[j] + '}$   Osaka area (' + short + ')', loc='left', color=INK, pad=2, fontsize=7.5)
    ax = fig.add_subplot(gs[:, 2]); ax.axis('off')
    ax.text(0.0, 0.985, r'$\bf{g}$   National health burden (main analysis)', transform=ax.transAxes, va='bottom', fontsize=7.5, color=INK)
    items = [('Excess high sleep disturbance in adults', '1.52 million', '95% interval 1.07–2.09 million\n(nighttime noise; counterfactual 40 dB)', BLUE),
             ('Deaths attributable to near-road NO$_2$', '6,460 per year', 'Age 30 years and over, RR 1.02 per 10 ' + UG + '\n(about 12,700 with RR 1.04)', ORANGE)]
    for i, (t, v, s, c) in enumerate(items):
        y = 0.86 - 0.24 * i; ax.add_patch(plt.Rectangle((0.05, y - 0.155), 0.012, 0.17, transform=ax.transAxes, color=c, lw=0))
        ax.text(0.09, y, t, transform=ax.transAxes, fontsize=7.2, color=INK2, va='top'); ax.text(0.09, y - 0.035, v, transform=ax.transAxes, fontsize=15, color=INK, va='top')
        ax.text(0.09, y - 0.092, s, transform=ax.transAxes, fontsize=6.5, color=INK2, va='top', linespacing=1.35)
    bx = ax.inset_axes([0.37, 0.14, 0.56, 0.17])
    vals = [15183, 63843]; labs = ['Sleep disturbance\n(YLD)', 'NO$_2$ deaths\n(YLL)']
    bx.barh([1, 0], vals, color=[BLUE, ORANGE], height=0.55)   # YLD 1,518,265 x 0.010; YLL 6,460 x 9.88
    for yv, vv in zip([1, 0], vals): bx.text(vv + 1200, yv, f'{vv / 1e3:.1f}', va='center', fontsize=7, color=INK)
    bx.set_yticks([1, 0]); bx.set_yticklabels(labs, fontsize=6.5); bx.set_xlim(0, 80000); bx.set_xticks([0, 20000, 40000, 60000, 80000]); bx.set_xticklabels(['0', '20', '40', '60', '80'])
    bx.set_xlabel('DALYs (thousand per year)'); bx.xaxis.grid(True, color=GRID, lw=0.5); bx.set_axisbelow(True)
    ax.text(0.05, 0.345, 'DALYs (limited to the two health\noutcomes assessed in this study)', transform=ax.transAxes, fontsize=6.5, color=INK2, va='bottom')
    ax.text(0.05, 0.062, 'YLD with a disability weight of 0.010 (WHO, 2024);\n106 thousand with 0.07 (WHO, 2011).', transform=ax.transAxes, fontsize=6, color=MUTED, va='top')
    save(fig, 'Fig1')


# ---------------------------------------------------------------- Fig. 2
def fig2():
    T = pd.read_csv(B + 'data/fig2_v8_by_quintile.csv'); T = T[T.unit == 'small_area']
    def g(var, st, col): X = T[(T['var'] == var) & (T.stratum == st)].sort_values('q'); return X.q.values, X[col].values
    fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 64 * MM)); plt.subplots_adjust(left=0.06, right=0.985, top=0.80, bottom=0.2, wspace=0.40)
    xt = ['Q1\nlowest', 'Q2', 'Q3', 'Q4', 'Q5\nhighest']
    def lines(ax, series, ylab, title, ymax, dy=None, xmax=7.3):
        for i, (lab, (x, y), c, ls, mk) in enumerate(series):
            ax.plot(x, y, color=c, ls=ls, lw=1.6, marker=mk, ms=4.5, mec='white', mew=0.8, zorder=3)
            ax.text(5.15, y[-1] + (dy[i] if dy else 0), f'{lab}\nratio {y[-1] / y[0]:.2f}', color=INK2, fontsize=6.3, va='center', ha='left', linespacing=1.15)
        ax.set_xticks(range(1, 6)); ax.set_xticklabels(xt); ax.set_xlim(0.7, xmax); ax.set_ylim(0, ymax); ax.set_yticks(np.arange(0, ymax + 0.1, 5 if ymax > 10 else 2))
        ax.set_ylabel(ylab); ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
        ax.set_title(title, loc='left', color=INK, fontsize=7.5, pad=14)
    lines(axs[0], [('Renter\nshare', g('p_rent', 'all', 'HSD_per1000adult'), INK, '-', 'o'), ('ADI', g('ADI_raw', 'all', 'HSD_per1000adult'), MUTED, '--', 's')],
          'Per 1,000 adults', r'$\bf{a}$   Excess high sleep disturbance', 20, [0.3, -0.3])
    lines(axs[1], [('Renter\nshare', g('p_rent', 'all', 'NO2D_std_per100k30'), INK, '-', 'o'), ('ADI', g('ADI_raw', 'all', 'NO2D_std_per100k30'), MUTED, '--', 's')],
          'Per 100,000 aged ≥30 years\n(common national death rate)', r'$\bf{b}$   Deaths attributable to near-road NO$_2$', 10, [0.45, -0.45])
    lines(axs[2], [('Mainly\ndetached', g('p_rent', 'detached_major', 'HSD_per1000adult'), VIOLET, '-', 'o'), ('Mainly\napartment', g('p_rent', 'apartment_major', 'HSD_per1000adult'), '#c98500', '-', 'D'),
                   ('All', g('p_rent', 'all', 'HSD_per1000adult'), MUTED, ':', '')], 'Per 1,000 adults', r'$\bf{c}$   Sleep disturbance by housing type', 25, [-1.2, 2.3, -0.9], 7.6)
    for ax in axs[:2]: ax.text(0, 1.02, 'Within-municipality quintiles of small areas', transform=ax.transAxes, fontsize=5.8, color=MUTED)
    axs[2].text(0, 1.02, 'Small areas stratified by main housing type', transform=axs[2].transAxes, fontsize=5.8, color=MUTED)
    for ax in axs: ax.set_xlabel('Within-municipality quintile of renter share or ADI' if ax is not axs[2] else 'Within-municipality quintile of renter share', fontsize=6.5)
    save(fig, 'Fig2')


# ---------------------------------------------------------------- Fig. 3
def fig3():
    D = pd.read_csv(B + 'data/fig3_curves.csv'); ALL = D[D.method == 'all'].red.iloc[0]
    fig, axs = plt.subplots(1, 2, figsize=(180 * MM, 72 * MM), sharey=True); plt.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.14, wspace=0.12)
    for ax, cost, xcol, lab, t in [(axs[0], 'len', 'len_pct', 'road length', r'$\bf{a}$   Budget measured as road length'), (axs[1], 'lane', 'lane_pct', 'lane length', r'$\bf{b}$   Budget measured as lane length')]:
        ax.axhline(ALL, color=INK, lw=0.8, ls=(0, (1, 1.5))); ax.text(99, ALL + 0.6, f'All sections treated {ALL:.1f}%', ha='right', va='bottom', fontsize=6.5, color=INK2)
        def line(sub, c, ls, name, lw=1.6):
            x = np.r_[0, sub[xcol].values]; y = np.r_[0, sub.red.values]
            if x[-1] < 100: x = np.r_[x, 100]; y = np.r_[y, y[-1]]
            ax.plot(x, y, color=c, ls=ls, lw=lw, label=name)
            return x, y
        g0 = line(D[(D.method == 'greedy') & (D.cost == cost) & (D.gamma == 0)], BLUE, '-', 'Efficient allocation (γ = 0)')
        line(D[(D.method == 'greedy') & (D.cost == cost) & (D.gamma == 16)], ORANGE, (0, (4, 1.5)), 'Equity-weighted allocation (γ = 16)')
        line(D[(D.method == 'rent_only') & (D.cost == cost)], AQUA, '-', 'Renter-share-only allocation', 1.3)
        Rr = D[D.method.str.startswith('random') & (D.cost == cost)]; Rm = Rr.groupby(Rr.groupby('method').cumcount()).agg({xcol: 'mean', 'red': 'mean'})
        line(Rm, GRAY, '-', 'Random allocation (mean of 100 orders)', 1.2)
        y10 = np.interp(10, *g0); ax.plot([10, 10], [0, y10], color=MUTED, lw=0.6); ax.scatter([10], [y10], s=18, color=BLUE, ec='white', lw=0.6, zorder=5)
        ax.annotate(f'10% of {lab}: {y10:.1f}%\n({y10 / ALL * 100:.0f}% of the all-sections benefit)', (10, y10), xytext=(30, -10), textcoords='offset points', fontsize=6.8, color=INK,
                    arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.5))
        ax.set_xlim(0, 100); ax.set_ylim(0, 33); ax.set_xlabel(f'Share of {lab} treated (%)'); ax.set_xticks([0, 10, 20, 40, 60, 80, 100])
        ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True); ax.set_title(t, loc='left', fontsize=7.5, color=INK, pad=4)
    axs[0].set_ylabel('Reduction in excess high sleep disturbance (%)')
    axs[0].text(97, 4, 'The curves for γ = 0 and γ = 16 almost overlap\n(22.1% and 21.5% at 10% of road length).\nFig. 4 shows the difference and the disparity.', ha='right', va='bottom', fontsize=6.5, color=INK2, linespacing=1.35)
    h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc='upper center', ncol=4, bbox_to_anchor=(0.525, 0.995), fontsize=6.5, handlelength=2.3, columnspacing=1.2)
    save(fig, 'Fig3')


# ---------------------------------------------------------------- Fig. 4
def fig4():
    A = pd.read_csv(B + 'data/v5_alloc_opt.csv'); A = A[A.budget == 10]
    V7 = pd.read_csv(B + 'data/v7_checks.csv').set_index('name'); V8 = pd.read_csv(B + 'data/v8_asj_checks.csv').set_index('name')
    V8b = pd.read_csv(B + 'data/v8b_asj_years.csv')
    RND = pd.read_csv(B + 'data/v10_random100.csv')
    def rnd(meas, col):
        v = RND[RND.measure == meas][{'ratio': 'ratio', 'diff': 'gap', 'reduction_pct': 'red'}[col]]; return v.mean(), v.quantile(0.025), v.quantile(0.975)
    base_ratio, base_gap = 1.7073, 7.2163
    G = [0, 0.5, 1, 2, 4, 8, 16]
    def Gs(meas): return [g for g in G if ((A.measure == meas) & (A.method == f'opt_g{g:g}')).any()]
    def curve(meas, col): return np.array([A[(A.measure == meas) & (A.method == f'opt_g{g:g}')][col].iloc[0] for g in Gs(meas)])
    def pt(meas, meth, col): return A[(A.measure == meas) & (A.method == meth)][col].iloc[0]
    fig = plt.figure(figsize=(180 * MM, 178 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.6, 1], left=0.075, right=0.98, top=0.955, bottom=0.06, wspace=0.24, hspace=0.36)
    axa = fig.add_subplot(gs[0, 0]); axb = fig.add_subplot(gs[0, 1]); axc = fig.add_subplot(gs[1, 0]); axd = fig.add_subplot(gs[1, 1])
    def frame(ax, title, xl, xlim, ylim, xbase, left=False):
        ax.axvline(xbase, color=MUTED, lw=0.7, ls=(0, (3, 2)))
        ax.annotate('Before\nmitigation', (xbase, ylim[1]), xytext=(-3 if left else 3, -1), textcoords='offset points', color=MUTED, fontsize=6.5, va='top', ha='right' if left else 'left', linespacing=1.1)
        ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_xlabel(xl); ax.set_ylabel('Reduction in excess high sleep disturbance (%)')
        ax.yaxis.grid(True, color=GRID, lw=0.5); ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
        ax.set_title(title, loc='left', fontsize=7.5, color=INK, pad=4)
    def trade(ax, meas, xcol, allx, ally, big=True, g16left=False, all_left=False):
        x = curve(meas, xcol); y = curve(meas, 'reduction_pct')
        ax.plot(x, y, color=LBLUE, lw=1.4, zorder=3); ax.scatter(x[1:-1], y[1:-1], s=10, color=LBLUE, ec='white', lw=0.5, zorder=4)
        s = 40 if big else 22; fs = 7 if big else 6.5
        ax.scatter([x[0]], [y[0]], s=s, color=BLUE, ec='white', lw=0.7, zorder=5); ax.scatter([x[-1]], [y[-1]], s=s, color=ORANGE, ec='white', lw=0.7, zorder=5)
        ax.annotate('Efficient (γ = 0)', (x[0], y[0]), xytext=(4, -11), textcoords='offset points', fontsize=fs, color=INK, ha='left')
        ax.annotate('Equity-weighted (γ = 16)', (x[-1], y[-1]), xytext=(-7, -3) if g16left else (-4, -11), textcoords='offset points', fontsize=fs, color=INK,
                    ha='right' if (g16left or not big) else 'center', va='center' if g16left else 'baseline')
        pts = [('Renter-share-only', pt(meas, 'old_rent', xcol), pt(meas, 'old_rent', 'reduction_pct'), AQUA, 'D'),
               ('Random (mean of 100 orders)', rnd(meas, xcol)[0], rnd(meas, 'reduction_pct')[0], GRAY, 's'),
               ('All sections treated\n(reference)', allx, ally, INK, '^')]
        rx_, ry_ = rnd(meas, xcol), rnd(meas, 'reduction_pct')
        ax.errorbar([rx_[0]], [ry_[0]], xerr=[[rx_[0] - rx_[1]], [rx_[2] - rx_[0]]], yerr=[[ry_[0] - ry_[1]], [ry_[2] - ry_[0]]], fmt='none', ecolor=GRAY, elinewidth=0.6, capsize=0, zorder=4)
        for lab, xx, yy, c, mk in pts:
            ax.scatter([xx], [yy], s=s * 0.8, marker=mk, color=c, ec='white', lw=0.6, zorder=5)
            isall = lab.startswith('All'); israd = lab.startswith('Random')
            L = (all_left and isall) or israd
            Aflag = all_left and isall; U = (israd and not big) or (Aflag and not big)
            if Aflag:
                ax.annotate('All sections\ntreated\n(reference)', (xx, yy), xytext=(6, 0), ha='left', va='center', textcoords='offset points', fontsize=6.8 if big else 6.5, color=INK2, linespacing=1.1)
                continue
            ax.annotate(lab, (xx, yy), xytext=(-5, 3) if U else ((-5, -3) if L else (5, -3)), ha='right' if L else 'left',
                        textcoords='offset points', fontsize=6.8 if big else 6.5, color=INK2, va='bottom' if U else 'top', linespacing=1.1)
        return x, y
    # a
    frame(axa, r'$\bf{a}$   Low-noise pavement: sleep benefit and disparity ratio', 'Disparity ratio (Q5/Q1 of renter share, rate per adult)', (1.40, 1.80), (0, 35), base_ratio)
    x, y = trade(axa, 'pav', 'ratio', 1.7085, 29.5867, all_left=True)
    axa.annotate('', xy=(x[-1], y[-1] + 1.0), xytext=(x[0], y[0] + 1.0), arrowprops=dict(arrowstyle='-|>', color=INK, lw=1.1, mutation_scale=9, connectionstyle='arc3,rad=0.12'))
    axa.text(0.03, 0.975, 'Central estimate: efficient → γ = 16', transform=axa.transAxes, fontsize=7, color=INK, va='top')
    axa.text(0.03, 0.895, 'Sleep benefit', transform=axa.transAxes, fontsize=7, color=INK2, va='top')
    axa.text(0.03, 0.855, '−2.7%', transform=axa.transAxes, fontsize=15, color=INK, va='top')
    axa.text(0.03, 0.765, '22.1% → 21.5%', transform=axa.transAxes, fontsize=7, color=INK2, va='top')
    axa.text(0.33, 0.895, 'Disparity ratio', transform=axa.transAxes, fontsize=7, color=INK2, va='top')
    axa.text(0.33, 0.855, '−0.11', transform=axa.transAxes, fontsize=15, color=INK, va='top')
    axa.text(0.33, 0.765, '1.63 → 1.52', transform=axa.transAxes, fontsize=7, color=INK2, va='top')
    axa.text(0.97, 0.11, 'Uncertainty (sections fixed to the main selection;\n2.5th–97.5th percentiles of within-run differences)\nSleep benefit −2.6 to −2.8%\nDisparity ratio −0.15 to −0.07',
             transform=axa.transAxes, fontsize=6.3, color=INK2, va='bottom', ha='right', multialignment='left', linespacing=1.35, bbox=dict(boxstyle='round,pad=0.45', fc='white', ec='#c3c2b7', lw=0.6))
    # b
    frame(axb, r'$\bf{b}$   Low-noise pavement: sleep benefit and absolute gap', 'Absolute gap (Q5 − Q1, per 1,000 adults)', (2.0, 7.7), (0, 35), base_gap, left=True)
    x, y = trade(axb, 'pav', 'diff', 5.1072, 29.5867, g16left=True)
    g16x, g16y = x[-1], y[-1]; rx, ry = pt('pav', 'old_rent', 'diff'), pt('pav', 'old_rent', 'reduction_pct')
    axb.add_patch(plt.Rectangle((g16x, 0), 7.7 - g16x, g16y, color='#f0efec', lw=0, zorder=1))
    axb.text(5.0, 9.5, 'Larger gap and smaller\nreduction than γ = 16', fontsize=6.5, color=MUTED, ha='left', va='top')
    axb.annotate('', xy=(g16x + 0.03, g16y - 1.2), xytext=(rx - 0.02, ry + 1.0), arrowprops=dict(arrowstyle='-|>', color=INK, lw=1.0, mutation_scale=9))
    axb.text(0.03, 0.975, 'Central estimate: γ = 16 and renter-share-only', transform=axb.transAxes, fontsize=7, color=INK, va='top')
    axb.text(0.03, 0.905, 'Absolute gap (per 1,000 adults)', transform=axb.transAxes, fontsize=7, color=INK2, va='top')
    axb.text(0.03, 0.875, '4.5 < 4.7', transform=axb.transAxes, fontsize=15, color=INK, va='top')
    axb.text(0.03, 0.79, 'Sleep benefit', transform=axb.transAxes, fontsize=7, color=INK2, va='top')
    axb.text(0.03, 0.76, '21.5% > 14.1%', transform=axb.transAxes, fontsize=15, color=INK, va='top')
    axb.text(0.03, 0.675, 'Smaller gap and larger benefit with γ = 16', transform=axb.transAxes, fontsize=7, color=INK2, va='top')
    gx = V7.loc['absgap_only']
    axb.scatter([gx.rent_gap], [gx.red], s=32, marker='o', color=VIOLET, ec='white', lw=0.6, zorder=5)
    axb.annotate('Gap-only', (gx.rent_gap, gx.red), xytext=(-4, -7), textcoords='offset points', fontsize=6.8, color=INK2, va='top')
    # c
    ax = axc
    rows = [('ASJ, year 0', V8.loc['asj_y0_reopt_g0_len'], V8.loc['asj_y0_reopt_g16_len']),
            ('Uniform 3 dB (main analysis)', V7.loc['uniform3_reopt_g0_len'], V7.loc['uniform3_reopt_g16_len']),
            ('ASJ, after 3 years', V8.loc['asj_y3_reopt_g0_len'], V8.loc['asj_y3_reopt_g16_len']),
            ('Rolling noise only', V7.loc['rolling3_reopt_g0_len'], V7.loc['rolling3_reopt_g16_len'])]
    b0 = V8b[(V8b.sel == 'uniform_sel_g0') & (V8b.method == 'health_avg')].iloc[0]; b16 = V8b[(V8b.sel == 'uniform_sel_g16') & (V8b.method == 'health_avg')].iloc[0]
    rows.append(('ASJ, 10-year average', b0, b16))
    frame(ax, r'$\bf{c}$   By pavement performance assumption (γ = 0 → 16)', 'Disparity ratio (Q5/Q1)', (1.40, 1.80), (0, 35), base_ratio)
    DY = {'ASJ, after 3 years': 1.9, 'Rolling noise only': 0, 'ASJ, 10-year average': -1.9}
    for lab, r0, r16 in rows:
        c = BLUE if 'main' in lab else INK2
        ax.annotate('', xy=(r16.rent_ratio, r16.red), xytext=(r0.rent_ratio, r0.red), arrowprops=dict(arrowstyle='-|>', color=c, lw=1.2, mutation_scale=7, shrinkA=2, shrinkB=2))
        ax.scatter([r0.rent_ratio], [r0.red], s=12, color=c, ec='white', lw=0.5, zorder=4)
        ax.text(r0.rent_ratio + 0.006, r0.red + DY.get(lab, 0), lab, fontsize=6.5, color=INK if c == BLUE else INK2, va='center')
    # d
    frame(axd, r'$\bf{d}$   50% reduction in nighttime heavy-vehicle traffic', 'Disparity ratio (Q5/Q1)', (1.40, 1.80), (0, 35), base_ratio)
    trade(axd, 'hv', 'ratio', 1.7172, 16.8302, big=False, all_left=True)
    save(fig, 'Fig4')


if __name__ == '__main__':
    which = sys.argv[1:] or ['1', '2', '3', '4']
    for k in which: {'1': fig1, '2': fig2, '3': fig3, '4': fig4}[k](); print('Fig', k, 'done')
