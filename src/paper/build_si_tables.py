"""Build the Supplementary Tables (markdown) of the paper from the result files in data/.
Output: si_tables.md (included by build.py into the Supplementary Information).
Ratios, gaps and relative losses are calculated from unrounded values."""
import os, re, numpy as np, pandas as pd
B = os.environ.get('F_ROOT','.') + '/data/'; HERE = os.path.dirname(os.path.abspath(__file__))
out = []; CHECK = {}


def f(x, d=1): return '–' if x is None or (isinstance(x, float) and np.isnan(x)) else f'{x:,.{d}f}'.replace('-', '−')
def table(head, rows): out.append('| ' + ' | '.join(head) + ' |\n|' + '---|' * len(head) + '\n' + '\n'.join('| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows) + '\n')
def title(n, t, note): out.append(f'**Supplementary Table {n} | {t}** {note}\n')
def after(t): out.append(t + '\n')


Q = pd.read_csv(B + 'v9_quintile_rates.csv').set_index(['measure', 'name'])
MC = pd.read_csv(B + 'v5_mc_summary.csv').drop_duplicates('stat').set_index('stat')
DR = pd.concat([pd.read_csv(B + 'v5_mc_draws_a.csv'), pd.read_csv(B + 'v5_mc_draws_b.csv')]); DR = DR[DR.draw >= 1] if (DR.draw >= 1).sum() >= 1000 else DR[DR.draw >= 0]
RND = pd.read_csv(B + 'v10_random100.csv')
V7 = pd.read_csv(B + 'v7_checks.csv'); V8 = pd.read_csv(B + 'v8_asj_checks.csv').set_index('name'); V8b = pd.read_csv(B + 'v8b_asj_years.csv').set_index(['sel', 'method'])
G6 = pd.read_csv(B + 'v6_greedy_checks.csv')
ALLRED = Q.loc[('pav', 'all'), 'red_pct']; AV0 = Q.loc[('pav', 'gamma0'), 'avoided']
CHECK['n_mc'] = len(DR)

# ---------------------------------------------------------------- Table 1
NAMES = [('none', 'Before mitigation'), ('gamma0', 'Efficient allocation (γ = 0)'), ('gamma0.5', 'Equity-weighted, γ = 0.5'), ('gamma1', 'Equity-weighted, γ = 1'), ('gamma2', 'Equity-weighted, γ = 2'),
         ('gamma4', 'Equity-weighted, γ = 4'), ('gamma8', 'Equity-weighted, γ = 8'), ('gamma16', 'Equity-weighted, γ = 16'), ('rent_only', 'Renter-share-only allocation'), ('gap_only', 'Gap-only allocation'), ('all', 'All sections treated')]
title(1, 'Allocation methods at 10% of road length (low-noise pavement, uniform 3 dB reduction).',
      'Excess high sleep disturbance (HSD) in adults. Q1–Q5 are the within-municipality quintiles of renter share (Q5 highest); adults in Q1–Q5 are 20.17, 20.76, 20.73, 20.85 and 22.56 million. The disparity ratio is the Q5/Q1 ratio of rates and the absolute gap is Q5 − Q1 (per 1,000 adults). Rates are shown to two decimals; ratios, gaps and relative losses were calculated from unrounded values. A section that crosses the 1% step is included as a whole, so the selected length slightly exceeds 10%.')
after('**a, Selected length, reduction and disparity.**')
rows = []
for k, lab in NAMES:
    r = Q.loc[('pav', k)]
    loss = '–' if k in ('none', 'all') else f(100 * (1 - r.avoided / AV0), 1)
    rows.append([lab, f(r.len_pct, 2), f(r.n_sections, 0), f(r.avoided / 1e3, 1), f(r.red_pct, 1), loss, f(r.ratio, 2), f(r.gap, 1)])
p = RND[RND.measure == 'pav']
rows.insert(-1, ['Random allocation (mean of 100 orders)', f(p.len_pct.mean(), 2), '–', f(p.red.mean() / 100 * Q.loc[('pav', 'none'), 'HSD0'] / 1e3, 1), f(p.red.mean(), 1), f(100 * (1 - p.red.mean() / Q.loc[('pav', 'gamma0'), 'red_pct']), 1), f(p.ratio.mean(), 2), f(p.gap.mean(), 1)])
table(['Allocation', 'Selected length (%)', 'Sections', 'Avoided excess HSD (thousand adults)', 'Reduction (%)', 'Loss relative to efficient allocation (%)', 'Disparity ratio', 'Absolute gap'], rows)
after('**b, Excess HSD per 1,000 adults after mitigation, by renter-share quintile.**')
table(['Allocation', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5'], [[lab] + [f(Q.loc[('pav', k), f'rate_q{i}'], 2) for i in range(1, 6)] for k, lab in NAMES])
after('**c, Avoided excess HSD (thousand adults), by renter-share quintile.**')
table(['Allocation', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q1–Q3', 'Q4–Q5'], [[lab] + [f(Q.loc[('pav', k), f'avoided_q{i}'] / 1e3, 1) for i in range(1, 6)] +
      [f(sum(Q.loc[('pav', k), f'avoided_q{i}'] for i in (1, 2, 3)) / 1e3, 1), f(sum(Q.loc[('pav', k), f'avoided_q{i}'] for i in (4, 5)) / 1e3, 1)] for k, lab in NAMES[1:]])
d13 = sum(Q.loc[('pav', 'gamma16'), f'avoided_q{i}'] - Q.loc[('pav', 'gamma0'), f'avoided_q{i}'] for i in (1, 2, 3)); d45 = sum(Q.loc[('pav', 'gamma16'), f'avoided_q{i}'] - Q.loc[('pav', 'gamma0'), f'avoided_q{i}'] for i in (4, 5))
CHECK['shift_q13_q45_total'] = (round(d13), round(d45), round(Q.loc[('pav', 'gamma16'), 'avoided'] - AV0))
after('**d, Differences between allocations within the same Monte Carlo run** (1,000 runs; sections fixed to the main selection; central estimate and 2.5th–97.5th percentiles of the within-run differences).')
def pct(s): return np.percentile(s, [2.5, 97.5])
def cen(a, b, col): return MC.loc[f'{col}_pav_{a}_b10', 'central'] - MC.loc[f'{col}_pav_{b}_b10', 'central']
rows = []
for a, b, lab in [('opt_g16', 'opt_g0', 'γ = 16 minus efficient'), ('opt_g16', 'old_rent', 'γ = 16 minus renter-share-only'), ('old_rent', 'opt_g0', 'Renter-share-only minus efficient')]:
    dred = DR[f'red_pav_{a}_b10'] - DR[f'red_pav_{b}_b10']; drat = DR[f'ratio_pav_{a}_b10'] - DR[f'ratio_pav_{b}_b10']; dgap = (DR[f'diff_pav_{a}_b10'] - DR[f'diff_pav_{b}_b10']) * 1e3
    rel = 100 * (1 - DR[f'red_pav_{a}_b10'] / DR[f'red_pav_{b}_b10'])
    crel = 100 * (1 - MC.loc[f'red_pav_{a}_b10', 'central'] / MC.loc[f'red_pav_{b}_b10', 'central'])
    def cell(c, s, d): lo, hi = pct(s); return f'{f(c, d)} ({f(lo, d)} to {f(hi, d)})'
    rows.append([lab, cell(cen(a, b, 'red'), dred, 2), cell(crel, rel, 1) if b == 'opt_g0' else '–', cell(cen(a, b, 'ratio'), drat, 2), cell(cen(a, b, 'diff') * 1e3, dgap, 2), f'{(drat < 0).mean() * 100:.0f}', f'{(dgap < 0).mean() * 100:.0f}'])
    CHECK[lab] = rows[-1][1:]
table(['Comparison', 'Reduction (percentage points)', 'Relative loss in sleep benefit (%)', 'Disparity ratio', 'Absolute gap (per 1,000 adults)', 'Runs with a smaller ratio (%)', 'Runs with a smaller gap (%)'], rows)

# ---------------------------------------------------------------- Table 2
title(2, 'Pavement performance assumptions (10% of road length).',
      'ASJ: equations for porous asphalt on ordinary roads in ASJ RTN-Model 2018 (light vehicles 5.7 − 7.3 lg(1 + y) dB, heavy vehicles 3.9 − 3.6 lg(1 + y) dB, y = years after paving). For year 0, after 3 years and rolling noise only, sections were reselected with each reduction. For the 10-year average, the sections of the main analysis were fixed, avoided cases, the disparity ratio and the absolute gap were calculated for each year y = 0–9, and the ten values were averaged. In the ASJ equation the reduction for light vehicles is negative from y = 6; the main calculation keeps the negative values, and the next row sets the lower limit to 0 dB.')
v7 = V7.set_index('name')
def r2(lab, a, g0, g16, sel):
    return [lab, f(a.red, 1), f(g0.red, 1), f(100 * g0.avoided / a.avoided, 1), f(g16.red, 1), f(100 * (1 - g16.avoided / g0.avoided), 1), f'{f(g0.rent_ratio, 2)} → {f(g16.rent_ratio, 2)}', f(g16.rent_ratio - g0.rent_ratio, 2), f'{f(g0.rent_gap, 1)} → {f(g16.rent_gap, 1)}', sel]
rows = [r2('Uniform 3 dB (main analysis)', v7.loc['uniform3_all'], v7.loc['uniform3_reopt_g0_len'], v7.loc['uniform3_reopt_g16_len'], 'Reselected'),
        r2('ASJ, year 0', V8.loc['asj_y0_all'], V8.loc['asj_y0_reopt_g0_len'], V8.loc['asj_y0_reopt_g16_len'], 'Reselected'),
        r2('ASJ, after 3 years', V8.loc['asj_y3_all'], V8.loc['asj_y3_reopt_g0_len'], V8.loc['asj_y3_reopt_g16_len'], 'Reselected'),
        r2('ASJ, 10-year average', V8b.loc[('all', 'health_avg')], V8b.loc[('uniform_sel_g0', 'health_avg')], V8b.loc[('uniform_sel_g16', 'health_avg')], 'Fixed (main)'),
        r2('ASJ, 10-year average, lower limit 0 dB', V8b.loc[('all', 'health_avg_floor0')], V8b.loc[('uniform_sel_g0', 'health_avg_floor0')], V8b.loc[('uniform_sel_g16', 'health_avg_floor0')], 'Fixed (main)'),
        r2('Rolling noise only, 3 dB', v7.loc['rolling3_all'], v7.loc['rolling3_reopt_g0_len'], v7.loc['rolling3_reopt_g16_len'], 'Reselected')]
table(['Assumption', 'All sections: reduction (%)', 'Efficient: reduction (%)', 'Share of all-sections benefit (%)', 'γ = 16: reduction (%)', 'Relative loss with γ = 16 (%)', 'Disparity ratio, γ = 0 → 16', 'Change in ratio', 'Absolute gap, γ = 0 → 16', 'Sections'], rows)
CHECK['t2_share'] = [r[3] for r in rows]; CHECK['t2_loss'] = [r[5] for r in rows]

# ---------------------------------------------------------------- Table 3
title(3, 'Budget unit and continuous treatment units.',
      'The budget is 10% of the total of each unit over all ordinary roads of the census (94,112 sections). Continuous units: adjacent sections of the same route (first seven digits of the census section number; end points within 0.0005°) were joined in order, and a unit was closed each time the joined length reached 2 km (or 5 km); the remaining sections were added to the last unit, and chains shorter than 2 km (5 km) were excluded. Eligible length and sections: 96.3% and 82,950 sections for 2 km units (42,028 units); 85.5% and 69,054 sections for 5 km units (20,624 units). The share of the all-sections benefit uses two denominators: all ordinary-road sections (29.6%) and all eligible sections (27.9% for 2 km, 23.6% for 5 km). With the 2 km and 5 km units, the treated length is 10.4% and 11.7% of the eligible length.')
A = G6[G6.part == 'A_cost'].set_index(['cost', 'gamma']); U = V7[V7.part == '3_unit'].set_index('name'); UD = pd.read_csv(B + 'v10_units_denominator.csv').set_index('name')
def r3(lab, g0, g16, red='red', ratio='ratio', gap=None, elig=None):
    row = [lab, f(g0['len_pct'], 1), f(g0[red], 1), f(100 * g0[red] / ALLRED, 1), '–' if elig is None else f(100 * g0[red] / elig, 1), f(g16[red], 1), f(100 * (1 - g16.avoided / g0.avoided), 1), f'{f(g0[ratio], 2)} → {f(g16[ratio], 2)}']
    return row
rows = [r3('Road length', A.loc[('len', 0.0)], A.loc[('len', 16.0)]), r3('Lane length', A.loc[('lane', 0.0)], A.loc[('lane', 16.0)]), r3('Paved area', A.loc[('area', 0.0)], A.loc[('area', 16.0)]),
        r3('Road length, units of 2 km or more', U.loc['units2km_g0_len'], U.loc['units2km_g16_len'], ratio='rent_ratio', elig=UD.loc['eligible_2km', 'red']),
        r3('Road length, units of 5 km or more', U.loc['units5km_g0_len'], U.loc['units5km_g16_len'], ratio='rent_ratio', elig=UD.loc['eligible_5km', 'red'])]
table(['Budget (10%)', 'Road length selected, γ = 0 (%)', 'Efficient: reduction (%)', 'Share of benefit, all sections (%)', 'Share of benefit, eligible sections (%)', 'γ = 16: reduction (%)', 'Relative loss with γ = 16 (%)', 'Disparity ratio, γ = 0 → 16'], rows)
CHECK['t3'] = [(r[0], r[3], r[4], r[6], r[7]) for r in rows]
st = V7[V7.part == '3_stats'].set_index('name')
after(f"Sections selected by the efficient allocation at 10% of road length without the unit constraint lay in {st.loc['free_g0_len', 'n_muni']:.0f} municipalities; {st.loc['free_g0_len', 'share_len_runs_ge2km']:.0f}% of the selected length was in continuous runs of 2 km or more and {st.loc['free_g0_len', 'share_len_runs_ge5km']:.0f}% in runs of 5 km or more.")

# ---------------------------------------------------------------- Table 4
title(4, 'Definitions of equity (10% of road length, low-noise pavement, uniform 3 dB).',
      'Weights were applied to the within-municipality quintiles of renter share or of the areal deprivation index (ADI). The gap-only allocation maximizes the fall in the Q5 rate minus the fall in the Q1 rate of renter share. Ratios and gaps are Q5/Q1 and Q5 − Q1 of excess HSD per 1,000 adults. For the ADI, a ratio below 1 means a lower rate in the most deprived quintile.')
E = V7[V7.part == '2_eq'].set_index('name')
E.loc['random', ['red', 'rent_ratio', 'rent_gap']] = [p.red.mean(), p.ratio.mean(), p.gap.mean()]
lab4 = [('eff_g0', 'Efficient allocation (γ = 0)'), ('rent_g4', 'Renter-share weight, γ = 4'), ('rent_g16', 'Renter-share weight, γ = 16'), ('adi_g4', 'ADI weight, γ = 4'), ('adi_g16', 'ADI weight, γ = 16'),
        ('absgap_only', 'Gap-only allocation'), ('rent_only_v4', 'Renter-share-only allocation'), ('all', 'All sections treated')]
rows = [['Before mitigation', '0.0', f(E.rent_ratio0.iloc[0], 2), f(E.rent_gap0.iloc[0], 1), f(E.adi_ratio0.iloc[0], 2), f(E.adi_gap0.iloc[0], 1)]] + [[l, f(E.loc[k, 'red'], 1), f(E.loc[k, 'rent_ratio'], 2), f(E.loc[k, 'rent_gap'], 1), f(E.loc[k, 'adi_ratio'], 2), f(E.loc[k, 'adi_gap'], 1)] for k, l in lab4]
table(['Allocation', 'Reduction (%)', 'Renter share: ratio', 'Renter share: gap', 'ADI: ratio', 'ADI: gap'], rows)
after('**b, Variables of the ADI and replacements made for the 2020 small-area tabulation.** Weights are those of the original index (ref. [[nakaya2011]]).')
table(['Variable (share)', 'Weight', 'Source in the 2020 small-area tabulation'], [
    ['Elderly couple households', '3.66', 'As in the original index'], ['Elderly single households', '8.73', 'As in the original index'],
    ['Single-parent households', '17.31', 'Replaced: nuclear-family households with a member under 18 minus couple-only and couple-with-children households with a member under 18 (includes single-father households; 2.1% of private households nationally). The original index uses single-mother households.'],
    ['Households in rented housing', '2.37', 'As in the original index'], ['Sales and service workers', '3.40', 'As in the original index'], ['Agriculture, forestry and fishery workers', '4.57', 'As in the original index'],
    ['Production, transport, construction and carrying workers', '4.79', 'As in the original index'], ['Unemployed persons', '14.97', 'Replaced: labor force minus employed persons']])

# ---------------------------------------------------------------- Table 5
title(5, 'Robustness of the selected sections and checks of the sequential selection (efficient allocation, 10% of road length).',
      'Loss: relative difference in the reduction between the sections of the main analysis evaluated under the changed exposure and the sections reselected under the changed exposure. Overlap: share of the selected length in common with the main analysis.')
Bp = G6[(G6.part == 'B_robust') & (G6.measure == 'pav')]; B20 = Bp[Bp.draw != 'Jcurve']; BJ = Bp[Bp.draw == 'Jcurve'].iloc[0]
SN = pd.read_csv(B + 'v6_greedy_secnoise.csv'); LS = pd.read_csv(B + 'v6_greedy_local.csv').set_index('method')
def rng(s, d=1): return f'{f(s.min(), d)}–{f(s.max(), d)}'
rows = [['Common exposure errors (prefecture SD 0.97 dB, municipality SD 1.37 dB), 20 draws', rng(B20.regret_pct, 2), rng(B20.overlap_len_pct, 1)],
        ['Independent error for each section, SD 2 dB, 10 draws', rng(SN[SN.sd_db == 2].regret_pct, 2), rng(SN[SN.sd_db == 2].overlap_len_pct, 1)],
        ['Independent error for each section, SD 3 dB, 10 draws', rng(SN[SN.sd_db == 3].regret_pct, 2), rng(SN[SN.sd_db == 3].overlap_len_pct, 1)],
        ['Japanese exposure–response curve', f(BJ.regret_pct, 2), f(BJ.overlap_len_pct, 1)]]
table(['Change in exposure or exposure–response', 'Loss (%)', 'Overlap (%)'], rows)
CHECK['t5'] = rows
after('**b, Checks of the sequential selection.**')
table(['Check', 'Reduction (%)', 'Selected length (%)', 'Result'], [
    ['Steps of 1% of road length (main analysis)', f(LS.loc['greedy_step0.01', 'red'], 2), f(LS.loc['greedy_step0.01', 'len_pct'], 2), '–'],
    ['Steps of 0.25%', f(LS.loc['greedy_step0.0025', 'red'], 2), f(LS.loc['greedy_step0.0025', 'len_pct'], 2), f"+{LS.loc['greedy_step0.0025', 'red'] - LS.loc['greedy_step0.01', 'red']:.2f} percentage points"],
    ['Swap of 0.5% of length between selected and unselected sections', f(LS.loc['swap_chunk0.005', 'red'], 2), f(LS.loc['swap_chunk0.005', 'len_pct'], 2), 'No improvement'],
    ['Swap of 0.1% of length', f(LS.loc['swap_chunk0.001', 'red'], 2), f(LS.loc['swap_chunk0.001', 'len_pct'], 2), 'No improvement'],
    ['Swap of 0.02% of length', f(LS.loc['swap_chunk0.0002', 'red'], 2), f(LS.loc['swap_chunk0.0002', 'len_pct'], 2), 'No improvement']])

# ---------------------------------------------------------------- Table 6
title(6, 'Monte Carlo settings and results, variance explained, receptors above 65 dB, and other exposure–response settings.',
      'The 95% intervals are the 2.5th–97.5th percentiles of 1,000 runs of the assumed uncertainty distributions; they are not confidence intervals based on sampling. Sections were fixed to the main selection.')
after('**a, Settings.**')
table(['Source of uncertainty', 'Distribution in each run'], [
    ['Exposure error common to an area', 'Ordinary roads: sum of normal errors with SD 0.97 dB for each prefecture and SD 1.37 dB for each municipality. Motorways: normal error with SD 2.87 dB for each prefecture. Estimated from the residuals of the comparison with measurements.'],
    ['Slope of the exposure–response relationship', 'Log odds ratio per 10 dB drawn from a normal distribution with mean ln(2.13) and the standard error from the 95% confidence interval 1.82–2.48 (ref. [[basner2018]]); a logistic curve through the WHO value at 40 dB with this slope rescales the excess at each receptor.'],
    ['Pavement noise reduction', 'Triangular distribution with minimum 1.5 dB, mode 3 dB and maximum 4.5 dB.']])
after('**b, Results and share of variance explained.** R² from separate regressions of each indicator on the drawn slope and its square, on the drawn pavement reduction and its square, and on both with their product (joint model). The single-factor R² values and the unexplained share come from different regressions and do not sum to 1.')
VS = pd.read_csv(B + 'v8_mc_variance_share.csv')
lab6 = ['Excess HSD in adults, national (million)', 'Disparity ratio before mitigation', 'Reduction with the efficient allocation (%)', 'Share of the all-sections benefit (%)', 'Relative loss with γ = 16 (%)', 'Change in disparity ratio (γ = 16 minus γ = 0)',
        'Heavy vehicles: relative loss with γ = 16 (%)', 'Heavy vehicles: change in disparity ratio']
cen6 = [Q.loc[('pav', 'none'), 'HSD0'] / 1e6, Q.loc[('pav', 'none'), 'ratio0'], Q.loc[('pav', 'gamma0'), 'red_pct'], 100 * AV0 / Q.loc[('pav', 'all'), 'avoided'], 100 * (1 - Q.loc[('pav', 'gamma16'), 'avoided'] / AV0),
        Q.loc[('pav', 'gamma16'), 'ratio'] - Q.loc[('pav', 'gamma0'), 'ratio'], 100 * (1 - Q.loc[('hv', 'gamma16'), 'avoided'] / Q.loc[('hv', 'gamma0'), 'avoided']), Q.loc[('hv', 'gamma16'), 'ratio'] - Q.loc[('hv', 'gamma0'), 'ratio']]
sc = [1e-6, 1, 1, 100, 100, 1, 100, 1]; dd = [2, 2, 1, 1, 1, 2, 1, 2]
rows = [[lab6[i], f(cen6[i], dd[i]), f'{f(VS.p2_5[i] * sc[i], dd[i])} to {f(VS.p97_5[i] * sc[i], dd[i])}', f(VS.slope[i], 2), f(VS.pavement[i], 2), f(VS.both[i], 2), f(VS.unexplained[i], 2)] for i in range(8)]
table(['Indicator', 'Central estimate', '2.5th–97.5th percentiles', 'R², slope only', 'R², pavement reduction only', 'R², joint model', 'Unexplained by joint model'], rows)
CHECK['t6b'] = [(r[0][:25], r[1], r[2]) for r in rows]
after('**c, Receptors above 65 dB, Japanese exposure–response curve and definition of night.** The WHO curve was derived for 40–65 dB. In the main analysis it is extrapolated above 65 dB (0.83% of adults, about 0.87 million; maximum 74.4 dB; 6.3% of excess HSD). Sections are those of the main analysis.')
C = pd.read_csv(B + 'v9_cap65.csv').set_index('name'); HS = pd.read_csv(B + 'v3_hsd_scenarios.csv').set_index(['measure', 'scen'])
def r6(lab, tot, r0, r16, ra, rat0, rat16): return [lab, f(tot / 1e6, 2), f(r0, 1), f(100 * r0 / ra, 1), f(r16, 1), f(100 * (1 - r16 / r0), 1), f'{f(rat0, 2)} → {f(rat16, 2)}']
rows = [r6('Main analysis (WHO curve, extrapolated above 65 dB)', Q.loc[('pav', 'none'), 'HSD0'], Q.loc[('pav', 'gamma0'), 'red_pct'], Q.loc[('pav', 'gamma16'), 'red_pct'], ALLRED, Q.loc[('pav', 'gamma0'), 'ratio'], Q.loc[('pav', 'gamma16'), 'ratio']),
        r6('WHO curve capped at the 65 dB value', C.loc['none', 'HSD0'], C.loc['gamma0', 'red_pct'], C.loc['gamma16', 'red_pct'], C.loc['all', 'red_pct'], C.loc['gamma0', 'ratio'], C.loc['gamma16', 'ratio']),
        r6('Japanese curve (ref. [[morinaga2025]])', MC.loc['H_base', 'jcurve'], MC.loc['red_pav_opt_g0_b10', 'jcurve'], MC.loc['red_pav_opt_g16_b10', 'jcurve'], MC.loc['red_pav_all_b100', 'jcurve'], MC.loc['ratio_pav_opt_g0_b10', 'jcurve'], MC.loc['ratio_pav_opt_g16_b10', 'jcurve']),
        ['Night defined as 23:00–07:00 (ordinary roads about +1.1 dB)', f(HS.loc[('att40', 'n2307'), 'thousand'] / 1e3, 2), '–', '–', '–', '–', '–']]
table(['Setting', 'Excess HSD in adults (million)', 'Efficient: reduction (%)', 'Share of all-sections benefit (%)', 'γ = 16: reduction (%)', 'Relative loss with γ = 16 (%)', 'Disparity ratio, γ = 0 → 16'], rows)
CHECK['t6c'] = rows; CHECK['cap_total_change_pct'] = 100 * (C.loc['none', 'HSD0'] / Q.loc[('pav', 'none'), 'HSD0'] - 1)

# ---------------------------------------------------------------- Table 7
title(7, 'Reduction of nighttime heavy-vehicle traffic by 50% (10% of road length) and the co-benefit for NO2.',
      'No shift of traffic to other hours or routes is assumed. Deaths attributable to near-road NO2 (age 30 years and over, relative risk 1.02 per 10 µg m⁻³) were recalculated after the fall in nighttime heavy-vehicle NOx emissions was applied to the near-road term of the NO2 model. Sections were selected by the reduction in excess HSD only.')
N = pd.read_csv(B + 'v5_no2_heavy.csv').set_index('scenario'); A5 = pd.read_csv(B + 'v5_alloc_opt.csv'); A5 = A5[(A5.budget == 10) & (A5.measure == 'hv')].set_index('method')
h = RND[RND.measure == 'hv']; allh = Q.loc[('hv', 'all')]
def r7(lab, red, ratio, gap, nk):
    return [lab, f(red, 1), f(100 * red / allh.red_pct, 1), f(ratio, 2), f(gap, 1), '–' if nk is None else f(N.loc[nk, 'deaths_avoided'], 0), '–' if nk is None else f(100 * N.loc[nk, 'deaths_avoided'] / N.loc['all_b100', 'deaths_avoided'], 0)]
rows = [r7('Efficient allocation (γ = 0)', Q.loc[('hv', 'gamma0'), 'red_pct'], Q.loc[('hv', 'gamma0'), 'ratio'], Q.loc[('hv', 'gamma0'), 'gap'], 'opt_g0_b10'),
        r7('Equity-weighted, γ = 1', Q.loc[('hv', 'gamma1'), 'red_pct'], Q.loc[('hv', 'gamma1'), 'ratio'], Q.loc[('hv', 'gamma1'), 'gap'], 'opt_g1_b10'),
        r7('Equity-weighted, γ = 4', Q.loc[('hv', 'gamma4'), 'red_pct'], Q.loc[('hv', 'gamma4'), 'ratio'], Q.loc[('hv', 'gamma4'), 'gap'], 'opt_g4_b10'),
        r7('Equity-weighted, γ = 16', Q.loc[('hv', 'gamma16'), 'red_pct'], Q.loc[('hv', 'gamma16'), 'ratio'], Q.loc[('hv', 'gamma16'), 'gap'], 'opt_g16_b10'),
        r7('Renter-share-only allocation', A5.loc['old_rent', 'reduction_pct'], A5.loc['old_rent', 'ratio'], A5.loc['old_rent', 'diff'], 'old_rent_b10'),
        r7('Random allocation (mean of 100 orders; NO2 for one order)', h.red.mean(), h.ratio.mean(), h.gap.mean(), 'old_random00_b10'),
        r7('All sections treated', allh.red_pct, allh.ratio, allh.gap, 'all_b100')]
table(['Allocation', 'Reduction in excess HSD (%)', 'Share of all-sections benefit (%)', 'Disparity ratio', 'Absolute gap', 'NO2 deaths avoided per year', 'Share of all-sections NO2 benefit (%)'], rows)
CHECK['t7_loss16'] = 100 * (1 - Q.loc[('hv', 'gamma16'), 'avoided'] / Q.loc[('hv', 'gamma0'), 'avoided'])

# ---------------------------------------------------------------- Table 8
title(8, 'Public input data.', 'All data are open and were used under their terms of use. Exposure estimates at the receptors come from the two earlier studies (D.M., manuscripts under review), which used the first five items.')
table(['Data', 'Content and year', 'Use in this study', 'Source'], [
    ['Population Census, small-area tabulation and boundaries', '2020; 210,877 small areas (cho-cho-aza); Table 3 (age), Table 7-1 (tenure of housing), Table 8-1 (type of building), Tables 6-1, 9 and 12 (households, labor force and occupation)', 'Receptor placement, population weights, renter share, housing type, areal deprivation index', 'Statistics Bureau of Japan, e-Stat (https://www.e-stat.go.jp)'],
    ['Road Traffic Census, section tables and hourly traffic tables', 'Fiscal year 2021; traffic volume by vehicle type, speed limit, number of lanes and carriageway width of each section', 'Noise sources, NOx emissions, lane length and paved area', 'Ministry of Land, Infrastructure, Transport and Tourism (https://www.mlit.go.jp/road/census/r3/)'],
    ['Census section geometry', 'Fiscal year 2021; line geometry of the census sections', 'Source locations, continuous treatment units', 'mlit-road-traffic-census-converter (MIT license); original data from the Ministry of Land, Infrastructure, Transport and Tourism'],
    ['Continuous roadside noise monitoring', 'Fiscal years 2019–2023; nighttime LAeq, coordinates and distance from the carriageway edge (12,648 records; 11,978 on ordinary roads)', 'Bias correction and cross-validation of the noise model; residual analysis by stratum', 'Ministry of the Environment, from the Environmental Observatory of the National Institute for Environmental Studies (https://tenbou.nies.go.jp)'],
    ['Continuous air quality monitoring', 'Fiscal year 2022; annual mean NO2 at monitoring stations', 'NO2 model of the earlier study', 'Ministry of the Environment, from the Environmental Observatory of the National Institute for Environmental Studies'],
    ['Vital Statistics', '2022; deaths by five-year age group', 'Baseline deaths, years of life lost', 'Ministry of Health, Labour and Welfare, e-Stat'],
    ['Abridged Life Table', '2022; life expectancy by age and sex', 'Years of life lost per death', 'Ministry of Health, Labour and Welfare'],
    ['National Land Numerical Information, land use subdivision mesh', '2021; 100 m grid', 'Building-land receptors', 'Ministry of Land, Infrastructure, Transport and Tourism (https://nlftp.mlit.go.jp/ksj/)']])

# ---------------------------------------------------------------- Table 9
title(9, 'Small-area regressions of excess HSD on renter share.',
      'Outcome: excess HSD per 1,000 adults in each small area (n = 208,743). Weighted least squares with population weights and municipality fixed effects; 95% confidence intervals from standard errors clustered by municipality. Renter share and apartment share are entered per 10 percentage points and per unit (0–1), respectively. Road density: length of census ordinary-road sections (km per km²). "Mainly apartment": small areas where detached houses are less than half of main households (47,841 small areas, 39% of the population; apartment buildings are half or more of main households for 95% of this population).')
HR = pd.read_csv(B + 'v8_housing_regression.csv')
TERM = {'rent10': 'Renter share (per 10 points)', 'road_dens': 'Road density (km per km²)', 'log_popdens': 'Log population density', 'p_apt': 'Apartment share (0–1)', 'apt_major': 'Mainly apartment (indicator)',
        'rent10_x_apt': 'Renter share × mainly apartment', 'slope_in_apartment_major': 'Slope of renter share in mainly-apartment areas', 'slope_in_detached_major': 'Slope of renter share in mainly-detached areas'}
MOD = {'M1': 'Model 1: renter share only', 'M2': 'Model 2: + road density and log population density', 'M3': 'Model 3: + apartment share', 'M4': 'Model 4: + mainly-apartment indicator',
       'M5': 'Model 5: renter share × mainly apartment, unadjusted', 'M6': 'Model 6: renter share × mainly apartment, adjusted for road density, log population density and apartment share'}
table(['Model', 'Term', 'Coefficient', '95% confidence interval'], [[MOD[r.model[:2]], TERM[r.term], f(r.coef, 2), f'{f(r.lo, 2)} to {f(r.hi, 2)}'] for r in HR.itertuples()])

# ---------------------------------------------------------------- Table 10
title(10, 'Receptors placed on building land.',
      'Receptors of each small area were placed again, with the same number and population, at random within the grid cells classed as building land in the 100 m land use mesh (2021); small areas without such cells (0.5% of the population) kept the uniform placement. Building land includes buildings other than housing. Fixed: sections selected in the main analysis, evaluated with the building-land receptors. Reselected: sections reselected with the building-land receptors.')
BL = pd.read_csv(B + 'v9_bldg_placement.csv').set_index('name'); M = Q.loc['pav']
def r11(lab, r, ref0=None, allr=None, ov=None):
    return [lab, f(r.HSD0 / 1e6 if lab.endswith('before mitigation') else r.HSD / 1e6, 2), f(r.red_pct, 1), '–' if allr is None else f(100 * r.avoided / allr.avoided, 1), '–' if ref0 is None else f(100 * (1 - r.avoided / ref0.avoided), 1), f(r.ratio, 2), f(r.gap, 1), f(r.rate_q1, 2), f(r.rate_q5, 2), '–' if ov is None else f(ov, 1)]
rows = [r11('Uniform placement (main analysis): before mitigation', M.loc['none']), r11('Uniform: efficient allocation', M.loc['gamma0'], None, M.loc['all']), r11('Uniform: γ = 16', M.loc['gamma16'], M.loc['gamma0'], M.loc['all']),
        r11('Uniform: renter-share-only', M.loc['rent_only'], M.loc['gamma0'], M.loc['all']), r11('Uniform: all sections treated', M.loc['all']),
        r11('Building land: before mitigation', BL.loc['none']), r11('Building land, fixed: efficient allocation', BL.loc['fixed_gamma0'], None, BL.loc['all']), r11('Building land, fixed: γ = 16', BL.loc['fixed_gamma16'], BL.loc['fixed_gamma0'], BL.loc['all']),
        r11('Building land, fixed: renter-share-only', BL.loc['fixed_rent_only'], BL.loc['fixed_gamma0'], BL.loc['all']), r11('Building land: all sections treated', BL.loc['all']),
        r11('Building land, reselected: efficient allocation', BL.loc['reopt_gamma0'], None, BL.loc['all'], BL.loc['reopt_gamma0', 'overlap_with_main']), r11('Building land, reselected: γ = 16', BL.loc['reopt_gamma16'], BL.loc['reopt_gamma0'], BL.loc['all'], BL.loc['reopt_gamma16', 'overlap_with_main'])]
table(['Placement and allocation', 'Excess HSD in adults (million)', 'Reduction (%)', 'Share of all-sections benefit (%)', 'Loss relative to efficient allocation (%)', 'Disparity ratio', 'Absolute gap', 'Q1 rate', 'Q5 rate', 'Overlap with main selection (% of length)'], rows)
CHECK['t11'] = [(r[0][:40], r[1], r[2], r[3], r[4], r[5], r[6]) for r in rows]

# ---------------------------------------------------------------- Table 11
title(11, 'Cross-validation residuals of the noise model at monitoring sites on ordinary roads, by stratum.',
      'Residual: estimate minus measurement (dB), leave-one-prefecture-out cross-validation, 11,978 records (fiscal years 2019–2023). Coordinates and distance from the carriageway edge were joined from the public monitoring data by prefecture, fiscal year, site number and municipality code (all records matched), and each site was assigned to the small area that contains it. Slopes of the residual on housing composition: 0.00 dB per 10 points of municipal apartment share (95% confidence interval −0.08 to 0.08); 0.13 dB per renter-share quintile of the small area (0.05 to 0.21), and 0.16 dB (0.07 to 0.25) with road-type fixed effects; standard errors clustered by municipality to allow for correlation between sites in the same municipality.')
RM = pd.read_csv(B + 'v9_cv_residual_by_muni_housing.csv'); RS = pd.read_csv(B + 'v9_cv_residual_by_smallarea_housing.csv'); RD = pd.read_csv(B + 'v9_cv_residual_by_distance.csv')
ROAD = {'3': 'National highways', '4': 'Principal local roads (prefectural)', '5': 'Principal local roads (designated-city)', '6': 'Prefectural roads', '7': 'Municipal roads of designated cities'}
STR = {'municipal apartment share': 'Apartment share of the municipality', 'municipal rent share': 'Renter share of the municipality', 'road type': 'Road type', 'predicted Lnight (dB)': 'Estimated Lnight (dB)',
       'small-area housing type': 'Housing type of the small area', 'within-municipality rent quintile of small area': 'Renter-share quintile of the small area (within municipality)', 'housing type x rent quintile': 'Housing type × renter-share quintile',
       'distance from roadway edge': 'Distance from the carriageway edge', 'all': 'All records'}
def lev(s, l):
    l = str(l)
    if s == 'road type': return ROAD[l[0]]
    l = l.replace('apartment-major', 'Mainly apartment').replace('detached-major', 'Mainly detached').replace('no households', 'No households').replace('no rent share', 'no renter share').replace('20 m+', '20 m or more').replace('-', '–')
    if l.startswith('['): a, b = l.strip('[)').split(', '); l = f'{float(a) * (100 if "share" in s else 1):g} to <{float(b) * (100 if "share" in s else 1):g}' + ('%' if 'share' in s else '')
    return 'All' if l == 'all' else l
rows = [[STR[r.stratum], lev(r.stratum, r.level), f(r.n, 0), f(r.mean_residual, 2), f(r.rmse, 2)] for r in pd.concat([RM[RM.level != 'all'], RS[RS.level != 'all'], RD, RM[RM.level == 'all']]).itertuples()]
table(['Stratum', 'Level', 'Records', 'Mean residual (dB)', 'RMSE (dB)'], rows)

# ---------------------------------------------------------------- Table 12
title(12, 'Disparity compared within municipalities.',
      'Municipalities (wards for designated cities; 1,895) that have small areas in both Q1 and Q5: 1,844 (99.8% of adults). For each municipality the Q1 and Q5 rates of excess HSD per 1,000 adults and their difference were calculated and then combined with two sets of weights. The ratio is the ratio of the weighted mean Q5 and Q1 rates. "Pooled" rows show the national pooling used in the main analysis. Intervals: 2.5th–97.5th percentiles of 2,000 bootstrap samples of municipalities (weights: adults of the municipality), with differences between allocations calculated in the same sample; they show the variation of the combined values when municipalities are resampled, not the Monte Carlo uncertainty of exposure or pavement effect.')
W = pd.read_csv(B + 'v10_within_muni_gap.csv'); SC = [('none', 'Before mitigation'), ('gamma0', 'Efficient allocation'), ('gamma16', 'γ = 16'), ('rent_only', 'Renter-share-only'), ('all', 'All sections treated')]
WL = [('adults_city', 'Within municipalities, weighted by adults of the municipality'), ('adults_q1q5', 'Within municipalities, weighted by adults in Q1 and Q5'), ('pooled_same_cities', 'Pooled, same 1,844 municipalities'), ('pooled_all(main)', 'Pooled, all municipalities (main analysis)')]
rows = []
for wk, wl in WL:
    for sk, sl in SC:
        r = W[(W.weight == wk) & (W.scenario == sk)].iloc[0]; rows.append([wl, sl, f(r.r1, 2), f(r.r5, 2), f(r.gap, 2), f(r.ratio, 2)])
table(['Aggregation', 'Scenario', 'Q1 rate', 'Q5 rate', 'Absolute gap', 'Ratio'], rows)
bs = W[W.weight == 'bootstrap_adults_city'].set_index('scenario'); cc = W[W.weight == 'city_count'].set_index('scenario')
wa = W[W.weight == 'adults_city'].set_index('scenario')
after('**b, Differences and bootstrap intervals (weights: adults of the municipality).**')
table(['Quantity', 'Estimate', '2.5th–97.5th percentiles'], [
    ['Absolute gap, before mitigation', f(wa.loc['none', 'gap'], 2), f"{f(bs.loc['gap_none', 'gap'], 2)} to {f(bs.loc['gap_none', 'ratio'], 2)}"],
    ['Absolute gap, efficient allocation', f(wa.loc['gamma0', 'gap'], 2), f"{f(bs.loc['gap_g0', 'gap'], 2)} to {f(bs.loc['gap_g0', 'ratio'], 2)}"],
    ['Absolute gap, γ = 16', f(wa.loc['gamma16', 'gap'], 2), f"{f(bs.loc['gap_g16', 'gap'], 2)} to {f(bs.loc['gap_g16', 'ratio'], 2)}"],
    ['Absolute gap, renter-share-only', f(wa.loc['rent_only', 'gap'], 2), f"{f(bs.loc['gap_rent', 'gap'], 2)} to {f(bs.loc['gap_rent', 'ratio'], 2)}"],
    ['Difference in gap, γ = 16 minus efficient', f(wa.loc['gamma16', 'gap'] - wa.loc['gamma0', 'gap'], 2), f"{f(bs.loc['dgap_g16_g0', 'gap'], 2)} to {f(bs.loc['dgap_g16_g0', 'ratio'], 2)}"],
    ['Difference in gap, γ = 16 minus renter-share-only', f(wa.loc['gamma16', 'gap'] - wa.loc['rent_only', 'gap'], 2), f"{f(bs.loc['dgap_g16_rent', 'gap'], 2)} to {f(bs.loc['dgap_g16_rent', 'ratio'], 2)}"],
    ['Difference in ratio, γ = 16 minus efficient', f(wa.loc['gamma16', 'ratio'] - wa.loc['gamma0', 'ratio'], 2), f"{f(bs.loc['dratio_g16_g0', 'gap'], 2)} to {f(bs.loc['dratio_g16_g0', 'ratio'], 2)}"]])
g = cc.loc['gamma16_lt_gamma0_gap']; s = cc.loc['same_gap_and_both_rates_equal']
after(f"Municipalities where the within-municipality absolute gap with γ = 16 was smaller than with the efficient allocation: {g.gap:.0f} ({g.r5:.0f}% of adults); larger: {g.r1:.0f}; the same: {g.ratio:.0f} ({s.r5:.0f}% of adults; differences below 10⁻⁹ per 1,000 adults were counted as the same, and in these municipalities the Q1 and Q5 rates were each identical under the two allocations).")
CHECK['t13'] = (g.gap, g.r5, g.r1, g.ratio)

# ---------------------------------------------------------------- Tables 14, 16
def bias_table(n, fn, pre, ttl, note, extra):
    title(n, ttl, note)
    Z = pd.read_csv(B + fn).set_index('name'); S = Z.loc[pre + '_summary']
    def rr(lab, k, base=None, allk=None):
        r = Z.loc[k]; av = r.HSD0 - r.HSD
        return [lab, f(r.red, 2), '–' if allk is None else f(100 * av / (Z.loc[allk].HSD0 - Z.loc[allk].HSD), 1), '–' if base is None else f(100 * (1 - av / (Z.loc[base].HSD0 - Z.loc[base].HSD)), 2), f(r.ratio, 2), f(r.gap, 1)]
    rows = [['Main model: before mitigation', '0.00', '–', '–', f(Z.loc['main_all', 'ratio0'], 2), f(Z.loc['main_all', 'gap0'], 1)], rr('Main model: efficient allocation', 'main_fixed_g0', None, 'main_all'), rr('Main model: γ = 16', 'main_fixed_g16', 'main_fixed_g0', 'main_all'),
            rr('Main model: renter-share-only', 'main_fixed_rent', 'main_fixed_g0', 'main_all'), rr('Main model: all sections treated', 'main_all'),
            ['Corrected model: before mitigation', '0.00', '–', '–', f(Z.loc[pre + '_all', 'ratio0'], 2), f(Z.loc[pre + '_all', 'gap0'], 1)],
            rr('Corrected model, fixed sections: efficient allocation', pre + '_fixed_g0', None, pre + '_all'), rr('Corrected model, fixed sections: γ = 16', pre + '_fixed_g16', pre + '_fixed_g0', pre + '_all'),
            rr('Corrected model, fixed sections: renter-share-only', pre + '_fixed_rent', pre + '_fixed_g0', pre + '_all'),
            rr('Corrected model, reselected: efficient allocation', pre + '_reopt_g0', None, pre + '_all'), rr('Corrected model, reselected: γ = 16', pre + '_reopt_g16', pre + '_reopt_g0', pre + '_all'), rr('Corrected model: all sections treated', pre + '_all')]
    table(['Model and allocation', 'Reduction (%)', 'Share of all-sections benefit (%)', 'Loss relative to efficient allocation (%)', 'Disparity ratio', 'Absolute gap'], rows)
    after(extra(Z, S)); CHECK[f't{n}'] = [(r[0][:45], r[1], r[2], r[3], r[4], r[5]) for r in rows]
def ex14(Z, S):
    a = Z.loc['roadtype_corrected_fixed_g0']; b = Z.loc['roadtype_corrected_reopt_g0']
    return (f"Excess HSD before mitigation with the corrected model: {Z.loc['roadtype_corrected_all', 'HSD0'] / 1e6:.2f} million adults. Efficient allocation: the fixed sections gave a reduction {b.red - a.red:.2f} percentage points smaller than the reselected sections "
            f"(relative difference {S.regret_fixed_g0:.2f}%). Overlap of the reselected sections with the main selection: {S.overlap_g0:.0f}% of length for the efficient allocation and {S.overlap_g16:.0f}% for γ = 16. "
            "National highways made up 32% of the selected length in the main analysis and 35% after reselection.")
bias_table(13, 'v10_roadtype_bias.csv', 'roadtype_corrected', 'Diagnostic sensitivity analysis: correction of the bias by road type.',
           'The mean cross-validation residual of each road type (Supplementary Table 11) was used as a correction: the energy contributed by sections of that type to each receptor was multiplied by 10^(−residual/10) (national highways +1.1 dB, principal local roads (prefectural) −0.4 dB, principal local roads (designated-city) +0.3 dB, prefectural roads −1.1 dB, municipal roads of designated cities −1.2 dB). Contributions of motorways were not changed. Fixed sections: those selected in the main analysis. Reductions are shown to two decimals; relative differences and losses were calculated from unrounded values.', ex14)

# ---------------------------------------------------------------- Table 14
title(14, 'Random allocation: distribution over 100 random orders (10% of road length).', 'Sections were put in random order (seeds 1000–1099) and selected until the cumulative length reached the budget.')
rows = []
for mk, ml in [('pav', 'Low-noise pavement, uniform 3 dB'), ('hv', 'Nighttime heavy-vehicle traffic −50%')]:
    x = RND[RND.measure == mk]
    for col, cl, d in [('red', 'Reduction in excess HSD (%)', 2), ('ratio', 'Disparity ratio', 3), ('gap', 'Absolute gap (per 1,000 adults)', 2)]:
        v = x[col]; rows.append([ml, cl, f(v.mean(), d), f(v.std(), d), f(v.min(), d), f(np.percentile(v, 2.5), d), f(np.percentile(v, 97.5), d), f(v.max(), d)])
table(['Measure', 'Indicator', 'Mean', 'SD', 'Minimum', '2.5th percentile', '97.5th percentile', 'Maximum'], rows)
x = RND[RND.measure == 'pav'].red; CHECK['t15_eff_over_random_min'] = Q.loc[('pav', 'gamma0'), 'red_pct'] / x.max()
after(f"The reduction with the efficient allocation (low-noise pavement, {Q.loc[('pav', 'gamma0'), 'red_pct']:.1f}%) was {Q.loc[('pav', 'gamma0'), 'red_pct'] / x.max():.1f} to {Q.loc[('pav', 'gamma0'), 'red_pct'] / x.min():.1f} times the reduction with each of the 100 random orders.")

def ex16(Z, S):
    return (f"Overlap of the reselected sections with the main selection: {S.overlap_g0:.1f}% of length for the efficient allocation and {S.overlap_g16:.1f}% for γ = 16.")
bias_table(15, 'v10_rentq_bias.csv', 'rentq_corrected', 'Diagnostic sensitivity analysis: correction of the bias by renter-share quintile of the small area.',
           'The mean cross-validation residual of each within-municipality renter-share quintile of the small area that contains the monitoring site (Supplementary Table 11; Q1 −0.43 dB, Q2 −0.11 dB, Q3 +0.13 dB, Q4 +0.05 dB, Q5 +0.15 dB) was used as a correction to the ordinary-road contribution at receptors of that quintile (for example, +0.4 dB in Q1 and −0.1 dB in Q5). Receptors in small areas without a renter share were not corrected. Fixed sections: those selected in the main analysis.', ex16)

# ---------------------------------------------------------------- Table 16
title(16, 'Simple roadside-population rule compared with the other allocations (10% of road length, low-noise pavement, uniform 3 dB).',
      'Simple rule: sections ranked by the number of adults at receptors within 200 m of the section divided by the section length, and selected in that order until the cumulative length reached the budget. It uses neither the noise contributions nor the exposure–response relationship. The second variant ranks sections by the number of adults within 200 m without dividing by length.')
PD = pd.read_csv(B + 'v10_popdensity_rule.csv').set_index('name')
LP = [('gamma0', 'Efficient allocation (γ = 0)'), ('gamma16', 'Equity-weighted allocation (γ = 16)'), ('adults200m_per_length', 'Simple rule: adults within 200 m per unit length'), ('adults200m_total', 'Adults within 200 m (not divided by length)'), ('rent_only', 'Renter-share-only allocation'), ('all', 'All sections treated')]
table(['Allocation', 'Sections', 'Reduction (%)', 'Share of all-sections benefit (%)', 'Loss relative to efficient allocation (%)', 'Disparity ratio', 'Absolute gap', 'Overlap with efficient allocation (% of length)'],
      [[l, f(PD.loc[k, 'n_sections'], 0), f(PD.loc[k, 'red'], 1), f(PD.loc[k, 'share_all'], 1), '–' if k in ('gamma0', 'all') else f(PD.loc[k, 'loss_vs_g0'], 1), f(PD.loc[k, 'ratio'], 2), f(PD.loc[k, 'gap'], 1), '–' if k == 'all' else f(PD.loc[k, 'overlap_g0_len'], 1)] for k, l in LP])

# ---------------------------------------------------------------- Table 17
title(17, 'Disability-adjusted life years (DALYs) for the two health outcomes assessed.',
      'Central estimates per year. The comparison is limited to the two outcomes of this study and does not cover all health effects of road traffic. YLD: excess HSD in adults × disability weight, treating the estimated state as lasting one year. YLL: deaths × 9.9 years, the mean remaining life expectancy at death for deaths at age 30 years and over (2022 abridged life table and 2022 deaths by five-year age group, both sexes combined).')
hsd = Q.loc[('pav', 'none'), 'HSD0']; D = 6460.0; DJ = MC.loc['H_base', 'jcurve']
table(['Outcome and assumption', 'Adults with excess HSD or deaths per year', 'DALYs per year'], [
    ['Excess HSD in adults, WHO curve, disability weight 0.010 (ref. [[whodw2024]])', f(hsd, 0), f(hsd * 0.010, 0)],
    ['Same, disability weight 0.006 to 0.015 (95% uncertainty interval of the weight)', f(hsd, 0), f'{f(hsd * 0.006, 0)} to {f(hsd * 0.015, 0)}'],
    ['Same, disability weight 0.07 used in an earlier assessment (ref. [[who2011]])', f(hsd, 0), f(hsd * 0.07, 0)],
    ['Excess HSD in adults, Japanese curve (ref. [[morinaga2025]]), disability weight 0.010', f(DJ, 0), f(DJ * 0.010, 0)],
    ['Deaths attributable to near-road NO2, age 30 years and over, relative risk 1.02 per 10 µg m⁻³ (ref. [[huangfu2020]]), 9.9 years per death', '6,460', '63,800'],
    ['Same, years of life lost summed by age group in each small area', '6,460', '64,300'],
    ['Same, relative risk 1.04 per 10 µg m⁻³ (ref. [[boogaard2022]])', 'about 12,700', 'about 125,000']])

txt = '\n'.join(out)
open(os.path.join(HERE, 'si_tables.md'), 'w').write(txt)
for k, v in CHECK.items(): print(k, '::', v)
