# Pipeline

The steps below are the order used for the paper. Run the scripts from the root of this repository (or set `F_ROOT` to the working directory). The scripts read and write under `data/`.

## Directories and environment variables

| Name | Meaning |
|---|---|
| `F_ROOT` | Working directory that holds `src/`, `data/` and `census2020_small_area/`. Default: current directory |
| `F_RECEPTORS` | Path of `receptors_no2_noise_2015_2023.parquet` (https://doi.org/10.5281/zenodo.24299015) |
| `F_NO2_SRC` | Path of `src/` of japan-traffic-no2 (for step 6) |
| `F_BUILDING_CELLS` | Path of `_building_cells.parquet`, the list of building-land cells made by `popalloc_building.py` of japan-traffic-no2 (for step 7) |
| `W` (argument, "noise data dir") | Data directory of japan-night-road-noise (made by its `get_data.sh`): `geo/traffic_census_2021_converted.parquet`, `census27.csv`, `national_validation_points.csv`, `estat_all/shp/r2ka*.shp` |
| `MAT` (argument, "M dir") | Output directory of step 2 |

## Random seeds (all fixed)

| Use | Seed | Script |
|---|---|---|
| Placement of receptors (all prefectures in order 01 to 47) | 20261106 | `src/noise_scen/noise_matrix.py` |
| Placement on building land | 20261123 with the prefecture number | `src/noise_scen/noise_matrix_bldg.py` |
| Monte Carlo simulation, 1,000 runs (two sets of 500; draw −1 is the Japanese curve, draw 0 the central values) | 20261122, 20261123 | `src/mc/mc_alloc_v5.py` |
| Random allocation, 100 orders | 1000–1099 | `src/targeted/v10_review.py rand` |
| Robustness of the selected sections (common errors, 20 exposures) | 606 | `src/targeted/greedy_v6.py` |
| Independent errors between sections | 707 | `src/targeted/greedy_v6_secnoise.py` |
| Bootstrap of municipalities (2,000 samples) | 20261123 | `src/targeted/v10_review.py muni` |

The selection of sections and the assignment of quintiles are deterministic.

## 1. Small-area indicators (2020 Population Census)

Download the small-area tabulation (Tables 3, 6-1, 7-1, 8-1, 9 and 12, one CSV per prefecture) from e-Stat into `census2020_small_area/{t03,t06_1,t07_1,t08_1,t09,t12}/`, and the 2022 deaths by age (statInfId 000040067346) into `data/`.

```
python src/build_adi.py        # data/adi_small_area_2020.csv (renter share p_rent, areal deprivation index ADI_raw)
python src/build_age.py        # data/age_small_area_2020.csv (share aged 20 and over, expected death rate)
python src/build_age2.py       # data/age2_small_area_2020.csv, data/age2_national.csv
```

## 2. Noise at receptors and the energy matrix

```
python src/noise_scen/noise_scenarios.py W src/noise_scen NAT        # receptor-level Lnight with uniform scenarios (NAT/rec_XX.parquet)
python src/heat_expo_receptor.py                                     # tropical nights at receptors (run in the working directory of japan-hotnight-sleep; read by burden_v2.py, results on tropical nights are not reported in the paper)
python src/burden_v2.py NAT                                          # data/receptor/receptor_burden_v2.parquet
python src/targeted/alloc_prep.py W                                  # data/alloc/section_scenarios.parquet (sections, lengths, renter-share-only and random selections)
python src/noise_scen/noise_matrix.py W src/noise_scen MAT           # MAT/M_XX.npz (receptors x sections), kc_XX.parquet, section_fl.npy
python src/targeted/road_length_by_sa.py W data/gen_road_length_by_sa.csv
```

## 3. Health burden and disparity (Figs. 1 and 2)

```
python src/burden_v3.py                      # data/v3_hsd_scenarios.csv, data/v3_by_quintile.csv
python src/fig2/fig2_v8.py                   # data/fig2_v8_by_quintile.csv, data/fig2_v8_gaps.csv
python src/figs/fig1_grid.py                 # data/fig1_grid_*.npz
python src/housing/housing_interaction_v2.py W   # data/v8_housing_regression.csv (Supplementary Table 9)
```

## 4. Allocation to sections (Figs. 3 and 4)

```
python src/targeted/greedy_alloc.py MAT MCREC        # data/v5_alloc_opt.csv, data/alloc/v5_opt_orders.parquet, MCREC/mcrec_XX.parquet
python src/figs/fig3_curves.py MAT W                 # data/fig3_curves.csv
python src/targeted/v9_quintiles.py MAT W main       # data/v9_quintile_rates.csv, v9_cap65.csv, v9_floor40_shares.csv (Supplementary Tables 1 and 6)
python src/targeted/v10_review.py MAT W curve 0 34   # random orders for Fig. 3; also "34 67" and "67 100" (data/v10_random_curves_*.csv)
python src/freeze/repro_check.py MAT W check         # reruns the main selection and compares it with the stored orders
```

The rows of the three `v10_random_curves_*.csv` files (methods `random0` to `random99`) were put into `data/fig3_curves.csv` in place of the random rows written by `fig3_curves.py`.

## 5. Uncertainty (Monte Carlo simulation)

```
python src/mc/mc_alloc_v5.py MCREC 500 20261122 a    # data/v5_mc_draws_a.csv
python src/mc/mc_alloc_v5.py MCREC 500 20261123 b    # data/v5_mc_draws_b.csv
python src/mc/mc_v5_summary.py                       # data/v5_mc_summary.csv
python src/mc/mc_variance_share.py                   # data/v8_mc_variance_share.csv
```

## 6. Heavy vehicles and NO2 co-benefit

```
python src/no2/no2_heavy.py W NO2H                   # NO2H/no2h_XX.parquet
python src/no2/no2_heavy_summary.py NO2H             # data/v5_no2_heavy.csv (Supplementary Table 7)
```

## 7. Sensitivity analyses

```
python src/targeted/greedy_v6.py MAT W               # budget unit, model dependence of the selected sections (data/v6_greedy_checks.csv)
python src/targeted/greedy_v6_local.py MAT W         # finer steps and local search (data/v6_greedy_local.csv)
python src/targeted/greedy_v6_secnoise.py MAT W      # independent errors between sections (data/v6_greedy_secnoise.csv)
python src/targeted/greedy_v7.py MAT W               # rolling noise only, other definitions of equity, continuous treatment units (data/v7_checks.csv)
python src/targeted/greedy_v8_asj.py MAT W           # ASJ RTN-Model 2018 porous asphalt (data/v8_asj_checks.csv)
python src/targeted/greedy_v8b_asj_years.py MAT W    # 10-year average (data/v8b_asj_years.csv)
python src/targeted/v10_units_denominator.py MAT W   # data/v10_units_denominator.csv
python src/targeted/v10_review.py MAT W muni         # within-municipality comparison (data/v10_within_muni_gap.csv)
python src/targeted/v10_review.py MAT W road         # correction by road type (data/v10_roadtype_bias.csv)
python src/targeted/v10_review.py MAT W rentq        # correction by renter-share quintile (data/v10_rentq_bias.csv)
python src/targeted/v10_review.py MAT W popden       # simple roadside-population rule (data/v10_popdensity_rule.csv)
python src/targeted/v10_review.py MAT W rand         # 100 random orders (data/v10_random100.csv)
python src/noise_scen/noise_matrix_bldg.py W src/noise_scen MATB     # receptors placed on building land
python src/targeted/v9_quintiles.py MATB W bldg      # data/v9_bldg_placement.csv (Supplementary Table 10)
python src/validation/cv_residual_strata.py W/national_validation_points.csv TENBOU W/estat_all/shp   # residuals by stratum (Supplementary Table 11)
```

`TENBOU` is the folder with the extracted roadside noise monitoring files of the Environmental Observatory (fiscal years 2019–2023).

## 8. Figures, tables and release files

```
python src/paper/make_figures_en.py                  # figures/Fig1-4
python src/paper/build_si_tables.py                  # Supplementary Tables 1-17
python src/release/build_release_data.py OUT         # small_area_burden_2020.csv, section_allocation.csv of the data deposit
```

## Reproducibility checks made for the paper

- The matrix was generated again and was identical in all 47 prefectures.
- The main selection (uniform 3 dB, 10% of road length) was run again: γ = 0 gave 16,496 sections (reduction 22.09%, disparity ratio 1.634) and γ = 16 gave 16,475 sections (21.49%, 1.521), identical to the stored orders. A stable sort gave the same selection (no ties at the boundary).
- The Monte Carlo runs with seed 20261122 were run again for draws −1 to 3 and matched the stored values.
- `src/paper/make_figures_en.py`, `src/paper/build_si_tables.py` and `src/mc/mc_v5_summary.py`, run from a fresh copy of this repository, reproduce the figures, the supplementary tables and the summary file.

Scripts of earlier, exploratory versions whose outputs are not used in the paper are not included.
